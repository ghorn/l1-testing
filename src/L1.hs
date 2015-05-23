{-# OPTIONS_GHC -Wall #-}

module L1 where

import Numeric.LinearAlgebra.HMatrix

{-

Basic equations

xdot = Am x + b (mu + theta^T x + sigma0)
y = c^T x

mu = F u


etahat = omegahat u + thetahat^T x + sigmahat

xhatdot = Am xhat + b etahat
yhat = c^T xhat

omegahatdot = Gamma Proj(omegahat, -xtilde^T P b u)
thetahatdot = Gamma Proj(thetahat, -xtilde^T P b x)
sigmahatdot = Gamma Proj(sigmahat, -xtilde^T P b)

u = -k D (etahat - kg r)

-}

rscale = flip scale

{-

The Proj(theta, err) operator.

Pages 291-294 describe this projection operator.  For a long time I
was confused by the "Gamma Proj(theta, err)" notation, and I
thought the Gamma belonged inside the projection operator.  It
turns out it's not quite just hard-bounding the parameter
estimates; it does keep them in the valid range, but this operator
needs to be smooth for the Lyapunov proofs.

-}
fproj :: a -> a -> q a -> a
fproj etheta thetamax theta =
  ((etheta + 1) * theta `dot` theta - maxsq) / (etheta * maxsq)
  where
    maxsq = thetamax `dot` thetamax

gradfproj :: a -> a -> q a -> q a
gradfproj etheta thetamax theta =
  2 * (etheta + 1) / (etheta + maxsq) `rscale` theta
  where
    maxsq = thetamax `dot` thetamax

proj :: a -> a -> q a -> x a -> x a
proj etheta thetamax theta signal
  | ft >= 0 && dfty > 0 = signal - df * dfty `rscale` ft
  | otherwise           = signal
  where
    ft = fproj etheta thetamax theta
    df = gradfproj etheta thetamax theta
    dfty = df `dot` signal

etheta0 :: Fractional a => a
etheta0 = 0.1

{-

Low-pass filter

-}
dstep :: w a -> u a -> y a -> y a
dstep w u y = ydot
  where
    ydot = w * (u - y)

{-

Discrete L1 controller step

-}

type L1States v a = (v a, a, v a, v a, v a)

l1step :: a -> a -> m a -> a -> a -> m a -> v a -> a -> v a -> L1States v a -> L1States v a 
l1step gamma kg p w dt am b r x (xhat, u, omegahat, thetahat, sigmahat) =
  (xhatdot, udot, omegahatdot, thetahatdot, sigmahatdot)
  where
    -- Compute error between reference model and true state
    xtilde = xhat - x
    -- Update parameter estimates.  The estimate derivatives we
    -- compute here will be used to form the next step's estimates; we
    -- use the values we receive as arguments for everything at this
    -- step.
    xtpb = (-xtilde) `dot` p #> b
    gp th sig = gamma `scale` proj etheta0 th sig
    omegahatdot = gp omegahat (xtpb `scale` u)
    thetahatdot = gp thetahat (xtpb `scale` x)
    sigmahatdot = gp sigmahat xtpb
    -- Update reference model state using the previous values.  The
    -- 'xhat' value we receive should be the model's prediction (using
    -- the previous xhat and xhatdot) for the true state 'x' at this
    -- timestep.
    eta = omegahat * u + thetahat `dot` x + sigmahat
    xhatdot = am #> xhat + b `rscale` eta
    -- Update the reference LPF state
    e = kg * r - eta
    udot = dstep w e u
