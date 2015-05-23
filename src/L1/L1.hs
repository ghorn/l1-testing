{-# OPTIONS_GHC -Wall #-}
{-# Language ScopedTypeVariables #-}
{-# Language FlexibleContexts #-}

module L1 where

import Casadi.CMatrix ( CMatrix )
import Casadi.Overloading ( SymOrd(..) )

import Dyno.Vectorize
import Dyno.View.M
import Dyno.View.Viewable
import Dyno.View.JV
import Dyno.View.View hiding ( J )

dot :: (View f, CMatrix a) => M f (JV Id) a -> M f (JV Id) a -> M (JV Id) (JV Id) a
dot x y = trans x `mm` y

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

p-}

--rscale :: Container c e
--rscale = flip scale

{-

The Proj(theta, err) operator.

Pages 291-294 describe this projection operator.  For a long time I
was confused by the "Gamma Proj(theta, err)" notation, and I
thought the Gamma belonged inside the projection operator.  It
turns out it's not quite just hard-bounding the parameter
estimates; it does keep them in the valid range, but this operator
needs to be smooth for the Lyapunov proofs.

-}
fproj :: (View x, Viewable a, CMatrix a) => S a -> S a -> M x (JV Id) a -> S a
fproj etheta thetamax theta =
  ((etheta + 1) * theta `dot` theta - maxsq) / (etheta * maxsq)
  where
    maxsq = thetamax `dot` thetamax

gradfproj :: (View x, Viewable a, CMatrix a) => S a -> S a -> M x (JV Id) a -> M x (JV Id) a
gradfproj etheta thetamax theta =
  (2 * (etheta + 1) / (etheta + maxsq)) `scale` theta
  where
    maxsq = thetamax `dot` thetamax

gt :: (Num a, SymOrd a) => a -> a -> a
gt x y = 1 - (x `leq` y)

lt :: (Num a, SymOrd a) => a -> a -> a
lt x y = 1 - (x `geq` y)

proj :: forall x a
        . (View x, Viewable a, CMatrix a, SymOrd (M x (JV Id) a))
        => S a -> S a -> M x (JV Id) a -> M x (JV Id) a -> M x (JV Id) a
proj etheta thetamax theta signal =
  signal - (1 - aOk) `scale` (ft `scale` (dfty `scale` df))
  where
    aOk :: S a
    aOk = (ft `geq` 0)*(dfty `gt` 0)
    ft :: S a
    ft = fproj etheta thetamax theta
    df :: M x (JV Id) a
    df = gradfproj etheta thetamax theta
    dfty :: S a
    dfty = df `dot` signal


--etheta0 :: Fractional a => a
--etheta0 = 0.1

{-

Low-pass filter

-}
dstep :: (CMatrix a, Viewable a) => S a -> S a -> S a -> S a
dstep w u y = ydot
  where
    ydot = w * (u - y)

scale :: (Viewable a, CMatrix a, View f, View g) => S a -> M f g a -> M f g a
scale s m = m `ms` (uncol s)

{-

Discrete L1 controller step

-}
data L1States x a =
  L1States
  { l1sXhat :: M x (JV Id) a
  , l1sU :: S a
  , l1sOmega :: S a -- M v (JV Id) a
  , l1sTheta :: M x (JV Id) a
  , l1sSigma :: S a -- M s (JV Id) a
  }

type S a = M (JV Id) (JV Id) a

l1step ::
  forall a x
  . (View x, Viewable a, CMatrix a)
  => S a
  -> S a
  -> S a
  -> S a
  -> S a
  -> S a
  -> M x x a -- p
  -> S a
  -> M x x a -- am
  -> M x (JV Id) a -- b
  -> S a
  -> M x (JV Id) a
  -> L1States x a -> L1States x a
l1step etheta0 omegamax sigmamax thetamax gamma kg p w am b r x
  (L1States xhat u omegahat thetahat sigmahat) =
    L1States xhatdot udot omegahatdot thetahatdot sigmahatdot
  where
    -- Compute error between reference model and true state
    xtilde :: M x (JV Id) a
    xtilde = xhat - x
    -- Update parameter estimates.  The estimate derivatives we
    -- compute here will be used to form the next step's estimates; we
    -- use the values we receive as arguments for everything at this
    -- step.
    xtpb :: S a
    xtpb = (-(trans xtilde)) `mm` p `mm` b

    gp :: View f => S a -> M f (JV Id) a -> M f (JV Id) a -> M f (JV Id) a
    gp somethingMax th sig = gamma `scale` proj etheta0 somethingMax th sig

    omegahatdot,sigmahatdot :: S a
    omegahatdot = gp omegamax omegahat (xtpb `scale` u)
    sigmahatdot = gp sigmamax sigmahat xtpb
    thetahatdot :: M x (JV Id) a
    thetahatdot = gp thetamax thetahat (xtpb `scale` x)
    -- Update reference model state using the previous values.  The
    -- 'xhat' value we receive should be the model's prediction (using
    -- the previous xhat and xhatdot) for the true state 'x' at this
    -- timestep.

    eta :: S a
    eta = omegahat * u + thetahat `dot` x + sigmahat

    xhatdot :: M x (JV Id) a
    xhatdot = am `mm` xhat + eta `scale` b
    -- Update the reference LPF state
    e :: S a
    e = kg * r - eta

    udot :: S a
    udot = dstep w e u
