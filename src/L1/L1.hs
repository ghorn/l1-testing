{-# OPTIONS_GHC -Wall #-}
{-# Language TypeFamilies #-}
{-# Language ScopedTypeVariables #-}
{-# Language FlexibleContexts #-}
{-# Language DeriveGeneric #-}
{-# Language DeriveFunctor #-}
{-# Language RecordWildCards #-}

module L1.L1
       ( L1States(..)
       , FullSystemState(..)
       , L1Params(..)
       , WQS(..)
       , prepareL1
       , integrate
       , integrate'
       ) where

import GHC.Generics ( Generic, Generic1 )

import qualified Data.Foldable as F
import Linear
import qualified Numeric.GSL.ODE as ODE
import qualified Numeric.LinearAlgebra.Data as D
import qualified Data.Vector as V
import qualified Data.Vector.Storable as SV

import Casadi.Overloading ( SymOrd(..) )
import Accessors
import Dyno.Vectorize

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
fproj :: (Metric x, Fractional a) => a -> a -> x a -> a
fproj etheta thetamax theta =
  ((etheta + 1) * theta `dot` theta - maxsq) / (etheta * maxsq)
  where
    maxsq = thetamax * thetamax

gradfproj :: (Functor x, Fractional a) => a -> a -> x a -> x a
gradfproj etheta thetamax theta =
  (2 * (etheta + 1) / (etheta*maxsq)) *^ theta
  where
    maxsq = thetamax * thetamax

gt :: (Num a, SymOrd a) => a -> a -> a
gt x y = 1 - (x `leq` y)

-- todo(move this and its friend to SymOrd (casadi bindings))
--lt :: (Num a, SymOrd a) => a -> a -> a
--lt x y = 1 - (x `geq` y)

proj :: forall x a
        . (Metric x, Floating a, SymOrd a) => a -> a -> x a -> x a -> x a
proj etheta thetamax theta y =
  y ^-^ correction *^ (ft *^ (ndfty *^ ndf))
  where
    correction :: a
    correction = (ft `geq` 0)*(dfty `gt` 0)
    ft :: a
    ft = fproj etheta thetamax theta
    df :: x a
    df = gradfproj etheta thetamax theta
    ndf :: x a
    ndf = (1 / sqrt (df `dot` df)) *^ df
    ndfty :: a
    ndfty = ndf `dot` y
    dfty :: a
    dfty = df `dot` y


--etheta0 :: Fractional a => a
--etheta0 = 0.1

{-

Low-pass filter

-}

{-

Discrete L1 controller step

-}


data L1States x a =
  L1States
  { l1sXhat :: x a
  , l1sU :: a
  , l1sWqsHat :: WQS x a
  } deriving (Functor, Generic, Generic1, Show)
instance Vectorize x => Vectorize (L1States x)
instance (Lookup a, Lookup (x a)) => Lookup (L1States x a)

data WQS x a =
  WQS
  { wqsOmega :: a
  , wqsTheta :: x a
  , wqsSigma :: a
  } deriving (Functor, Generic, Generic1, Show)
instance Vectorize x => Vectorize (WQS x)
instance (Lookup a, Lookup (x a)) => Lookup (WQS x a)

data FullSystemState x a =
  FullSystemState
  { ffsX :: x a
  , ffsWQS :: WQS x a
  } deriving (Functor, Generic, Generic1)
instance Vectorize x => Vectorize (FullSystemState x)

data FullL1State x a =
  FullL1State
  { controllerState :: L1States x a
  , systemState :: FullSystemState x a
  } deriving (Functor, Generic, Generic1)
instance Vectorize x => Vectorize (FullL1State x)

data L1Params x a =
  L1Params
  { l1pETheta0 :: a
  , l1pOmegaBounds :: (a, a)
  , l1pSigmaBounds :: (a, a)
  , l1pThetaBounds :: (x a, a)
  , l1pGamma :: a
  , l1pKg :: a
  , l1pK :: a
  , l1pP :: x (x a)
  , l1pW :: a
  , l1pKsp :: x (x a)
  } deriving Functor


prepareL1 ::
  forall x u a
  . (Metric x, F.Foldable x, u ~ Id, Floating a, SymOrd a)
  => L1Params x a
  -> x (x a) -- dx/dx
  -> x a     -- dx/du
  -> FullSystemState x a -> L1States x a -> a -> L1States x a
prepareL1 l1params dxdx dxdu = retFun
  where
    retFun :: FullSystemState x a -> L1States x a -> a -> L1States x a
    retFun ffs l1States r =
      fullOde FullL1State { controllerState = l1States, systemState = ffs } r

    fullOde :: FullL1State x a -> a -> L1States x a
    fullOde (FullL1State controllerState systemState) r =
      ddtL1States l1params dxdx dxdu r (ffsX systemState) controllerState


ddtL1States ::
  forall a x
  . (Additive x, F.Foldable x, Metric x, Floating a, SymOrd a)
  => L1Params x a
  -> x (x a) -- am
  -> x a -- b
  -> a
  -> x a
  -> L1States x a -> L1States x a
ddtL1States L1Params{..} am b r xestimate l1states =
  L1States xhatdot udot wqsDot
  where
    wqsHat :: WQS x a
    L1States xhat u wqsHat = l1states

    wqsDot :: WQS x a
    wqsDot =
      WQS
      { wqsOmega = omegahatdot
      , wqsTheta = thetahatdot
      , wqsSigma = sigmahatdot
      }

    omegahat :: a
    omegahat = wqsOmega wqsHat
    thetahat :: x a
    thetahat = wqsTheta wqsHat
    sigmahat :: a
    sigmahat = wqsSigma wqsHat

    -- Compute error between reference model and true state
    xtilde :: x a
    xtilde = xhat ^-^ xestimate
    -- Update parameter estimates.  The estimate derivatives we
    -- compute here will be used to form the next step's estimates; we
    -- use the values we receive as arguments for everything at this
    -- step.
    xtpbg :: a
    xtpbg = - l1pGamma * (xtilde `dot` (l1pP !* b))

    gp :: (Additive f, Metric f) => (f a, a)-> f a -> f a -> f a
    gp (scenter, snorm) th sig = unshift ret0
      where
        ret0 = proj l1pETheta0 snorm (shift th) (shift sig)

        shift z   = z ^-^ scenter
        unshift z = z ^+^ scenter

    -- same as gp but for scalars, wrap with Id and call gp
    gp' :: (a, a)-> a -> a -> a
    gp' (scenter, snorm) th sig = unId $ gp (Id scenter, snorm) (Id th) (Id sig)


    omegahatdot,sigmahatdot :: a
    omegahatdot = gp' l1pOmegaBounds omegahat (xtpbg * u)
    sigmahatdot = gp' l1pSigmaBounds sigmahat xtpbg
    thetahatdot :: x a
    thetahatdot = gp l1pThetaBounds thetahat (xtpbg *^ xestimate)
    -- Update reference model state using the previous values.  The
    -- 'xhat' value we receive should be the model's prediction (using
    -- the previous xhat and xhatdot) for the true state 'x' at this
    -- timestep.

    eta :: a
    eta = omegahat * u + thetahat `dot` xestimate + sigmahat

    xhatdot :: x a
    xhatdot = am !* xhat ^+^ eta *^ b ^-^ l1pKsp !* xtilde
    -- Update the reference LPF state
    e :: a
    e = l1pKg * r - eta

    udot :: a
    udot = l1pK * e
--    udot = l1pW * (l1pK * e - u)

integrate :: Vectorize x => (x Double -> x Double) -> Double -> x Double -> x Double
integrate f h x0 = devectorize $ sv $ last sol
  where
    vs :: V.Vector Double -> SV.Vector Double
    vs = SV.fromList .  V.toList
    sv :: SV.Vector Double -> V.Vector Double
    sv =  V.fromList . SV.toList

    sol = D.toRows $
          ODE.odeSolveV
          ODE.MSAdams
          h 1e-7 1e-5 f'
          (vs (vectorize x0))
          (SV.fromList [0, h])
    f' :: Double -> SV.Vector Double -> SV.Vector Double
    f' _ x = vs $ vectorize $ f (devectorize (sv x))


integrate' :: Vectorize x
              => (Double -> x Double -> x Double)
              -> Double -> [Double] -> x Double -> [x Double]
integrate' f h times x0 = map (devectorize . sv) sol
  where
    vs :: V.Vector Double -> SV.Vector Double
    vs = SV.fromList .  V.toList
    sv :: SV.Vector Double -> V.Vector Double
    sv =  V.fromList . SV.toList

    sol = D.toRows $
          ODE.odeSolveV
          ODE.MSAdams
          h 1e-7 1e-5 f'
          (vs (vectorize x0))
          (SV.fromList times)
    f' :: Double -> SV.Vector Double -> SV.Vector Double
    f' t x = vs $ vectorize $ f t (devectorize (sv x))
