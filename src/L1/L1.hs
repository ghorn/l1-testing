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

import qualified Numeric.GSL.ODE as ODE
import qualified Numeric.LinearAlgebra.Data as D
import qualified Data.Vector as V
import qualified Data.Vector.Storable as SV
import qualified Numeric.LinearAlgebra.HMatrix as HMat

import Casadi.CMatrix ( CMatrix )
import Casadi.SX ( SX )
import Casadi.DMatrix ( DMatrix )
import Casadi.MX ( MX )
import Casadi.Overloading ( SymOrd(..) )

import Accessors

import Dyno.Vectorize
import Dyno.View.M ( M, ms, mm, trans, uncol, col, hsplitTup, fromHMat )
import Dyno.View.Viewable
import Dyno.View.JV
import Dyno.View.Fun
import Dyno.View.FunJac
import Dyno.View.View

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
    maxsq = thetamax * thetamax

gradfproj :: (View x, Viewable a, CMatrix a) => S a -> S a -> M x (JV Id) a -> M x (JV Id) a
gradfproj etheta thetamax theta =
  (2 * (etheta + 1) / (etheta*maxsq)) `scale` theta
  where
    maxsq = thetamax * thetamax

gt :: (Num a, SymOrd a) => a -> a -> a
gt x y = 1 - (x `leq` y)

-- todo(move this and its friend to SymOrd (casadi bindings))
--lt :: (Num a, SymOrd a) => a -> a -> a
--lt x y = 1 - (x `geq` y)

proj :: forall x a
        . (View x, Viewable a, CMatrix a, SymOrd (M x (JV Id) a))
        => S a -> S a -> M x (JV Id) a -> M x (JV Id) a -> M x (JV Id) a
proj etheta thetamax theta y =
  y - correction `scale` (ft `scale` (ndfty `scale` ndf))
  where
    correction :: S a
    correction = (ft `geq` 0)*(dfty `gt` 0)
    ft :: S a
    ft = fproj etheta thetamax theta
    df :: M x (JV Id) a
    df = gradfproj etheta thetamax theta
    ndf :: M x (JV Id) a
    ndf = (1 / sqrt (df `dot` df)) `scale` df
    ndfty :: S a
    ndfty = ndf `dot` y
    dfty :: S a
    dfty = df `dot` y


--etheta0 :: Fractional a => a
--etheta0 = 0.1

{-

Low-pass filter

-}

scale :: (Viewable a, CMatrix a, View f, View g) => S a -> M f g a -> M f g a
scale s m = m `ms` (uncol s)

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
  { controllerState :: J (JV (L1States x)) a
  , systemState :: J (JV (FullSystemState x)) a
  } deriving (Generic, Generic1)
instance Vectorize x => View (FullL1State x)

data L1Params x a =
  L1Params
  { l1pETheta0 :: S a
  , l1pOmegaBounds :: (S a, S a)
  , l1pSigmaBounds :: (S a, S a)
  , l1pThetaBounds :: (M x (JV Id) a, S a)
  , l1pGamma :: S a
  , l1pKg :: S a
  , l1pK :: S a
  , l1pP :: M x x a
  , l1pW :: S a
  } deriving Functor


type S a = M (JV Id) (JV Id) a

prepareL1 ::
  forall x u
  . (Vectorize x, u ~ Id)
  => (FullSystemState x (J (JV Id) SX) -> J (JV Id) SX -> x (J (JV Id) SX))
  -> L1Params (JV x) MX -- todo(symbolic leak)
  -> IO (FullSystemState x Double -> L1States x Double -> Double -> IO (L1States x Double))
prepareL1 userOde l1params = do
  let f :: JacIn (JTuple (JV x) (JV u)) (J (JV (WQS x))) SX -> JacOut (JV x) (J JNone) SX
      f (JacIn xu wqs) = JacOut (catJV' xdot) (cat JNone)
        where
          xdot :: x (J (JV Id) SX)
          xdot = userOde ffs u

          ffs :: FullSystemState x (J (JV Id) SX)
          ffs =
            FullSystemState
            { ffsX = splitJV' x
            , ffsWQS = splitJV' wqs
            }

          JTuple x u = split xu

  sxF <- toSXFun "user ode" f
  userJac <- toFunJac sxF

  let fullOde :: JTuple (FullL1State x) (JV Id) MX -> J (JV (L1States x)) MX
      fullOde (JTuple fullL1States r) = l1dot
        where
          l1dot = ddtL1States l1params dx'dx dx'du (col r) (col xestimate) (catJV' controllerState')

          FullL1State controllerState'' systemState'' = split fullL1States
          controllerState'@(L1States xhat u wqsHat) = splitJV' controllerState''
          systemState' = splitJV' systemState''

          xestimate = catJV' (ffsX systemState')
          
          --jacIn :: JacIn (JTuple f0 (JV Id)) (WQS x) (J (JV Id) MX)
          jacIn = JacIn (cat (JTuple (catJV' xhat) (catJV' (Id u))))
                        (catJV' wqsHat)
          dx'dxu :: M (JV x) (JTuple (JV x) (JV u)) MX
          Jac dx'dxu _ _ = call userJac jacIn

          dx'dx :: M (JV x) (JV x) MX
          dx'du :: M (JV x) (JV u) MX
          dx'dx = fromHMat $ HMat.fromLists [[0, 1], [-1, -1.4]]
          dx'du = fromHMat $ HMat.fromLists [[0], [1]]
          --(dx'dx, dx'du) = hsplitTup dx'dxu

  fullOdeMX <- toMXFun "full ode with l1" fullOde

  let retFun :: FullSystemState x Double -> L1States x Double -> Double
                -> IO (L1States x Double)
      retFun ffs l1States r = do
        let fullL1States :: FullL1State x DMatrix
            fullL1States =
              FullL1State
              { controllerState = v2d $ catJV l1States
              , systemState = v2d $ catJV ffs
              }
            input = JTuple (cat fullL1States) (v2d (catJV (Id r)))
        ret <- eval fullOdeMX input
        return (splitJV (d2v ret))

  return retFun


ddtL1States ::
  forall a x
  . (Vectorize x, Viewable a, CMatrix a)
  => L1Params (JV x) a
  -> M (JV x) (JV x) a -- am
  -> M (JV x) (JV Id) a -- b
  -> S a
  -> M (JV x) (JV Id) a
  -> J (JV (L1States x)) a -> J (JV (L1States x)) a
ddtL1States L1Params{..} am b r xestimate l1states =
  catJV' $ L1States (splitJV' (uncol xhatdot)) (unId (splitJV' (uncol udot))) wqsDot
  where
    L1States xhat0 u0 wqsHat' = splitJV' l1states

    xhat :: M (JV x) (JV Id) a
    xhat = col (catJV' xhat0)

    u :: M (JV Id) (JV Id) a
    u = col u0

    wqsDot :: WQS x (J (JV Id) a)
    wqsDot =
      WQS
      { wqsOmega = unId $ splitJV' (uncol omegahatdot)
      , wqsTheta =        splitJV' (uncol thetahatdot)
      , wqsSigma = unId $ splitJV' (uncol sigmahatdot)
      }

    wqsHat :: WQS x (J (JV Id) a)
    wqsHat = wqsHat'
    omegahat :: M (JV Id) (JV Id) a
    omegahat = col (catJV' (Id (wqsOmega wqsHat)))
    thetahat :: M (JV x) (JV Id) a
    thetahat = col (catJV' (wqsTheta wqsHat))
    sigmahat :: M (JV Id) (JV Id) a
    sigmahat = col (catJV' (Id (wqsSigma wqsHat)))

    -- Compute error between reference model and true state
    xtilde :: M (JV x) (JV Id) a
    xtilde = xhat - xestimate
    -- Update parameter estimates.  The estimate derivatives we
    -- compute here will be used to form the next step's estimates; we
    -- use the values we receive as arguments for everything at this
    -- step.
    xtpbg :: S a
    xtpbg = l1pGamma `scale` ((-(trans xtilde)) `mm` l1pP `mm` b)

    gp :: View f => (M f (JV Id) a, S a)-> M f (JV Id) a -> M f (JV Id) a -> M f (JV Id) a
    gp (scenter, snorm) th sig = unshift ret0
      where
        ret0 = proj l1pETheta0 snorm (shift th) (shift sig)

        shift z   = z - scenter
        unshift z = z + scenter


    omegahatdot,sigmahatdot :: S a
    omegahatdot = gp l1pOmegaBounds omegahat (xtpbg `scale` u)
    sigmahatdot = gp l1pSigmaBounds sigmahat xtpbg
    thetahatdot :: M (JV x) (JV Id) a
    thetahatdot = gp l1pThetaBounds thetahat (xtpbg `scale` xestimate)
    -- Update reference model state using the previous values.  The
    -- 'xhat' value we receive should be the model's prediction (using
    -- the previous xhat and xhatdot) for the true state 'x' at this
    -- timestep.

    eta :: S a
    eta = omegahat * u + thetahat `dot` xestimate + sigmahat

    xhatdot :: M (JV x) (JV Id) a
    xhatdot = am `mm` xhat + eta `scale` b
    -- Update the reference LPF state
    e :: S a
    e = l1pKg * r - eta

    udot :: S a
    udot = l1pW * (l1pK * e - u)

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
