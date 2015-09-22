{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes #-}

module L1.DT
       ( DTState(..), StateSpace(..), DTMIMOParams(..), DisturbanceParams(..), FDParams(..)
       , takeStep, stepControl
       ) where

import Debug.Trace
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.HMatrix ((#>), (<>), (!)) -- (===), (|||)
import Numeric.LinearAlgebra.Data (Matrix, Vector)

import Foreign.Storable (Storable)
--import Data.List (zipWith4)

import qualified Numeric.GSL.ODE as ODE

{- Data types -}
data DTState a =
  DTState
  { stateXHat :: Vector a
  , stateU :: Vector a
  , stateD :: Vector a
  , stateSigmaHat :: Vector a
  , stateEta2 :: Vector a
  } deriving Show

data StateSpace a = StateSpace {
    ssA :: Matrix a
  , ssB :: Matrix a
  , ssC :: Matrix a
  , ssD :: Matrix a
  } deriving (Show, Read)

data DTMIMOParams a =
  DTMIMOParams
  { dmpKg :: Matrix a
--  , dmpSigmaBounds :: (Vector a, Vector a)
  , dmpK :: Matrix a
  , dmpRef :: StateSpace a
  , dmpBum :: Matrix a
  , dmpD :: StateSpace a
  , dmpP :: Matrix a
  , dmpKsp :: Matrix a
  , dmpNStates :: Int
  , dmpNInputs :: Int
  , dmpOmega0 :: Matrix a
  } deriving Show

{- Helper functions -}

--(<#) :: H.Numeric a => Vector a -> Matrix a -> Vector a
--v <# m = H.tr m #> v -- todo(this should be a function)
--
--calcNStates :: (H.Element a, H.Container Vector a, Num a) => DTMIMOParams a -> Int
--calcNStates DTMIMOParams{..} =
--  dmpNStates +           -- xhat
--  dmpNInputs +           -- u
--  dsize dmpD +  -- D(s)
--  dmpNInputs^(2::Int) +  -- omegaHat
--  dmpNInputs*dmpNStates + -- thetaHat
--  dmpNInputs             -- sigmaHat
--  where
--    dsize = fst . H.size . ssA

stepControl :: forall a.
               ( H.Numeric a, Storable a, Show a
               , H.Field a, Num a, Ord a
               , Floating a, Num (Vector a)) =>
               a -> DTMIMOParams a -> Vector a -> DTState a -> Vector a -> DTState a
stepControl ts DTMIMOParams{..} r state xestimate = nextState
  where
    -- Marshalling and unmarshalling
    DTState
      { stateXHat = xhat
      , stateU = u
      , stateD = d
      , stateSigmaHat = sigmaHat
      } = state
    
--    [xhat, u, d, sigmahat] =
--      H.takesV [ dmpNStates
--               , dmpNInputs
--               , fst . H.size . ssA $ dmpD
--               , dmpNStates
--               ] state

    [sigmamhat, sigmaumhat] =
      H.takesV [ dmpNInputs
               , dmpNStates - dmpNInputs
               ] sigmaHat

    nextState =
      DTState
      { stateXHat = nextXHat
      , stateU = nextU
      , stateD = nextD
      , stateSigmaHat = nextSigmaHat
      , stateEta2 = stateEta2 state -- ignoring this
      }
    -- nextState = -- H.vjoin [xhat', u', d', sigmahat']

    -- Feedback controller
    xtilde = xhat - xestimate

    -- State predictor (Euler step)
    xhatdot = ssA dmpRef #> xhat
              + ssB dmpRef #> (dmpOmega0 #> u + sigmamhat)
              + dmpBum #> sigmaumhat
              - dmpKsp #> xtilde
    nextXHat = xhat + ts `H.scale` xhatdot

    -- Parameter estimation
    expmAdt = H.expm $ ts `H.scale` ssA dmpRef
    b = H.fromBlocks [[ssB dmpRef, dmpBum]]
    nextSigmaHat = lins lhs rhs
      where
        lhs = (expmAdt - H.ident dmpNStates) <> b
        rhs = - (ssA dmpRef #> expmAdt #> xtilde)

    -- Matrix-vector linear solve
    lins a = H.flatten . H.linearSolveSVD a . H.reshape 1

    -- Reference filter
    etahat2m = lins (ssB dmpRef) (dmpBum #> sigmaumhat)
    --etahat2m =
    --nextEtaHat2 = etahat2m -- zero of the right dimensions

    eta = dmpOmega0 #> u + sigmamhat + etahat2m + dmpKg #> r
    ddot = ssA dmpD #> d + ssB dmpD #> dmpK #> eta
    nextD = d + ts `H.scale` ddot

    nextU = dmpK #> ssC dmpD #> d


takeStep
  :: DTMIMOParams Double
  -> StateSpace Double
  -> DisturbanceParams Double
  -> (Double -> Vector Double)
  -> Double
  -> (Double, DTState Double, Vector Double)
  -> (Double, DTState Double, Vector Double)
takeStep cp refSys distParams ref dt (t, dtState0, x0) = (t + dt, dtStateNext, xnext)
  where
    u0 = stateU dtState0
    dtStateNext = stepControl dt cp (ref t) dtState0 xestimate
    xestimate = H.subVector 0 3 x0 -- for now

    xnext :: Vector Double
    xnext = last (H.toRows xs)

    xs :: Matrix Double
    xs = int (ddtSystem refSys distParams u0) x0 (H.fromList [t, t + dt])
      where
        int = ODE.odeSolveV ODE.RKck 1e-4 1e-7 1e-7

{- Here is something like the MIMO system tested in section 3.2. -}
ddtSystem
  :: forall a
     . (Show a, Fractional a, Floating a, Num (Vector a), H.Numeric a, H.Indexable (Vector a) a) =>
     StateSpace a -> DisturbanceParams a -> Vector a -> a -> Vector a -> Vector a
ddtSystem refSys distParams u t state = state'
  where
    
    [x, xz, zuv] = H.takesV [3, 2, 2] state

    state' = H.vjoin [x', xz', zuv']

    -- This is the modeled dynamics.
    x' = (ssA refSys + dpA distParams) #> x
         + ssB refSys #> (dpOmega distParams #> u)
         + (f (dpFParams distParams)) x z

    -- This is the unmodeled dynamics
    [xz0, xz1] = H.toList xz
    xz' = H.fromList [xz1, -xz0 + 0.8 * (1 - xz0*xz0) * xz1]

    -- This is how the modeled state couples into the unmodeled
    -- dynamics
    y = H.fromList [1, -2, 1] `H.dot` x
    ydot = H.fromList [1, -2, 1] `H.dot` x'

    -- This is just a scalar calculation
    z = 0.1 * (xz0 - xz1) + zu

    -- This seems to be a hidden state.  It's a scalar, but it's
    -- second-order, so we store it as 'zuv' ('zu vector') composed of
    -- the quantity of interest and its derivative.
    [zu, zu'] = H.toList zuv
    zu'' = y - ydot - 8 * zu' - zu
    zuv' = H.fromList [zu', zu'']
           
    f FDParams{..} xx zz =
      H.fromList
      [ fdpk1 / 3 * (xx `H.dot` xx) + tanh (fdpk2 / 2) * (xx ! 0) + fdpk3 * zz
      , fdpk4 / 2 * (xx ! 1) / (cosh (xx ! 1))
        + fdpk5 / 5 * (xx ! 2)**2
        + fdpk6 / 2 * (1 - exp (-fdplambda * t))
        + fdpk7 / 2 * zz
      , fdpk8 * (xx ! 2) * cos (fdpwu * t) + fdpk9 * zz**2
      ]
       

data DisturbanceParams a =
  DisturbanceParams
  { dpA :: Matrix a
  , dpFParams :: FDParams a
  , dpOmega :: Matrix a
  } deriving Show

data FDParams a =
  FDParams
  { fdpk1 :: a
  , fdpk2 :: a
  , fdpk3 :: a
  , fdpk4 :: a
  , fdpk5 :: a
  , fdpk6 :: a
  , fdpk7 :: a
  , fdpk8 :: a
  , fdpk9 :: a
  , fdplambda :: a
  , fdpwu :: a
  } deriving Show

--{- Here is the old SISO system from chapter 2.2. -}
--
--ddtOldSystem
--  :: (Floating r, Num (Vector r), H.Indexable (Vector r) r,
--      H.Numeric r) =>
--     r -> Vector r -> Vector r -> Vector r
--ddtOldSystem t u x = xdot
--  where
--    -- Dynamics
--    p = x H.! 0
--    v = x H.! 1
--    xdot = H.vjoin [H.scalar v, accel]
--    accel = omega #> (u + H.scalar (m*g*armLen*(cos p) / 2) + sigma + theta #> x)
--
--    -- Disturbances
--    omega = H.scalar $ 1.8 + 0.5*sin(3*t)
--    theta = H.fromLists [[ 0.3 + sin (0.7*pi*t) + cos (3.4 * pi * t)
--                         , -4.8 + 0.4 * sin (2.8*pi*t)]]
--    sigma = H.scalar $
--            3 * (cos p) + 5 * sin (0.85 * pi * t) + cos (7.7*pi/5*t) +
--            2 * cos (5.1 * pi / 3.3 * t) + 3 * sin (1.4 * pi / 1.3 * t) - 4
--
--    -- Parameters
--    m = 1
--    g = 9.8
--    armLen = 0.5
--
--oldReference :: (Floating r, H.Container c r) => r -> c r
--oldReference t = H.scalar $ 2 + cos (2 * t / pi + pi)
--
--oldParams =
--  MIMOParams
--  { mpETheta0 = 0.1
--  , mpOmegaBounds = fromBounds (H.fromLists [[0.1]]) (H.fromLists [[2]])
--  , mpSigmaBounds = fromBounds (H.scalar (-50)) (H.scalar 50)
--  , mpThetaBounds = (H.fromLists [[0, 0]], 5)
--  , mpGamma = 100e3
--  , mpKg = H.scalar 1
--  , mpK = d'omega
--  , mpRef = StateSpace {
--      ssA = H.fromLists [[0, 1], [-1, -1.4]]
--    , ssB = H.fromLists [[0], [1]]
--    , ssC = H.fromLists [[1, 0]]
--    , ssD = H.scalar 0
--    }
--  , mpD = StateSpace {
--      ssA = H.fromLists [[0, 1], [-(d'omega**2), (-2*d'omega*d'zeta)]]
--    , ssB = H.fromLists [[0], [d'omega**2]]
--    , ssC = H.fromLists [[1, 0]]
--    , ssD = H.scalar 0
--    }
--  , mpP = H.fromLists [[1.41429, 0.50000], [0.50000, 0.7142]]
--  , mpKsp = H.fromLists [[0, 0], [0, 1]]
--  , mpNStates = nstate
--  , mpNInputs = ninput
--  }
--  where
--    nstate = 2
--    ninput = 1
--    d'omega = 60
--    d'zeta = 1
--
--oldSystemState0 = H.fromList [-0.2, -0.3]
--
--oldControlState0 = H.vjoin [xhat0, u0, d0, wqsHat0]
--  where
--    xhat0 = H.fromList [0.2, 0.1]
--    u0 = H.scalar 0
--    d0 = H.konst 0 2
--    wqsHat0 = H.vjoin [omegaHat0, thetaHat0, sigmaHat0]
--    omegaHat0 = H.scalar 0.2
--    thetaHat0 = H.konst (-0.5) 2
--    sigmaHat0 = H.scalar (-3.3)
