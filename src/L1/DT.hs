{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE RankNTypes #-}

module L1.DT where

import qualified Numeric.LinearAlgebra.HMatrix as H

import Numeric.LinearAlgebra.HMatrix ((#>), (<>), (===), (|||))
import Numeric.LinearAlgebra.Data (Matrix, Vector)

import Foreign.Storable (Storable)
import Data.List (zipWith4)

import qualified Numeric.GSL.ODE as ODE

{- Data types -}

data StateSpace a = StateSpace {
    ssA :: Matrix a
  , ssB :: Matrix a
  , ssC :: Matrix a
  , ssD :: Matrix a
  } deriving (Show, Read)

data DTMIMOParams a = DTMIMOParams {
    dmpSigmaBounds :: (Vector a, Vector a)
  , dmpKg :: Matrix a
  , dmpK :: Matrix a
  , dmpRef :: StateSpace a
  , dmpRefUnmodeled :: StateSpace a
  , dmpD :: StateSpace a
  , dmpP :: Matrix a
  , dmpKsp :: Matrix a
  , dmpNStates :: Int
  , dmpNInputs :: Int
  , dmpOmega0 :: Matrix a
  }

type State a = Vector a

{- Helper functions -}

(<#) :: H.Numeric a => Vector a -> Matrix a -> Vector a
v <# m = H.tr m #> v

calcNStates :: (H.Element a, H.Container Vector a) => MIMOParams a -> Int
calcNStates MIMOParams{..} = mpNStates +           -- xhat
                             mpNInputs +           -- u
                             dsize mpD +  -- D(s)
                             mpNInputs^(2::Int) +  -- omegaHat
                             mpNInputs*mpNStates + -- thetaHat
                             mpNInputs             -- sigmaHat
  where
    dsize = fst . H.size . ssA

{-

Feedback controller -- a direct MIMO extension of the controller from
section 2.2.

-}
stepControl :: forall a.
               ( H.Numeric a, Storable a
               , H.Field a, Num a, Ord a
               , Floating a, Num (Vector a)) =>
               a -> DTMIMOParams a -> Vector a -> State a -> Vector a -> State a
stepControl ts DTMIMOParams{..} r state xestimate = statedot
  where
    -- Marshalling and unmarshalling
    [xhat, u, d, sigmahat] =
      H.takesV [ dmpNStates
               , dmpNInputs
               , fst . H.size . ssA $ dmpD
               , dmpNStates
               ] state
    [sigmamhat, sigmaumhat] =
      H.takesV [ dmpNInputs
               , dmpNStates - dmpNInputs
               ] sigmahat

    statedot = H.vjoin [xhat', u', d', sigmahat']

    -- Feedback controller
    xtilde = xhat - xestimate

    -- State predictor (Euler step)
    xhatdot = ssA dmpRef #> xhat
            + ssB dmpRef #> (dmpOmega0 #> u + sigmamhat)
            + ssB dmpRefUnmodeled #> sigmaumhat
            - dmpKsp #> xtilde
    xhat' = xhat + ts `H.scale` xhatdot

    -- Parameter estimation
    eps = H.expm $ ts `H.scale` ssA dmpRef
    b = H.fromBlocks [[ssB dmpRef, ssB dmpRefUnmodeled]]
    lhs = (eps - H.ident dmpNStates) <> b
    rhs = - (ssA dmpRef #> eps #> xtilde)
    sigmahat' = lins lhs rhs

    -- Matrix-vector linear solve
    lins a = H.flatten . H.linearSolveSVD a . H.reshape 1

    -- Reference filter
    etahat2m = lins (ssB dmpRef) (ssB dmpRefUnmodeled #> sigmaumhat)

    eta = dmpOmega0 #> u + sigmamhat + etahat2m + dmpKg #> r
    ddot = ssA dmpD #> d + ssB dmpD #> eta
    d' = d + ts `H.scale` ddot

    u' = dmpK #> ssC dmpD #> d


{- Integration tools -}

ddt ref ddtS ns ddtC cp t fullstate = fullstatedot
  where
    [systemState, ctrlState] = H.takesV [ns, ctrlSize] fullstate
    fullstatedot = H.vjoin [systemStatedot, ctrlStatedot]

    nstate = mpNStates cp
    ninput = mpNInputs cp
    ctrlSize = calcNStates cp

    -- Time-varying reference
    refvec = ref t

    -- Controller dynamics
    ctrlStatedot = ddtC cp refvec ctrlState $
                   H.subVector 0 ns systemState
    u = H.subVector nstate ninput ctrlState

    -- System dynamics
    systemStatedot = ddtS t u systemState

integrate cp ref dt ctrl sys c0 x0 tf =
  where


integrate' cp ref dt ctrl sys c0 x0 t tf =
  where
    u = H.subVector nstate ninput c0

    nstate = dmpNStates cp
    ninput = dmpNInputs cp
    refvec = ref t

    x = int (ddtSystem u) x0 (H.scalar $ tf - t)
    c = ctrl refvec c0 x0
    int = ODE.odeSolveV ODE.RKck 1e-4 1e-7 1e-7

{- Here is something like the MIMO system tested in section 3.2. -}

ddtSystem u t x = xdot
  where
    omega _ = H.fromLists [[0.28, -0.2], [0.0, 0.25]]
    theta _ = H.fromLists [[0.1, -0.2, 0.04], [2.1, -0.5, 1.4]]
    sigma _ = H.fromList [0.1, -0.1]

    xdot = ssA refSys #> x + ssB refSys #> (omega t #> u + theta t #> x + sigma t)

state0 = H.vjoin [systemState0, ctrlState0]

reference _ = H.konst 0 2

refSys = StateSpace am b c d
  where
  am = H.fromLists [[-1, 0, 0], [0, 0, 1], [0, -1, -1.8]]
  b = H.fromLists [[1, 0], [0, 0], [1, 1]]
  c = H.fromLists [[1, 0, 0], [0, 1, 0]]
  d = H.konst 0 (2, 2)

dSys = StateSpace am b c $ H.konst 0 (2, 2)
  where
    d'omega = 8
    d'zeta = 1.1
    sisoAm = H.fromLists [[0, 1], [-(d'omega**2), (-2*d'omega*d'zeta)]]
    sisoB = H.fromLists [[0], [d'omega**2]]
    sisoC = H.fromLists [[1, 0]]
    am = H.diagBlock $ replicate 2 sisoAm
    b = H.diagBlock $ replicate 2 sisoB
    c = H.diagBlock $ replicate 2 sisoC

ctrlParams = DTMIMOParams {
    dmpKg = H.fromLists [[1, 0], [-1, 1]]
  , dmpK = H.diagl [8, 8]
  , dmpRef = refSys
  , dmpD = dSys
  , dmpP = H.fromLists [[0.5, 0, 0], [0, 1.45556, -0.5], [0, -0.5, 0.55556]]
  , dmpKsp = H.konst 0 (nstate, nstate)
  , dmpNStates = nstate
  , dmpNInputs = ninput
  }
  where
    nstate = 3
    ninput = 2

systemState0 = H.fromList [0, 0, 1]

ctrlState0 = H.vjoin [xhat0, u0, d0, wqsHat0]
  where
    xhat0 = systemState0
    u0 = H.konst 0 2
    d0 = H.konst 0 4
    wqsHat0 = H.vjoin [H.flatten omegahat0, H.flatten thetahat0, sigmahat0]
    omegahat0 = H.fromLists [[0.25, -0.2], [-0.2, 0.25]]
    thetahat0 = H.konst 0 (2, 3)
    sigmahat0 = H.konst 0 2

{- Here is the old SISO system from chapter 2.2. -}

ddtOldSystem t u x = xdot
  where
    -- Dynamics
    p = x H.! 0
    v = x H.! 1
    xdot = H.vjoin [H.scalar v, accel]
    accel = omega #> (u + H.scalar (m*g*armLen*(cos p) / 2) + sigma + theta #> x)

    -- Disturbances
    omega = H.scalar $ 1.8 + 0.5*sin(3*t)
    theta = H.fromLists [[ 0.3 + sin (0.7*pi*t) + cos (3.4 * pi * t)
                         , -4.8 + 0.4 * sin (2.8*pi*t)]]
    sigma = H.scalar $
            3 * (cos p) + 5 * sin (0.85 * pi * t) + cos (7.7*pi/5*t) +
            2 * cos (5.1 * pi / 3.3 * t) + 3 * sin (1.4 * pi / 1.3 * t) - 4

    -- Parameters
    m = 1
    g = 9.8
    armLen = 0.5

oldReference t = H.scalar $ 2 + cos (2 * t / pi + pi)

oldParams = MIMOParams {
    mpETheta0 = 0.1
  , mpOmegaBounds = fromBounds (H.fromLists [[0.1]]) (H.fromLists [[2]])
  , mpSigmaBounds = fromBounds (H.scalar (-50)) (H.scalar 50)
  , mpThetaBounds = (H.fromLists [[0, 0]], 5)
  , mpGamma = 100e3
  , mpKg = H.scalar 1
  , mpK = d'omega
  , mpRef = StateSpace {
      ssA = H.fromLists [[0, 1], [-1, -1.4]]
    , ssB = H.fromLists [[0], [1]]
    , ssC = H.fromLists [[1, 0]]
    , ssD = H.scalar 0
    }
  , mpD = StateSpace {
      ssA = H.fromLists [[0, 1], [-(d'omega**2), (-2*d'omega*d'zeta)]]
    , ssB = H.fromLists [[0], [d'omega**2]]
    , ssC = H.fromLists [[1, 0]]
    , ssD = H.scalar 0
    }
  , mpP = H.fromLists [[1.41429, 0.50000], [0.50000, 0.7142]]
  , mpKsp = H.fromLists [[0, 0], [0, 1]]
  , mpNStates = nstate
  , mpNInputs = ninput
  }
  where
    nstate = 2
    ninput = 1
    d'omega = 60
    d'zeta = 1

oldSystemState0 = H.fromList [-0.2, -0.3]

oldControlState0 = H.vjoin [xhat0, u0, d0, wqsHat0]
  where
    xhat0 = H.fromList [0.2, 0.1]
    u0 = H.scalar 0
    d0 = H.konst 0 2
    wqsHat0 = H.vjoin [omegaHat0, thetaHat0, sigmaHat0]
    omegaHat0 = H.scalar 0.2
    thetaHat0 = H.konst (-0.5) 2
    sigmaHat0 = H.scalar (-3.3)
