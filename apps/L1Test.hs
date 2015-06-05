{-# OPTIONS_GHC -Wall #-}
{-# Language FlexibleContexts #-}
{-# Language ScopedTypeVariables #-}
{-# Language DeriveGeneric #-}
{-# Language DeriveFunctor #-}

module Main
       ( main
       , RoboX(..)
       , SimStates(..)
       , ddtRoboX
       ) where

import GHC.Generics ( Generic, Generic1 )

import Linear

import Dyno.Vectorize
import Accessors

import L1.L1

data RoboX a =
  RoboX
  { xPos :: a
  , xVel :: a
  } deriving (Functor, Generic, Generic1, Show)
instance Vectorize RoboX
instance Lookup a => Lookup (RoboX a)

ddtRoboX :: Floating a => FullSystemState RoboX a -> a -> RoboX a
ddtRoboX (FullSystemState x@(RoboX p v) (WQS omegaBar theta sigmaBar)) u  =
  RoboX
  { xPos = v
  , xVel = (usat u + m*g*armLen*(cos p)/2 + sigmaBar + x `dot` theta) * omegaBar
  }
  where
    -- This approximation is a bit more linear in the mid-range than a
    -- plain logistic function.
    usat t = saturate 0 (maxTorque*2) 0.3 t - t*4.4 / (10 + (t / sqrt(10))**4)
    m = 1
    g = 9.8
    armLen = 0.5
    maxTorque = 8
    saturate center range slope t = center + range * (-0.5 + 1 / (1 + exp (-slope * t)))

fromBounds :: Fractional a => a -> a -> (a, a)
fromBounds lb ub = (0.5 * (lb + ub), 0.5*(ub - lb))

l1params :: L1Params RoboX RoboX Double
l1params =
  L1Params
  { l1pETheta0 = 0.1
  , l1pOmegaBounds = fromBounds 0.1 2
  , l1pSigmaBounds = fromBounds (-50) 50
  , l1pThetaBounds = (RoboX 0 0, 5)
  , l1pGamma = 100e3
  , l1pKg = 1
  , l1pK = d'omega
  -- A second-order filter for D gives a third-order input filter
  -- overall, when combined with the outer "r - eta" loop.
  , l1pDA = RoboX (RoboX 0               1                  )
                  (RoboX (-(d'omega**2)) (-2*d'omega*d'zeta))
  , l1pDB = RoboX 0 (d'omega**2)
  -- Here is a first-order D instead for comparison.
  --
  -- , l1pDA = RoboX (RoboX (-d'omega) 0) (RoboX 0 0)
  -- , l1pDB = RoboX (d'omega) 0
  , l1pDC = RoboX 1 0
  , l1pP = RoboX (RoboX 1.41429 0.50000)
                 (RoboX 0.50000 0.7142)
  , l1pKsp = RoboX (RoboX 0 0)
                   (RoboX 0 1)
  }
  where
    -- This parameter specifies the bandwidth for both the outer loop
    -- and the inner filter D(s) of the input low-pass filter.
    d'omega = 60
    -- This parameter specifies the damping coefficient for D(s)
    d'zeta = 1

data SimStates d x a =
  SimStates
  { ssX :: x a
  , ssL1 :: L1States d x a
  } deriving (Functor, Generic, Generic1, Show)
instance (Vectorize d, Vectorize x) => Vectorize (SimStates d x)
instance (Lookup a, Lookup (d a), Lookup (x a)) => Lookup (SimStates d x a)

main :: IO ()
main = do
  let reference t = cos (2 * t / pi) + pi
      dxdx = RoboX (RoboX 0 1)
                   (RoboX (-1) (-1.4))
      dxdu = RoboX 0 1
      lol = prepareL1 l1params dxdx dxdu
      x0 :: RoboX Double
      x0 = RoboX (-0.2) (-0.3)

      l0 :: L1States RoboX RoboX Double
      l0 =
        L1States
        { l1sXhat = RoboX (0.2) (0.1)
        , l1sU = 0
        , l1sD = fill 0
        , l1sWqsHat = WQS {
            wqsOmega = 0.2
          , wqsTheta = fill (-0.5)
          , wqsSigma = -3.3
          }
        }

      dfdt :: Double -> SimStates RoboX RoboX Double -> SimStates RoboX RoboX Double
      dfdt t (SimStates x l1) = SimStates x' l1'
        where
          fss =
            FullSystemState
            { ffsX = x
            , ffsWQS = wqs
            }
          wqs =
            WQS
            { wqsOmega = 1.8 + 0.5*sin(3*t)
            -- todo(mp): breaks when it's zero
            , wqsTheta =
              RoboX
              { xPos = 0.3 + sin (0.7*pi*t) + cos (3.4 * pi * t)
              , xVel = -4.8 + 0.4 * sin (2.8*pi*t)
              }
            -- todo(mp): breaks when it's zero
            , wqsSigma = 2 * (cos p) + 5 * sin (0.85 * pi * t) + cos (7.7*pi/5*t) +
                         2 * cos (5.1 * pi / 3.3 * t) + 3 * sin (1.4 * pi / 1.3 * t)
            }
          p = xPos x
          r = reference t
          l1' = lol fss l1 r
          x' = ddtRoboX fss (l1sU l1)

      simTimes = [0,0.0002..25]
      refs = map reference simTimes
      sols :: [SimStates RoboX RoboX Double]
      sols = integrate' dfdt 0.01 simTimes (SimStates x0 l0)
  putStrLn $ unlines $
    toMatlab "ret" sols ++ [" r = [" ++ unwords (map show refs) ++ "];"]
  putStrLn $ "time = " ++ show simTimes ++ ";"
  return ()


toMatlab :: (Vectorize f, Lookup (f Double)) => String -> [f Double] -> [String]
toMatlab topName xs = map (uncurry (fieldToMatlab xs)) at
  where
    at = flatten $ accessors (fill 0)
    fieldToMatlab xzus name get = topName ++ "." ++ name ++ " = " ++ show (map get xzus) ++ ";"
