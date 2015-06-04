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

--data RoboU a =
--  RoboU
--  { uTorque :: a
--  } deriving (Functor, Generic, Generic1)
--instance Vectorize RoboU

--data RoboTheta a =
--  RoboTheta
--  { qF1 :: a
--  , qF2 :: a
--  } deriving (Functor, Generic, Generic1)
--instance Vectorize RoboTheta

ddtRoboX :: Floating a => FullSystemState RoboX a -> a -> RoboX a
ddtRoboX (FullSystemState x@(RoboX p v) (WQS omegaBar theta sigmaBar)) u  =
--ddtRoboX x@(RoboX p v) u theta sigmaBar omegaBar =
  RoboX
  { xPos = v
  , xVel = (u + m*g*armLen*(cos p)/2 + sigmaBar + x `dot` theta) * omegaBar
  }
  where
    m = 1
    g = 9.8
    armLen = 0.5

fromBounds :: Fractional a => a -> a -> (a, a)
fromBounds lb ub = (0.5 * (lb + ub), 0.5*(ub - lb))

l1params :: L1Params RoboX Double
l1params =
  L1Params
  { l1pETheta0 = 0.1
  , l1pOmegaBounds = fromBounds 0.1 2
  , l1pSigmaBounds = fromBounds (-50) 50
  , l1pThetaBounds = (RoboX 0 0, 5)
  , l1pGamma = 100e3
  , l1pKg = 1
  , l1pK = 60
  , l1pP = RoboX (RoboX 1.41429 0.50000)
                 (RoboX 0.50000 0.7142)
  , l1pW = 1
  }

data SimStates x a =
  SimStates
  { ssX :: x a
  , ssL1 :: L1States x a
  } deriving (Functor, Generic, Generic1, Show)
instance Vectorize x => Vectorize (SimStates x)
instance (Lookup a, Lookup (x a)) => Lookup (SimStates x a)

main :: IO ()
main = do
  let dxdx = RoboX (RoboX 0 1)
                   (RoboX (-1) (-1.4))
      dxdu = RoboX 0 1
      lol = prepareL1 l1params dxdx dxdu
      x0 :: RoboX Double
      x0 = RoboX (-0.2) (-0.3)

      l0 :: L1States RoboX Double
      l0 =
        L1States
        { l1sXhat = RoboX (0.4) (0.5)
        , l1sU = 0
        , l1sWqsHat = wqs0
        }
      wqs0 :: WQS RoboX Double
      wqs0 = WQS { wqsOmega = 0.2
                 , wqsTheta = fill (-0.5) -- todo(mp): breaks when it's zero
                 , wqsSigma = -3.3      -- todo(mp): breaks when it's zero
                 }

      --dfdt :: WQS RoboX Double -> Double -> SimStates RoboX Double -> SimStates RoboX Double
      --dfdt wqs r (SimStates x l1) = unsafePerformIO $ do
      dfdt :: Double -> SimStates RoboX Double -> SimStates RoboX Double
      dfdt t (SimStates x l1) = SimStates x' l1'
        where
          fss =
            FullSystemState
            { ffsX = x
            , ffsWQS = wqs
            }
          --wqs = wqs0
          wqs =
            WQS
            { wqsOmega = 1.1 + 0.3*sin(3*t)
            , wqsTheta =
              RoboX
              { xPos = sin (0.5*pi*t) + cos (pi * t)
              , xVel = -1 + 0.1 * sin (3*pi*t)
              }
            , wqsSigma = (cos p) + 2 * sin (pi * t) + cos (7*pi/5*t)
            }
          p = xPos x
          r = cos (2*t/pi)
          -- r = 0
          l1' = lol fss l1 r
          x' = ddtRoboX fss (l1sU l1)

--  print $ dfdt wqs0 reference (SimStates x0 l0)

--      simTimes = [0,0.01..2]
      simTimes = [0,0.0002..10]
      sols :: [SimStates RoboX Double]
      sols = integrate' dfdt 0.01 simTimes (SimStates x0 l0)
--  mapM_ print sols
  putStrLn $ unlines $ toMatlab "ret" sols
  putStrLn $ "time = " ++ show simTimes ++ ";"
  return ()


toMatlab :: (Vectorize f, Lookup (f Double)) => String -> [f Double] -> [String]
toMatlab topName xs = map (uncurry (fieldToMatlab xs)) at
  where
    at = flatten $ accessors (fill 0)
    fieldToMatlab xzus name get = topName ++ "." ++ name ++ " = " ++ show (map get xzus) ++ ";"
