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

import System.IO.Unsafe ( unsafePerformIO )

import qualified Numeric.LinearAlgebra.HMatrix as HMat
import Linear

--import Casadi.DMatrix ( DMatrix )
import Casadi.MX ( MX )
import Dyno.Vectorize
import Dyno.View.JV
import Dyno.View.M ( fromHMat )
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

l1params :: L1Params (JV RoboX) MX
l1params =
  L1Params
  { l1pETheta0 = 0.1
  , l1pOmegaBounds = fromBounds 0.1 2
  , l1pSigmaBounds = fromBounds (-50) 50
  , l1pThetaBounds = (0, 5)
  , l1pGamma = 100e3
  , l1pKg = 1
  , l1pK = 60
  , l1pP = fromHMat $
           HMat.fromLists
--           [ [0.812564928052004, -0.359997120023040]
--           , [-0.359997120023040, 0.308568960019748]
--           ]
           [ [1.41429,   0.50000]
           , [ 0.50000,   0.7142]
           ]
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
  lol <- prepareL1 ddtRoboX l1params
  let x0 :: RoboX Double
      x0 = RoboX 0 0

      l0 :: L1States RoboX Double
      l0 =
        L1States
        { l1sXhat = x0
        , l1sU = 0
        , l1sWqsHat = wqs0
        }
      reference = 0.0

      wqs0 :: WQS RoboX Double
      wqs0 = fill 1
--      wqs0 = WQS { wqsOmega = 1
--                 , wqsTheta = RoboX 0 0
--                 , wqsSigma = 0
--                 }

  let dfdt :: WQS RoboX Double -> Double -> SimStates RoboX Double -> SimStates RoboX Double
      dfdt wqs r (SimStates x l1) = unsafePerformIO $ do
        let fss =
              FullSystemState
              { ffsX = x
              , ffsWQS = wqs
              }
        (l1', _x') <- lol fss l1 r
        let x' = ddtRoboX fss (l1sU l1)
        return (SimStates x' l1')

--  print $ dfdt wqs0 reference (SimStates x0 l0)

--      simTimes = [0,0.01..2]
      simTimes = [0,0.001..25]
      sols :: [SimStates RoboX Double]
      sols = integrate' (dfdt wqs0 reference) 0.01 simTimes (SimStates x0 l0)
--  mapM_ print sols
  putStrLn $ unlines $ toMatlab "ret" sols
  putStrLn $ "time = " ++ show simTimes ++ ";"
  return ()


toMatlab :: (Vectorize f, Lookup (f Double)) => String -> [f Double] -> [String]
toMatlab topName xs = map (uncurry (fieldToMatlab xs)) at
  where
    at = flatten $ accessors (fill 0)
    fieldToMatlab xzus name get = topName ++ "." ++ name ++ " = " ++ show (map get xzus) ++ ";"
