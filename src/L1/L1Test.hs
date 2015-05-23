{-# OPTIONS_GHC -Wall #-}
{-# Language ScopedTypeVariables #-}
{-# Language DeriveGeneric #-}
{-# Language DeriveFunctor #-}

module L1.L1Test
       ( RoboX(..)
       , RoboU(..)
       , ddtRoboX
       ) where

import GHC.Generics ( Generic, Generic1 )

import Linear

import Dyno.Vectorize

data RoboX a =
  RoboX
  { xPos :: a
  , xVel :: a
  } deriving (Functor, Generic, Generic1)
instance Vectorize RoboX

data RoboU a =
  RoboU
  { uTorque :: a
  } deriving (Functor, Generic, Generic1)
instance Vectorize RoboU

--data RoboTheta a =
--  RoboTheta
--  { qF1 :: a
--  , qF2 :: a
--  } deriving (Functor, Generic, Generic1)
--instance Vectorize RoboTheta

ddtRoboX :: Floating a => RoboX a -> RoboU a -> RoboX a -> a -> a -> RoboX a
ddtRoboX x@(RoboX p v) (RoboU u) theta sigmaBar omegaBar =
  RoboX
  { xPos = v
  , xVel = (u + m*g*r*(cos p)/2 + sigmaBar + x `dot` theta) * omegaBar
  }
  where
    m = 1
    g = 9.8
    r = 0.5
