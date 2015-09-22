{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards #-}

module Main ( main ) where

import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.Data ( Vector )

import L1.DT

dSys :: StateSpace Double
dSys =
  StateSpace
  { ssA = H.fromLists
          [ [ 0.00000,  0.00000, 0.00000, 0.00000, 0.00000 ]
          , [-0.01000,  0.00000, 0.00000, 0.00000, 28.00000]
          , [ 0.00000, 10.00000, 0.00000, 0.00000, 27.80000]
          , [ 0.00000,  0.00000, 100.000, 0.00000, 101.90000]
          , [ 0.00000,  0.00000, 0.00000, -100.00, -167.00000]
          ]
  , ssB = H.fromLists [ [-0.014142, -0.14142]
                      , [0.000000, 0.000000]
                      , [0.000000, 0.000000]
                      , [0.000000, 0.000000]
                      , [0.000000, 0.000000]
                      ]

  , ssC = H.fromLists [ [ 0.00000, 0.00000, 0.00000, 0.00000, -0.070711 ]
                      , [ 0.00000, 0.00000, 0.00000, 0.00000, -0.070711 ]
                      ]
  , ssD = H.fromLists [ [ 0, 0 ]
                      , [ 0, 0 ]
                      ]
  }

-- Case 1 from the textbook
distParams0 :: DisturbanceParams Double
distParams0 =
  DisturbanceParams
  { dpA = H.fromLists [ [0.2, -0.2, -0.3]
                      , [-0.2, -0.2, 0.6]
                      , [-0.1, 0, -0.9]
                      ]
  , dpOmega = H.fromLists [ [0.6, 0.2]
                          , [0.2, 1.2]]
  , dpFParams = FDParams (-1) 1 0 1 0 0.2 1 0.6 (-0.7) 0.3 5
  }
              
refSys :: StateSpace Double
refSys = StateSpace am b c d
  where
  am = H.fromLists [ [-1,  0,  0]
                   ,  [0,  0,  1]
                   ,  [0, -1, -1.8]
                   ]
  b = H.fromLists [[1, 0], [0, 0], [1, 1]]
  c = H.fromLists [[1, 0, 0], [0, 1, 0]]
  d = H.konst 0 (2, 2)

ctrlState0 :: DTState Double
ctrlState0 =
  DTState
  { stateXHat = H.konst 0 3
  , stateU = H.konst 0 2
  , stateD = H.konst 0 5
  , stateSigmaHat = H.konst 0 3
  , stateEta2 = H.scalar 0 -- error "FUCK OFF"
  }

reference :: Double -> Vector Double
reference _ = H.fromList [-0.8, 0.8]

systemState0 :: Vector Double
systemState0 = H.fromList [0, 0, 1, 0, 0, 0, 0]

ctrlParams :: DTMIMOParams Double
ctrlParams =
  DTMIMOParams
  { dmpKg = H.fromLists [[1, 0], [-1, 1]]
  , dmpBum = H.nullspace $ H.tr' (ssB (dmpRef ctrlParams))
  , dmpK = H.diagl [8, 8]
  , dmpRef = refSys
  , dmpD = dSys
  , dmpP = H.fromLists [[0.5, 0, 0], [0, 1.45556, -0.5], [0, -0.5, 0.55556]]
  , dmpKsp = H.konst 0 (nstate, nstate)
  , dmpNStates = nstate
  , dmpNInputs = ninput
  , dmpOmega0 = H.ident 2
  }
  where
    nstate = 3
    ninput = 2

toCol DTState{..} = H.vjoin [stateXHat, stateU, stateD, stateSigmaHat]

main :: IO ()
main = do
  H.saveMatrix "butts" "%.9f" (H.fromRows (map third soln))
  H.saveMatrix "times" "%.6f" $ H.asRow times'
  H.saveMatrix "reference" "%.9f" $ H.fromColumns $ map reference times
  H.saveMatrix "ctrl" "%.9f" . H.fromColumns . map (toCol . second) $ soln
  where
    dt = 1e-2
    tf = 50
    count = length times
    times = [0, dt .. tf]
    times' = H.fromList times

    -- soln = integrate (ddt oldReference ddtOldSystem 2 ddtControl oldParams)
    --        (H.vjoin [oldSystemState0, oldControlState0]) times'

--    soln = integrate ctrlParams refSys distParams0 reference 5e-4 ddtSystem  state0 times'
    third (_,_,x) = x
    second (_, x, _) = x
    soln = take count $ iterate (takeStep ctrlParams refSys distParams0 reference dt)
           (0, ctrlState0, systemState0)
