{-# OPTIONS_GHC -Wall #-}

import qualified Numeric.LinearAlgebra.HMatrix as H

import L1.DT

main = do
  H.saveMatrix "butts" "%.9f" soln
  H.saveMatrix "times" "%.6f" $ H.asRow times'
  H.saveMatrix "reference" "%.9f" $ H.asRow . H.vjoin $ map reference times
  where
    times = [0, 5e-4 .. 15]
    times' = H.fromList times

    -- soln = integrate (ddt oldReference ddtOldSystem 2 ddtControl oldParams)
    --        (H.vjoin [oldSystemState0, oldControlState0]) times'

    soln = integrate (ddt reference ddtSystem 3 ddtControl ctrlParams) state0 times'
