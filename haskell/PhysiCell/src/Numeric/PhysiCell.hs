{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE OverloadedStrings   #-}

module Numeric.PhysiCell where

import qualified Language.C.Inline.Cpp as C

data UserData;

C.context $ C.cppCtx `mappend` C.cppTypePairs [
  ("UserData"  , [t|UserData|])
  ]

C.include "<iostream>"
C.include "<vector>"
C.include "<array>"
C.include "<tuple>"
C.include "<stdexcept>"
C.include "test.h"
C.include "<math.h>"
C.include "<stdio.h>"
C.include "Main.h"

C.include "../../../../core/PhysiCell.h"
C.include "../../../../modules/PhysiCell_standard_modules.h"

C.include "../../../../sample_projects/heterogeneity/custom_modules/heterogeneity.h"

-- instance Storable UserData where
--   poke _ _    = error "poke"
--   peek _      = error "peek"
--   sizeOf _    = error "sizeOf"
--   alignment _ = error "alignment"

runPhysiCell :: IO ()
runPhysiCell = do
  print "Hello world"
