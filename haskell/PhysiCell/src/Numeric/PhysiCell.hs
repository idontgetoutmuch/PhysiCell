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
C.include "<omp.h>"
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

C.include "heterogeneity.h"

-- instance Storable UserData where
--   poke _ _    = error "poke"
--   peek _      = error "peek"
--   sizeOf _    = error "sizeOf"
--   alignment _ = error "alignment"


runPhysiCell :: IO ()
runPhysiCell = do
  flag <- [C.block|  bool {
          XML_status = false;
          XML_status = load_PhysiCell_config_file( "./PhysiCell_settings.xml" );
          return XML_status;
        } |]
  print flag

  [C.block| void {
          using namespace PhysiCell;
          std::cout << "\nThreads: " << PhysiCell_settings.omp_num_threads << std::endl;
          SeedRandom();
          std::string time_units = "min";
          setup_microenvironment();
        } |]

  print "Tests finished"
