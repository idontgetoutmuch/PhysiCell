{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE OverloadedStrings   #-}

module Main where

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

C.include "../../../core/PhysiCell.h"
C.include "../../../modules/PhysiCell_standard_modules.h"

C.include "heterogeneity.h"

-- instance Storable UserData where
--   poke _ _    = error "poke"
--   peek _      = error "peek"
--   sizeOf _    = error "sizeOf"
--   alignment _ = error "alignment"


main :: IO ()
main = do
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
          double mechanics_voxel_size = 30;
          Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
          create_cell_types();
          setup_tissue();

          set_save_biofvm_mesh_as_matlab( true );
          set_save_biofvm_data_as_matlab( true );
          set_save_biofvm_cell_data( true );
          set_save_biofvm_cell_data_as_custom_matlab( true );

          char filename[1024];
          sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
          save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

          microenvironment.simulate_diffusion_decay( diffusion_dt );
        } |]

  print "Tests finished"
