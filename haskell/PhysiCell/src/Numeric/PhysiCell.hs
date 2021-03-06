{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE OverloadedStrings   #-}

module Numeric.PhysiCell (
  runPhysiCell
  ) where

import qualified Language.C.Inline.Cpp as C
import           Foreign.Ptr (Ptr)
import           Control.Monad

import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM


data Test;
data Microenvironment;

C.context $ C.cppCtx `mappend` C.cppTypePairs [
    ("Test::Test", [t|Test|])
  , ("BioFVM::Microenvironment", [t|Microenvironment|])
  ]

C.include "<iostream>"
C.include "<omp.h>"
C.include "<vector>"
C.include "<array>"
C.include "<tuple>"
C.include "<stdexcept>"
C.include "<math.h>"
C.include "<stdio.h>"

C.include "../../../../core/PhysiCell.h"
C.include "../../../../modules/PhysiCell_standard_modules.h"

-- FIXME
C.include "../cancer_biorobots.h"

-- instance Storable UserData where
--   poke _ _    = error "poke"
--   peek _      = error "peek"
--   sizeOf _    = error "sizeOf"
--   alignment _ = error "alignment"

runPhysiCell :: IO ()
runPhysiCell = do
  flag <- [C.block|  bool {
          static bool XML_status = false;
          XML_status = load_PhysiCell_config_file( "./cancer_biorobots.xml" );
          return XML_status;
        } |]
  print flag

  when (flag == 0) (error "Problem with config file")

  pt <- [C.block| BioFVM::Microenvironment* {
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

          static int oxygen_index  = microenvironment.find_density_index( "oxygen" );
          std::cout << "\nOxygen Index: " << oxygen_index << std::endl;

          static int chemoattractant_index  = microenvironment.find_density_index( "chemoattractant" );
          std::cout << "\nChemoattractant Index: " << chemoattractant_index << std::endl;

	  std::cout << "total agents: " << all_cells->size() << std::endl;

          display_simulation_status( std::cout );

          sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

          save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

          microenvironment.simulate_diffusion_decay( diffusion_dt );

          ((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );

          PhysiCell_globals.current_time += diffusion_dt;

          display_simulation_status( std::cout );

          std::cout << default_microenvironment_options.calculate_gradients << std::endl;

          std::cout << microenvironment.mesh.voxels.size() << std::endl;

          std::cout << "\nx dimension: " << microenvironment.mesh.x_coordinates.size()
                    << "  y dimension: " << microenvironment.mesh.y_coordinates.size()
                    << "  z dimension: " << microenvironment.mesh.z_coordinates.size()
                    << std::endl;

          return &microenvironment;
        } |] :: IO (Ptr Microenvironment)

  xSz <- [C.block| int {
          static int oxygen_index  = (*$(BioFVM::Microenvironment* pt)).find_density_index( "oxygen" );
          std::cout << "\nOxygen Index: " << oxygen_index << std::endl;
          std::cout << (*$(BioFVM::Microenvironment* pt))(1)[oxygen_index] << std::endl;
          std::cout << (*$(BioFVM::Microenvironment* pt)).number_of_voxels() << std::endl;
          static int xSz = microenvironment.mesh.x_coordinates.size();
          return xSz;
        } |]

  print xSz
  print "Tests finished"
