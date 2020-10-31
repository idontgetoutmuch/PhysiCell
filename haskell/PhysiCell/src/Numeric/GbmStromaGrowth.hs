{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes         #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE OverloadedStrings   #-}

module Numeric.GbmStromaGrowth (
  runGbmStromaGrowth
  ) where

import qualified Language.C.Inline.Cpp as C
import           Foreign.Ptr (Ptr)
import           Control.Monad

import           Data.Int

import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM


data Test
data Microenvironment
data Cells

C.context $ C.cppCtx `mappend` C.cppTypePairs [
    ("BioFVM::Microenvironment", [t|Microenvironment|])
  , ("Test::Cells", [t|Cells|])
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
C.include "./test.h"

C.include "../GBM_stroma_growth.h"


numCancerCells :: Ptr Cells -> IO C.CInt
numCancerCells pt = do
  [C.block| int {
   return (*$(Test::Cells* pt)).size();
 } |]

runGbmStromaGrowth :: IO ()
runGbmStromaGrowth = do

  flag <- [C.block|  bool {
          using namespace PhysiCell;
          static bool XML_status = false;
          XML_status = load_PhysiCell_config_file( "./GBM_stroma_growth.xml" );
          return XML_status;
        } |]
  print flag

  when (flag == 0) (error "Problem with config file")

  pt <- [C.block| BioFVM::Microenvironment* {
          using namespace PhysiCell;

          // OpenMP setup
	  omp_set_num_threads(PhysiCell_settings.omp_num_threads);
          std::cout << "\nThreads: " << PhysiCell_settings.omp_num_threads << std::endl;

          // From PhysiCell but not in GBM
          SeedRandom();

          // time setup
          std::string time_units = "min";

          /* Microenvironment setup */
          setup_microenvironment();
          double last_time_updated = 0; // keep tabs of the last time the radius was increased
          double number_of_times_updated = 1; // the number of times the radius has been increased

          /* PhysiCell setup */

          // set mechanics voxel size, and match the data structure to BioFVM
          double mechanics_voxel_size = 30;

          Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
          create_cell_types();
          setup_tissue_circle();

          set_save_biofvm_mesh_as_matlab( true );
          set_save_biofvm_data_as_matlab( true );
          set_save_biofvm_cell_data( true );
          set_save_biofvm_cell_data_as_custom_matlab( true );

          char filename[1024];
          sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
          save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

          PhysiCell_SVG_options.length_bar = 200;

          // for simplicity, set a pathology coloring function

          std::vector<std::string> (*cell_coloring_function)(Cell*) = colouring;

          sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
          SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

          microenvironment.simulate_diffusion_decay( diffusion_dt );

          static int oxygen_index  = microenvironment.find_density_index( "oxygen" );
          std::cout << "\nOxygen Index: " << oxygen_index << std::endl;

          static int wall_index  = microenvironment.find_density_index( "wall" );
          std::cout << "\nWall Index: " << wall_index << std::endl;

          display_simulation_status( std::cout );

	  std::cout << "total agents: " << all_cells->size() << std::endl;

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

  rt <- [C.block| Test::Cells* {
         return all_cells;
        } |] :: IO (Ptr Cells)

  nCanCells <- numCancerCells rt

  print nCanCells
  print "Tests finished"
