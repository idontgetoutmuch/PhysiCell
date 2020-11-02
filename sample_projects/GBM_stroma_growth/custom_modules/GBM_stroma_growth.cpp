/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./GBM_stroma_growth.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>

Cell_Definition cancer_cell; 
Cell_Definition stroma_cell; 

// create stroma cell definition
void create_stroma_cells( void )
{
	stroma_cell = cell_defaults;
	
	stroma_cell.name = "stroma cell";
	stroma_cell.type = 2;
	
	// turn off proliferation 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live); 
	stroma_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0; 

	// cell morphology
	stroma_cell.phenotype.geometry.radius = parameters.doubles("stroma_radius");//7.5;
	stroma_cell.phenotype.volume.total = 1767;
	stroma_cell.phenotype.volume.fluid_fraction = 0.75;
	stroma_cell.phenotype.volume.fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.total;
	stroma_cell.phenotype.volume.solid = stroma_cell.phenotype.volume.total-stroma_cell.phenotype.volume.fluid;
	stroma_cell.phenotype.volume.nuclear = 500;
	stroma_cell.phenotype.volume.nuclear_solid = 125;
	stroma_cell.phenotype.volume.nuclear_fluid = stroma_cell.phenotype.volume.nuclear - stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.cytoplasmic = stroma_cell.phenotype.volume.total - stroma_cell.phenotype.volume.nuclear;
	stroma_cell.phenotype.volume.cytoplasmic_fluid = stroma_cell.phenotype.volume.fluid_fraction*stroma_cell.phenotype.volume.cytoplasmic;
	stroma_cell.phenotype.volume.cytoplasmic_solid = stroma_cell.phenotype.volume.cytoplasmic-stroma_cell.phenotype.volume.cytoplasmic_fluid;
	stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 2.53;
	stroma_cell.phenotype.volume.target_solid_cytoplasmic = stroma_cell.phenotype.volume.cytoplasmic_solid;
	stroma_cell.phenotype.volume.target_solid_nuclear = stroma_cell.phenotype.volume.nuclear_solid;
	stroma_cell.phenotype.volume.target_fluid_fraction = stroma_cell.phenotype.volume.fluid_fraction;
	stroma_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = stroma_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	
	// update phenotype function (this is empty as they currently don't do anything)
	stroma_cell.functions.update_phenotype = stroma_function;
	
	return;
}

// create GBM tumour cell definition
void create_cancer_cells( void )
{
	
	cancer_cell = cell_defaults;
		
	// cell morphology
	cancer_cell.phenotype.geometry.radius = 10.75;
	cancer_cell.phenotype.volume.total = 4/3*3.1416*(cancer_cell.phenotype.geometry.radius)*(cancer_cell.phenotype.geometry.radius)*(cancer_cell.phenotype.geometry.radius);//5203.7;
	cancer_cell.phenotype.volume.fluid_fraction = 0.75;
	cancer_cell.phenotype.volume.fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.total;
	cancer_cell.phenotype.volume.solid = cancer_cell.phenotype.volume.total-cancer_cell.phenotype.volume.fluid;
	cancer_cell.phenotype.volume.nuclear = 740;
	cancer_cell.phenotype.volume.nuclear_solid = 185;
	cancer_cell.phenotype.volume.nuclear_fluid = cancer_cell.phenotype.volume.nuclear - cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.cytoplasmic = cancer_cell.phenotype.volume.total - cancer_cell.phenotype.volume.nuclear;
	cancer_cell.phenotype.volume.cytoplasmic_fluid = cancer_cell.phenotype.volume.fluid_fraction*cancer_cell.phenotype.volume.cytoplasmic;
	cancer_cell.phenotype.volume.cytoplasmic_solid = cancer_cell.phenotype.volume.cytoplasmic-cancer_cell.phenotype.volume.cytoplasmic_fluid;
	cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic/cancer_cell.phenotype.volume.nuclear;//6.0321;
	cancer_cell.phenotype.volume.target_solid_cytoplasmic = cancer_cell.phenotype.volume.cytoplasmic_solid;
	cancer_cell.phenotype.volume.target_solid_nuclear = cancer_cell.phenotype.volume.nuclear_solid;
	cancer_cell.phenotype.volume.target_fluid_fraction = cancer_cell.phenotype.volume.fluid_fraction;
	cancer_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = cancer_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio;
	cancer_cell.phenotype.volume.calcified_fraction = 0; 
	cancer_cell.phenotype.volume.calcification_rate = 0;
	cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.35*cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength;
	
	//update phenoypte function (updates proliferation rate of each cell and uptakes their movement)
	cancer_cell.functions.update_phenotype = cancer_cell_proliferation_movement;
	
	cancer_cell.name = "cancer cell";
	cancer_cell.type = 1; 
	
	return;
}

// create the two cell type definitions by creating a general definition
void create_cell_types( void )
{
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	//cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;//0.00073549;//0.0000064812*(1-9740/512720);
	
 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 

	// initialise cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.0;//0.05;
		
	// add variables to track stop or go GBM cell persistence time
	cell_defaults.custom_data.add_variable( "persistence_time", "dimensionless", 0.0 ); // how long cells will persist in move or stop phenotype
	cell_defaults.custom_data.add_variable( "cell_motility_type", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "rep_rate", "dimensionless", 0.0 );
				
	// turn off secretion from these cells (oxygen and wall which is the migratory domain)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	static int wall_index = microenvironment.find_density_index( "wall");
	
	cell_defaults.phenotype.secretion.secretion_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[wall_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[wall_index] = 10; 
		
	// update cell and phenotype
	cell_defaults.functions.update_phenotype = cancer_cell_proliferation_movement;
	cell_defaults.phenotype.motility.is_motile = true; 
		
	cell_defaults.name = "holder cell"; 
	cell_defaults.type = 0; 
	
	//create cell types
	create_cancer_cells();
	create_stroma_cells();
		
	return; 
}

void setup_microenvironment( void )
{
	// make sure not override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	initialize_microenvironment(); 	

	// initialise oxygen (which is a redundant variable for now) and wall (which represents the migratory and non migratory domain)
	static int oxygen_index = microenvironment.find_density_index( "oxygen");
	static int wall_index = microenvironment.find_density_index( "wall");
	
	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		
		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")+10)*(parameters.doubles("tumour_radius")+10))//1280*1280)//420*420)//
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[wall_index] = 1;
			
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("tumour_radius")-10)*(parameters.doubles("tumour_radius")-10))//420*420)//1260*1260)//1260*1260)400*400)//
		{
			
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[wall_index] = 3.5;
		}
		else
		{	
			microenvironment(n)[oxygen_index] = 0;
			microenvironment(n)[wall_index] = 5;
		}	
	}
	return; 
}

// this function initialises cells in the domain
void setup_tissue_circle( void )
{
	double Radius = parameters.doubles("tumour_radius");
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	// setting up distributions for movement and persistance of cells
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
			
	double GBM_NO = parameters.ints("initial_GBM_cells");
	double stroma_NO = parameters.ints("initial_stroma_cells");
			
	// place and initialise GBM cells 	
	for( int i=0; i<GBM_NO; i++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		
		pCell = create_cell( cancer_cell );		
		pCell->assign_position( x , y , 0.0 );
		
		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
				
		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();
			
			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pCell->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		} 
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP
	
		double go_stop_var = UniformRandom();
		
		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}
			
		}
	}
	
	// place and initialise stroma cells
	for( int l=0; l<stroma_NO; l++ )
	{	
		double R = sqrt(UniformRandom())*Radius;
		double alp = UniformRandom()*2*3.141;
		x = R*cos(alp);
		y = R*sin(alp);
		pCell = create_cell( stroma_cell );
		pCell->assign_position( x , y , 0.0 );
	}
	
	return;
}

// This function determines the proliferation rate of GBM cells and calls movement function
void cancer_cell_proliferation_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	double max_pressure = parameters.doubles("maximum_proliferation_pressure");
	
	//Determines the proliferation rate of each cell
	if(pCell->type ==1)
	{	
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
		
		if( pCell->state.simple_pressure<=max_pressure)
		{
			pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("GBM_cell_proliferation_rate");
		}
		else if( pCell->state.simple_pressure>max_pressure) 
		{
			pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;
		}
	}
	
	// calls GBM cell movement function
	cell_movement( pCell, phenotype, dt);
	return;
	
}

// determines the movement of each GBM cell (i.e. stop or go and it's speed)
void cell_movement( Cell* pCell, Phenotype& phenotype, double dt ) 
{

	static int wall_index = microenvironment.find_density_index( "wall" ); 
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	
	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
	
	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; // 1 = go, 2 = stop
	
	// discretized distribution for cell 'go' persistence times
	std::vector<double> go_times_cumul(8);
    go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;
	
	// discretized distribution for cell persistence times
	std::vector<double> persistence_times_vec(8);
    persistence_times_vec[0] = 0;
	persistence_times_vec[1] = 30;
	persistence_times_vec[2] = 60;
	persistence_times_vec[3] = 90;
	persistence_times_vec[4] = 120;
	persistence_times_vec[5] = 150;
	persistence_times_vec[6] = 180;
	persistence_times_vec[7] = 240;
	
	// discretized distribution for cell speed
	std::vector<double> speed_cumul(12);
    speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;
	
	std::vector<double> speed_vec(12);
    speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;
	
	if( wall_amount<2 & pCell->type == 1)
	{
		pCell->phenotype.motility.migration_speed = 0;
	}
	else		
	{
		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
		{
			// assign new type (stop = 2, or go = 1)
			double new_type_rand = UniformRandom();
			if(new_type_rand<=0.5)// GO
			{
				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type
				
				double speed_var = UniformRandom();
			
				for( int k=0; k<12; )
				{
					if( speed_var> speed_cumul[k] )
					{k++;}
					else
					{
						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
						k = 12;
					}
				}
			} 
			else
			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type

			// assign persistence time - needs to be a real time!
			double go_stop_var = UniformRandom();
			for( int j=0; j<8; )
			{
				if( go_stop_var> go_times_cumul[j] )
				{j++;}
					else
				{
					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
					j = 8;
				}
			}
		}
			
	}
	return;
}

void stroma_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	return;
}

// this function plots the cells
std::vector<std::string> colouring( Cell* pCell )
{
	std::vector< std::string > output( 4, "darkgrey" ); 
		
	if( pCell->type == 1 && pCell->phenotype.death.dead==false)
	{
		output[0] = "rgb(104, 55, 99)";//"orchid";//"rgb(255,230,230)";
		output[1] = "rgb(104, 55, 99)";
		output[2] = "rgb(85, 50, 70)";//"plum";//"rgb(255,230,230)";
		output[3] = "rgb(85, 50, 70)";

		return output; 
					
	}
	else if( pCell->phenotype.death.dead == true )
	{ 
		output[0] = "rgb(255, 255, 224)";
		output[1] = "rgb(255, 255, 224)";
		output[2] = "rgb(255, 228, 181)";
		output[3] = "rgb(255, 228, 181)";
			
		return output; 
	}
	else if( pCell->type == 2)
	{
		output[0] = "rgb(234, 172, 199)";//"rgb(255,230,230)";
		output[1] = "rgb(234, 172, 199)";
		output[2] = "rgb(243, 186, 211)";//"rgb(255,230,230)";
		output[3] = "rgb(243, 186, 211)";
		
		return output; 
	}
	return output;
}
