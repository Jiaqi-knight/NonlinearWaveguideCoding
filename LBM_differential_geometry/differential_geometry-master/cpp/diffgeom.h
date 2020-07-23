/*
 * =====================================================================================
 * 
 * Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
 * Author: Luca Di Stasio <luca.distasio@gmail.com>
 *                        <luca.distasio@ingpec.eu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * =====================================================================================
 */

#ifndef DIFFGEOM_H
#define DIFFGEOM_H

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <typeinfo>
#include <type_traits> // C++0x
//#include <tr1/type_traits> // C++03, use std::tr1
#include <vector>

using namespace std;

//============================================================================//
//============================================================================//
/*
      A class to perform differential geometry computation on a given map
*/
//============================================================================//
//============================================================================//

//===================================================
//==================  HEADER  =======================
//===================================================

template<class T = double>
class diffgeom {

//===================================================  
//                  Variables
//===================================================
private:
  
// I/O management
  string date,                                                // Current date and time
         WorkDir,                                             // Working directory (for input and output files, not code)
	       InpDir,                                              // Input files directory
	       OutDir,                                              // Output files directory
	       InpFileName,                                         // Input file name with data and settings to run the simulation
	       OutFileName,                                         // Output filename to save results
	       RunMode,                                             // Mode of operation. Options: gui: start graphical user interface; shell: start interactive shell mode in current console; file: read data from file and start computation, without user interaction. 
	       EndMapManFileName,                                   // CSV file with coordinates of points describing the co-domain manifold
	       BaseMapManFileName;                                  // CSV file with coordinates of points describing the domain manifold
 
// Input parameters
  int Dim,                                                    // Space dimensions 
      Nx,                                                     // Number of points in x-direction
      Ny,                                                     // Number of points in x-direction
      Nz,                                                     // Number of points in x-direction
      A;                                                      // Accuracy of finite difference approximations
	 
// Control parameters
  vector<bool> computed_quantities;                           // Keep track of geometrical quantities computed
  unsigned int float_operations;                              // Number of floating point operations performed

// Output quantities
  vector<vector<double> > comp_nodes;                         // List of nodes' coordinates in computational space, ordered by helical numbering
  vector<vector<double> > phys_nodes;                         // List of nodes' coordinates in physical space, ordered by helical numbering
  vector<vector<double> > cov_base;                           // Covariant base vectors
  vector<vector<double> > contr_base;                         // Contravariant base vectors
  vector<vector<double> > cov_metric;                         // Covariant metric tensor
  vector<vector<double> > contr_metric;                       // Contravariant metric tensor
  vector<vector<double> > first_christoffel;                  // First Christoffel symbols
  vector<vector<double> > second_christoffel;                 // Second Christoffel symbols
  vector<vector<double> > contr_Riemann;                      // Contravariant Riemann curvature tensor
  vector<vector<double> > cov_Riemann;                        // Covariant Riemann curvature tensor
  vector<vector<double> > cov_Ricci;                          // Covariant Ricci curvature tensor
  double Ricci_scalar;                                        // Ricci scalar curvature

  
//===================================================  
//                      Methods
//===================================================  
public:
  
  // Constructor (default)
  diffgeom();
  
  // Constructor when specific run-mode is provided
  diffgeom(string input);
  
  //Destructor
  ~diffgeom();
  
  // Functions
  const string currentDateTime();                                                // Get current system date and time

  void print_intro(),                                                            // Print header to screen
       set_IO(),                                                                 // Set I/O parameters interactively
       get_inpfilename(const string&),                                           // Get the input file name from string variable in main, which depends on command line input
       edit_inpfile(),                                                           // Edit the input file
       read_inpfile(),                                                           // Read the input file
       print_inputdata(),                                                        // Print data read from input file
       initialize(const string&),                                                // Initialize the simulation
       read_coordinates(),                                                       // Read coordinates from .csv files
       get_coordinates(vector<vector<double> >, vector<vector<double> >),        // Get coordinates from API
       compute_covbase(),                                                        // Compute covariant base vectors
       compute_contrbase(),                                                      // Compute contravariant base vectors
       compute_covmetric(),                                                      // Compute covariant metric tensor
       compute_contrmetric(),                                                    // Compute contravariant metric tensor
       compute_firstChristoffel(),                                               // Compute first Christoffel symbols
       compute_secondChristoffel(),                                              // Compute second Christoffel symbols
       compute_contrRiemann(),                                                   // Compute contravariant Riemann curvature tensor
       compute_covRiemann(),                                                     // Compute covariant Riemann curvature tensor
       compute_covRicci(),                                                       // Compute covariant Ricci tensor
       compute_Ricciscalar();                                                    // Compute Ricci scalar curvature
  
  vector<double> compute_metricdet(),                                            // Compute the determinant of the covariant metric tensor
                 compute_invmetricdet();                                         // Compute the determinant of the contravariant metric tensor
                 
  vector<vector<double> > compute_covveccomp(vector<vector<double> >),           // Compute covariant components of given vector
                          compute_contrveccomp(vector<vector<double> >),         // Compute contravariant components of given vector
                          compute_vecfromcovcomp(vector<vector<double> >),       // Compute vector given its covariant components
                          compute_vecfromcontrcomp(vector<vector<double> >);     // Compute vector given its contravariant components
       
  void save_geometry();                                                          // Save computed geometric quantities to file
  
  vector<vector<double> > provide_compdomain(),                                  // Provide computational domain to external program
                          provide_physdomain(),                                  // Provide physical domain to external program
                          provide_covbase(),                                     // Provide covariant base vectors to external program
                          provide_contrbase(),                                   // Provide contravariant base vectors to external program
                          provide_covmetric(),                                   // Provide covariant metric tensor to external program
                          provide_contrmetric(),                                 // Provide contravariant metric tensor to external program
                          provide_firstChristoffel(),                            // Provide first Christoffel symbols to external program
                          provide_secondChristoffel(),                           // Provide second Christoffel symbols to external program
                          provide_contrRiemann(),                                // Provide contravariant Riemann curvature tensor to external program
                          provide_covRiemann(),                                  // Provide covariant Riemann curvature tensor to external program
                          provide_covRicci();                                    // Provide covariant Ricci tensor to external program
  double provide_Ricciscalar();                                                  // Provide Ricci scalar curvature to external program
};