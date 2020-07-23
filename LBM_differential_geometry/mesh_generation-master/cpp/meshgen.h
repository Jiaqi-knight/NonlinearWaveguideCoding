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

#ifndef MESHGEN_H
#define MESHGEN_H

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
              A class for the generation of 2D and 3D meshes
*/
//============================================================================//
//============================================================================//

//===================================================
//==================  HEADER  =======================
//===================================================

template<class T = double>
class meshgen {

//===================================================  
//                  Variables
//===================================================
private:
 
// Input parameters
  int Dim,                                                    // Space dimensions 
      Nx,                                                     // Number of points in x-direction
      Ny,                                                     // Number of points in x-direction
      Nz;                                                     // Number of points in x-direction
	 
// Control parameters
  unsigned int float_operations;                              // Number of floating point operations performed

// Output quantities
  vector<vector<double> > comp_nodes;                         // List of nodes' coordinates in computational space, ordered by helical numbering
  vector<vector<double> > phys_nodes;                         // List of nodes' coordinates in physical space, ordered by helical numbering

  
//===================================================  
//                      Methods
//===================================================  
public:
  
  // Constructor (default)
  meshgen();
  
  // Constructor (default)
  meshgen(int D, int X, int Y, int Z = 1);
  
  //Destructor
  ~meshgen();
  
  
};

#endif