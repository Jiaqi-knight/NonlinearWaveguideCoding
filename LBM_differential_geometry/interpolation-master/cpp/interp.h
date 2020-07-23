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

#ifndef INTERP_H
#define INTERP_H

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
    A class providing tools for 1D, 2D, 3D Lagrange and Hermite interpolation
*/
//============================================================================//
//============================================================================//


//===================================================
//==================  HEADER  =======================
//===================================================

template<class T = double>
class interp {

    //===================================================  
    //                  Variables
    //===================================================
    private:

    // Control parameters
    bool measure;                                                               // True --> count floating point and integer operations, false otherwise
    unsigned int floats;                                                        // Number of floating point operations performed
    unsigned int ips;                                                           // Number of instructions (integer operations) performed

    // Output quantities
    int D,                                                                      // Space dimension
        N,                                                                      // Number of elements of x vector (interpolation domain values)
        M,                                                                      // Number of elements of y vector (interpolation domain values)
        L,                                                                      // Number of elements of z vector (interpolation domain values)
        P,                                                                      // Number of elements of xval vector  (values for functional evaluation of interpolated function)
        Q,                                                                      // Number of elements of yval vector  (values for functional evaluation of interpolated function)
        R,                                                                      // Number of elements of zval vector  (values for functional evaluation of interpolated function)
        K;                                                                      // Maximum order of derivatives
    
    vector<T> x,                                                                // x vector (interpolation domain values)
              y,                                                                // y vector (interpolation domain values)
              z,                                                                // z vector (interpolation domain values)
              xval,                                                             // xval vector  (values for functional evaluation of interpolated function)
              yval,                                                             // yval vector  (values for functional evaluation of interpolated function)
              zval;                                                             // zval vector  (values for functional evaluation of interpolated function)
  
    vector<vector<vector<vector<T> > > > w;                                     // {w_i = f(x_i,y_i,z_i)
                                                                                // dw/dx|_i = df/dx(x_i,y_i,z_i)
                                                                                // d^2w/dx^2|_i = d^2f/dx^2(x_i,y_i,z_i)
                                                                                // ...
                                                                                // d^Kw/dx^K|_i = d^Kf/dx^K(x_i,y_i,z_i)}
    
    vector<vector<vector<T> > > f;                                              // f_k = f~(xval_k,yval_k,zval_k)
  
    //===================================================  
    //                      Methods
    //===================================================  
    public:
  
    // Constructor (default)
    interp();
    
    // Constructor (1D)
    interp(T xinp[], T winp[], T xvalinp[], int Ninp, int Pinp, int Kinp);
    
    // Constructor (1D measure performances)
    interp(T xinp[], T winp[], T xvalinp[], int Ninp, int Pinp, int Kinp, bool meas);
    
    // Constructor (1D)
    interp(vector<T> xinp, vector<vector<T> > winp, vector<T> xvalinp);
    
    // Constructor (1D measure performances)
    interp(vector<T> xinp, vector<vector<T> > winp, vector<T> xvalinp, bool meas);
    
    // Constructor (2D)
    interp(T xinp[], T yinp[], T winp[], T xvalinp[], T yvalinp[], int Ninp, int Minp, int Pinp, int Qinp, int Kinp);
    
    // Constructor (2D measure performances)
    interp(T xinp[], T yinp[], T winp[], T xvalinp[], T yvalinp[], int Ninp, int Minp, int Pinp, int Qinp, int Kinp, bool meas);
    
    // Constructor (2D)
    interp(vector<T> xinp, vector<T> yinp, vector<vector<vector<T> > > winp, vector<T> xvalinp, vector<T> yvalinp);
    
    // Constructor (2D measure performances)
    interp(vector<T> xinp, vector<T> yinp, vector<vector<vector<T> > > winp, vector<T> xvalinp, vector<T> yvalinp, bool meas);
    
    // Constructor (3D)
    interp(T xinp[], T yinp[], T zinp[], T winp[], T xvalinp[], T yvalinp[], T zvalinp[], int Ninp, int Minp, int Linp, int Pinp, int Qinp, int Rinp, int Kinp);
    
    // Constructor (3D measure performances)
    interp(T xinp[], T yinp[], T zinp[], T winp[], T xvalinp[], T yvalinp[], T zvalinp[], int Ninp, int Minp, int Linp, int Pinp, int Qinp, int Rinp, int Kinp, bool meas);
  
    // Constructor (3D)
    interp(vector<T> xinp, vector<T> yinp, vector<T> zinp, vector<vector<vector<vector<T> > > > winp, vector<T> xvalinp, vector<T> yvalinp, vector<T> zvalinp);
    
    // Constructor (3D measure performances)
    interp(vector<T> xinp, vector<T> yinp, vector<T> zinp, vector<vector<vector<vector<T> > > > winp, vector<T> xvalinp, vector<T> yvalinp, vector<T> zvalinp, bool meas);
  
  
    // Destructor
    ~interp();
  
    // General tools
    
    vector<vector<T> > NewtonDividedDifferences(vector<T> xinp, vector<T> yinp);            // Build Lagrange interpolator using Newton's divided differences
    
    T LagrangeInterpolant(vector<T> xinp, T xvalinp, unsigned int index);                   // Build Lagrange interpolator using Lagrange's form
    vector<T> LagrangeInterpolant(vector<T> xinp, vector<T> xvalinp, unsigned int index);   // Build Lagrange interpolator using Lagrange's form
    
    T HermiteInterpolant(vector<T> xinp, T xvalinp, unsigned int index);                   // Build Hermite interpolator
    
    void Lagrange1D(),                                                                      // 1D Lagrange interpolation
         Lagrange2D(),                                                                      // 2D Lagrange interpolation
         Lagrange3D(),                                                                      // 3D Lagrange interpolation
         Hermite1D(),                                                                       // 1D Hermite interpolation
         Hermite2D(),                                                                       // 2D Hermite interpolation
         Hermite3D(),                                                                       // 3D Hermite interpolation
         interpolation(int type = 0);                                                       // Perform interpolation of data (type = 0 --> Lagrange interpolation, default; type = 1 --> Hermite interpolation)
    
    //API
    
    vector<T> f1D();                                                                        // Provide 1D interpolated function to external program                               
    vector<vector<T> > f2D();                                                               // Provide 2D interpolated function to external program
    vector<vector<vector<T> > > f3D();                                                      // Provide 3D interpolated function to external program
  
    // Overloading of operators
    
    template<class U>
    friend ostream& operator<<(ostream& output, const interp<U>& rhs);
  
};

#endif