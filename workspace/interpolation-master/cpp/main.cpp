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

//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/parsers.hpp>
//#include <boost/program_options/variables_map.hpp>

#include <iostream>
#include <fstream>

#include "interp.h"
#include "interp.cpp"

using namespace std;
//using namespace boost;
//using namespace boost::program_options;

//==================  MAIN  =======================
int main(int argc, char** argv){
    
    int N, M, L, P, Q, R, K;
    
    double x1[] = {-1.5,-1.7,2.5};
    double y1[] = {pow(-1.5,2),pow(-1.7,2),pow(2.5,2)};
    
    double x2[] = {-1.5,-1.7,1.1,2.5};
    double y2[] = {pow(-1.5,2),pow(-1.7,2),pow(1.1,2),pow(2.5,2)};
    
    double xval[10];
    for(unsigned int i=0; i<10; i++){
        xval[i] = -1.4 + i*(2.4+1.4)/9;
    }
    
    return 0;
}