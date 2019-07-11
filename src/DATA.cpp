#include "DATA.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * DATA.cpp
 *
 * VERSION 1.0
 * LAST UPDATED: 2016-02-19
 *
 * ~~~~~~~~
 * license:
 * ~~~~~~~~
 * Copyright (C) 2016 Prentice Lab
 *
 * This file is part of the SPLASH model.
 *
 * SPLASH is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * SPLASH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ~~~~~~~~~
 * citation:
 * ~~~~~~~~~
 * T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
 * Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
 * algorithms for simulating habitats (SPLASH): Robust indices of radiation
 * evapotranspiration and plant-available moisture, Geoscientific Model
 * Development, 2016 (in progress)
 *
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * This script reads a CSV file with input for the SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. Created read_csv and read_txt class member functions [15.02.17]
 * 02. Added year to class member variables [15.02.17]
 * 03. Added DATA header file to include [15.02.19]
 *
 * //////////////////////////////////////////////////////////////////////// */

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Constructors:
   ///////////////////////////////////////////////////////////////////// */
DATA::DATA()
    : num_lines(0), year(0)
{
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Function Definitions
   ///////////////////////////////////////////////////////////////////// */
void DATA::read_csv(string fname, int y) {
    /* ***********************************************************************
    Name:     DATA.read_csv
    Input:    - string, file name (fname)
              - int, year (y)
    Output:   None.
    Features: Reads all three daily input variables (sf, tair, and pn) for
              a single year from a CSV file that includes a headerline.
    *********************************************************************** */
    string line;
    string sf_str, tair_str, pn_str;
    double sf_val, tair_val, pn_val;
    size_t pos;
    vector<int> posit;
    int commas;

    file_name = fname;

    // Open stream to file:
    ifstream my_file( fname.c_str() );
    if (my_file.is_open()){
        int i = 0;
        while (!my_file.eof()){
            getline(my_file, line);

            // Count the number of commas in the line:
            commas = 0;
            pos = line.find_first_of(",");
            if(pos!=string::npos){
                commas++;
            }
            while (pos != string::npos){
                pos = line.find_first_of(",",pos+1);
                if(pos != string::npos){
                    commas++;
                }
            }

            if (commas == 2){
                // Save the positions for each comma:
                posit.resize(3, 0);
                pos = line.find_first_of(",");
                int j = 0;
                while (pos != string::npos)
                {
                    posit[j]=pos;
                    pos = line.find_first_of(",",pos+1);
                    j++;
                }

                // Skip headerline:
                if (i != 0) {
                    // Extract substrings based on comma locations:
                    sf_str = line.substr(0, posit[0]);
                    tair_str = line.substr((posit[0]+1),(posit[1]-posit[0]-1));
                    pn_str = line.substr((posit[1]+1),string::npos);

                    // Convert strings to doubles:
                    sf_val = atof(sf_str.c_str());
                    tair_val = atof(tair_str.c_str());
                    pn_val = atof(pn_str.c_str());

                    // Append doubles to vectors:
                    sf_vec.push_back(sf_val);
                    tair_vec.push_back(tair_val);
                    pn_vec.push_back(pn_val);
                }
            } // end if commas==2

            i++;
        } // end while file != eof
        num_lines = (i - 1);  // subtract one for last iteration
        num_lines--;          // subtract one more for the headerline
    } else {
        cout << "! Could not read " << fname << endl;
    }
    my_file.close();

    // Set class variable for year:
    if (y == -1){
        if (num_lines == 366){
            year = 2000;
        } else if (num_lines == 365) {
            year = 2001;
        }
    } else {
        year = y;
    }
}

void DATA::read_txt(string fname, string var, int y){
    /* ***********************************************************************
    Name:     DATA.read_txt
    Input:    - string, file name (fname)
              - string, variable name (var)
              - int, year (y)
    Output:   None.
    Features: Reads plain text file (no header) into one of daily input
              arrays.
    *********************************************************************** */
    string line;
    double val;

    // Add file name to list of file names:
    file_names.push_back(fname);

    // Open stream to file:
    ifstream my_file( fname.c_str() );
    if ( my_file.is_open() ){
        int i = 0;
        while ( getline(my_file, line) ){
            // Read and convert string to double:
            val = atof( line.c_str() );

            // Append to appropriate array:
            if (var == "pn"){
                pn_vec.push_back(val);
            } else if (var == "sf"){
                sf_vec.push_back(val);
            } else if (var == "tair"){
                tair_vec.push_back(val);
            }

            i++;
        } // end while file != eof
        num_lines = i;
    } else {
        cout << "! Could not read " << fname << endl;
    }
    my_file.close();

    // Set class variable for year:
    if (y == -1){
        if (num_lines == 366){
            year = 2000;
        } else if (num_lines == 365) {
            year = 2001;
        }
    } else {
        year = y;
    }
}

vector<double> DATA::get_all_sf() {
    /* ***********************************************************************
    Name:     DATA.get_all_sf
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of sun hour fractions
    *********************************************************************** */
    return sf_vec;
}

vector<double> DATA::get_all_tair() {
    /* ***********************************************************************
    Name:     DATA.get_all_tair
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of air temperature
    *********************************************************************** */
    return tair_vec;
}

vector<double> DATA::get_all_pn() {
    /* ***********************************************************************
    Name:     DATA.get_all_pn
    Input:    None
    Output:   vector<double>
    Features: This function returns the vector of precipitation
    *********************************************************************** */
    return pn_vec;
}

double DATA::get_one_sf(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_sf
    Input:    int (n)
    Output:   double
    Features: This function returns the sun hour fraction for a given index.
    *********************************************************************** */
    return sf_vec[n];
}

double DATA::get_one_tair(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_tair
    Input:    int (n)
    Output:   double
    Features: This function returns the air temperature for a given index.
    *********************************************************************** */
    return tair_vec[n];
}

double DATA::get_one_pn(int n){
    /* ***********************************************************************
    Name:     DATA.get_one_pn
    Input:    int (n)
    Output:   double
    Features: This function returns the precipitation for a given index.
    *********************************************************************** */
    return pn_vec[n];
}

int DATA::nlines(){
    /* ***********************************************************************
    Name:     DATA.nlines
    Input:    None
    Output:   int
    Features: This function returns the number of lines in the input file.
    *********************************************************************** */
    return num_lines;
}

int DATA::get_year(){
    /* ***********************************************************************
    Name:     DATA.year
    Input:    None
    Output:   int
    Features: This function returns the year of the input file.
    *********************************************************************** */
    return year;
}
