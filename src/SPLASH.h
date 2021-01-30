#include <vector>

#include "etr.h"
#include "smr.h"
#include "DATA.h"
#include "EVAP.h"
using namespace std;
/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SPLASH.h
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
 * This is the header file for the C++ SPLASH class.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. moved header includes from cpp file here [15.02.19]
 * 02. created header guard [15.02.19]
 * 03. added etr structure [15.02.19]
 * 04. added smr structure [15.02.19]
 * 05. added DATA header to include list [15.02.19]
 * 06. added quick_run & spin_up functions [15.02.19]
 * 07. removed constants; now in global.h [16.01.22]
 * 08. made wn_vec a private variable [16.02.06]
 * 09. added print daily wn function [16.02.06]
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SPLASH_H
#define SPLASH_H

#include <Rcpp.h>

class SPLASH {
    private:
        // Static variables:
        double lat;                       // latitude, degrees
        double elv;                       // elevation, meters
        double econ, pw, rn_d;           // variable substitutions
        double visc;                      // viscosity Pa s
        // Daily status variables:
        EVAP evap;                        // daily EVAP class
        etr dvap;                         // daily etr struct
        smr dsoil;                        // daily smr struct
        double precip;                    // daily precipitation, mm

        // Daily soil moisture
        std::vector<double> wn_vec;

    public:
        // Constructors:
        SPLASH(double a, double b);

        // Functions:
        void quick_run(int n, int y, double wn, double sw_in, double tc,
                       double pn, smr &dsm, double slop, double asp, double snow, double snowfall, vector <double> &soil_info,
                       double qin, double td, double nd);
        void run_one_day(int n, int y, double wn, double sw_in, double tc,
                         double pn, smr &dsoil, double slop, double asp,double snow, double snowfall, vector <double> &soil_info,double qin, double td, double nd);
        Rcpp::List spin_up(int n, int y, vector<double> &sw_in, vector <double> &tair, vector <double> &pn, double slop, double asp, vector <double> &snowfall, vector <double> &soil_info);        
        Rcpp::List run_one_year(int n, int y, vector<double> &sw_in, vector <double> &tair, vector <double> &pn, vector <double> &wn_vec, double slop, double asp, 
                                vector <double> &snow, vector <double> &snowfall, vector <double> &soil_info, vector <double> &qin_vec,vector <double> &td_vec,vector <double> &nds);
        double dcos(double x);
        double dsin(double x);
        double dtan(double x);
        double moist_surf(double depth,double z,double bub_p,double wn,double SAT,double RES,double lambda);
        double inf_GA(double bub_press,double theta_i,double Ksat,double SAT,double lambda,double P,double tdur,double slop);
        double get_wtd(double wn, double depth, double bub_press,double theta_s,double theta_r,double lambda);
        double get_lat();
        double get_elv();
        void print_daily_wn();
        void print_vals();
};
#endif
