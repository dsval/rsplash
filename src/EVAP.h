#include "etr.h"
#include "srad.h"
#include "SOLAR.h"

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * EVAP.h
 *
 * VERSION 1.0
 * LAST UPDATED: 2016-09-11
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
 * This is the C++ header file for the EVAP class.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. removed kCw and kWm from constants (moved to SPLASH class) [15.02.17]
 * 02. changed get_cn to get_cond [15.02.17]
 * 03. added all necessary inputs from cpp file [15.02.19]
 * 04. created header guard [15.02.19]
 * 05. added etr struct to include list [15.02.19]
 * 06. removed the unnecessary get() functions [15.02.19]
 * 07. removed constants; now in global.h [16.01.22]
 * 08. addressed specific heat limitation [16.09.11]
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef EVAP_H
#define EVAP_H
class EVAP {
    private:
        double lat;                // latitude, degrees
        double elv;                // elevation, m
        SOLAR solar;               // SOLAR class
        double ru, rv, rw, rx;     // variable substitutions
        double hn, hi;             // hour angles, degrees
        double rnl;                // longwave radiation flux, W/m^2
        double rn_d;               // daily net radiation, J/m^2
        double rnn_d;              // daily nighttime net radiation, J/m^2
        double s;                  // slope sat. vap. press. temp. curve, Pa/K
        double lv;                 // enthalpy of vaporization, J/kg
        double pw;                 // density of water, kg/m^3
        double g;                  // psychrometric constant, Pa/K
        double patm;               // atmospheric pressure, Pa
        double econ;               // water-to-energy conversion, m^3/J
        double cn;                 // daily condensation, mm
        double eet_d;              // daily equilibrium ET, mm/d
        double pet_d;              // daily potential ET, mm/d
        double aet_d;              // daily actual ET, mm/d
        double cos_hi;             // cosine of hour angle, hi
        double snowmelt;           // daily snowmelt mm/d
        double visc;               // viscosity Pa s
        srad d_sr;                 // daily srad struct
        etr d_etr;                 // daily etr struct

        // Functions:
        double dcos(double x);
        double dsin(double x);
        double sat_slope(double tc);
        double enthalpy_vap(double tc);
        double elv2pres(double z);
        double density_h2o(double tc, double p);
        float calc_viscosity_h2o(float tc, float p);
        double psychro(double tc, double p);
        double specific_heat(double tc);
        
    public:
        // Constructors
        EVAP();
        EVAP(double a, double b);

        // General purpose functions:
        void calculate_daily_fluxes(double sw, int n, int y, double sw_in,
                                    double tc, double slop, double asp,double snow);
        etr get_vals();
        void display();
};
#endif
