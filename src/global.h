/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * global.h
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
 * Header file for the global constants.
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef GLOBAL_H
#define GLOBAL_H

namespace Global {
    // forward declarations of global constants
    extern const double A;        // constant for Rnl
    extern const double alb_sw;   // shortwave albedo
    extern const double alb_vis;  // visible light albedo
    extern const double b;        // constant for Rnl
    extern const double c;        // cloudy transmittivity
    extern const double d;        // angular coefficient of transmittivity
    extern const double fFEC;     // from flux to energy conversion, umol/J
    extern const double kfus;     // latent heat of fusion J/Kg (Monteith & Unsworth, 1990)
    extern const double G;        // gravitational acceleration, m/s^2
    extern const double Gsc;      // solar constant, W/m^2
    extern const double L;        // temperature lapse rate, K/m
    extern const double Ma;       // molecular weight of dry air, kg/mol
    extern const double Mv;       // molecular weight of water vapor, kg/mol
    extern const double Po;       // standard atmosphere, Pa
    extern const double R;        // universal gas constant, J/mol/K
    extern const double To;       // base temperature, K
    extern const double w;        // entrainment factor
    extern const double Cw;       // supply constant, mm/hr
    extern const double Wm;       // soil moisture capacity, mm
    extern const double PI;       // pi
    extern const double pir;      // degrees to radians
    extern const double fluidity; // fluidity, mm^2 (density*gravity)/viscosity
    // Paleoclimate variables:
    extern const double e;        // eccentricity
    extern const double eps;      // obliquity, degrees
    extern const double omega;    // longitude of perihelion, degrees
}
#endif
