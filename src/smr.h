/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * smr.h
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
 * algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
 * transpiration and plant-available moisture, Geoscientific Model Development,
 * 2015 (in progress)
 *
 * ~~~~~~~~~~~~
 * description:
 * ~~~~~~~~~~~~
 * Header file for the smr ('s'oil 'm'oisture and 'r'unoff) struct.
 *
 * ~~~~~
 * todo:
 * ~~~~~
 * 1. Consider adding an index (for day or month of year)
 * 2. Consider adding precip and AET for annual water balance checks, i.e.:
 *    SUM_yr(Pn + Cn) = SUM_yr(AET + RO) --- or make another struct for this?
 *
 * //////////////////////////////////////////////////////////////////////// */
#ifndef SMR_H
#define SMR_H
struct smr {
    double sm;      // soil moisture, mm
    double ro;      // runoff, mm
    double swe;     // snow water equivalent, mm
    double bflow;   // baseflow, mm
    double sqout;   // drainge from the previous day
    double tdr;     // day left to drain totally water from upslope
};
#endif
