#include <cmath>
#include <stdio.h>

#include "global.h"
#include "EVAP.h"
#include "SOLAR.h"

using namespace std;

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * EVAP.cpp
 *
 * VERSION 1.0-r2
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
 * This script includes the definitions for the EVAP class, which calculates
 * the daily quantities of radiation, evaporation, and condensation for the
 * SPLASH model.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. finished class function definitions and debugged [15.02.07]
 * 02. removed kCw and kWm from constants (moved to SPLASH class) [15.02.17]
 * 03. changed get_cn to get_cond [15.02.17]
 * 04. added EVAP header file to include list [15.02.19]
 * 05. updated R and To [15.08.22]
 * 06. included global.h [16.01.22]
 * 07. updated documentation [16.05.27]
 * 08. addressed specific heat limitation [16.09.11]
 *
 * //////////////////////////////////////////////////////////////////////// */

/*
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Class Constructors:
/////////////////////////////////////////////////////////////////////
*/
EVAP::EVAP()
    : lat(0.0), elv(0.0), solar(SOLAR())
{
}

EVAP::EVAP(double a, double b){
    lat = a;
    elv = b;
    solar = SOLAR(lat, elv);
}

void EVAP::calculate_daily_fluxes(double sw, int n, int y, double sw_in,
                                  double tc, double slop, double asp, double snow, double nd){
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate radiation values
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    solar.calculate_daily_fluxes(n, y, sw_in, tc, slop, asp,snow, nd);
    d_sr = solar.get_vals();
    ru = d_sr.ru;
    rv = d_sr.rv;
    rw = d_sr.rw;
    rnl = d_sr.rnl;
    hn = d_sr.hn;
    rn_d = d_sr.rn_d;
    rnn_d = d_sr.rnn_d;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1.1. Assume water temperature 0.0 if air temperature < 0.0
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double tw = 0.0;
    if(tc<0.0){
        tw = 0.0; 
    }else{
        tw = tc;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate water-to-energy conversion (econ), m^3/J
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    patm = elv2pres(elv);
    s = sat_slope(tc);
    lv = enthalpy_vap(tc);
    pw = density_h2o(tw, patm);
    g = psychro(tc, patm);
    econ = s/(lv*pw*(s + g));
    //econ = s/(lv*pw*(s + 0.24*g));
    visc = calc_viscosity_h2o(tw,patm);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate daily condensation (wc), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cn = (1.0e3)*econ*abs(rnn_d);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Estimate daily equilibrium evapotranspiration (eet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eet_d = (1.0e3)*econ*rn_d;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Estimate daily potential evapotranspiration (pet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pet_d = (1.0 + Global::w)*eet_d;
    //pet_d = eet_d;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6. Calculate variable substitute (rx), (mm/hr)/(W/m^2)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rx = (3.6e6)*(1.0 + Global::w)*econ;
    // maximum instatntaneous demand (pet_max), mm/hr
    double pet_max = rx*((rw*(ru+rv))- rnl) ;
    //assume cw = pet_max
    //sw *= pet_max;
    //rx = (3.6e6)*econ;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 9. Estimate daily water supply, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //9.1. Soil matric potential to MPa
    //double psi_m_mpa = psi_m * 0.00000980665;
    // //9.2. Transpiration potential to mol/m2
    // double T_mol = pet_d / (18/pw);
    // //9.3. Leaf water potential [atm]
    // double psi_l = T_mol*(0.082*(273.15+tw)*log(0.995));
    // // potential to MPa
    // psi_l *= 0.101325;
    // //9.4. Minimum Resistance soil to leaf asuming field capacity
    // double Rp = (0.033-psi_l)/T_mol;
    //9.5. water supply mol/day/m2
    //double sw = (psi_m_mpa-psi_l)/Rp;
    //9.5. water supply mm/day
    //sw *= (18/pw);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 9. Estimate daily water supply, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /*
    //9.1. Soil matric potential to MPa
    double psi_m_mpa = psi_m * 0.00000980665;
    // //9.2. Leaf water potential at critical leaf RWC [MPa]
    // molar volume
    double v_m = 18/pw;
    double psi_l_c = ((0.082*(273.15+tw)*log(0.90))/(v_m))*0.101325;
    // //9.4. Minimum Resistance soil to leaf asuming field capacity and Leaf water potential at critical leaf RWC [MPa] 
    double Rp = (-0.033-psi_l_c)/pet_d;

    // //9.3. Leaf water potential assuming RWC at 98% [MPa]
    double psi_l = ((0.082*(273.15+tw)*log(0.98))/(v_m))*0.101325;
        
    //9.5. water supply mm/day
    double sw = (psi_m_mpa-psi_l)/Rp;
    */

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate the intersection hour angle (hi), degrees
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cos_hi = sw/(rw*rv*rx) + rnl/(rw*rv) - ru/rv;
    if (cos_hi >= 1.0) {
        // Supply exceeds demand:
        hi = 0.0;
    } else if (cos_hi <= -1.0) {
        // Supply limits demand everywhere:
        hi = 180.0;
    } else {
        hi = acos(cos_hi);
        hi /= Global::pir;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 8. Estimate daily snowmelt, mm and energy available for sublimation assume first snow melt then evaporate
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (tc >= 1.0) {
		snowmelt = min(snow,(rn_d/(pw*Global::kfus*1.1))*1000.0);
	}else{
		snowmelt=0.0;
	}
	double melt_enrg = (snowmelt/1000)*pw*Global::kfus;
	// energy available after snowmelt
    double AE = rn_d-melt_enrg;
    double sublimation = min(snowmelt,(AE*econ)*1000.0);
    melt_enrg += ((sublimation/1000.0)/econ);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 9. Estimate daily actual evapotranspiration (aet_d), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aet_d = sw*hi*Global::pir;
    aet_d += rx*rw*rv*(dsin(hn) - dsin(hi));
    aet_d += (rx*rw*ru - rx*rnl)*(hn - hi)*Global::pir;
    aet_d *= (24.0/Global::PI);
    //aet_d = sw*pet_d;
    aet_d -= (melt_enrg*econ*1000.0);
    if (aet_d < 0.0){
        aet_d = 0.0;
    }
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Function Definitions
   ///////////////////////////////////////////////////////////////////// */
double EVAP::dcos(double x){
    /* ***********************************************************************
    Name:     EVAP.dcos
    Input:    double, angle, degrees (x)
    Output:   double, cos(x*pi/180)
    Features: Calculates the cosine of an angle given in degrees
    *********************************************************************** */
    double v = cos(x*Global::pir);
    return v;
}

double EVAP::dsin(double x){
    /* ***********************************************************************
    Name:     EVAP.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = sin(x*Global::pir);
    return v;
}

double EVAP::sat_slope(double tc){
    /* ***********************************************************************
    Name:     EVAP.sat_slope
    Input:    double, air temperature (tc), degrees C
    Output:   double, slope of sat vap press temp curve (s)
    Features: Calculates the slope of the sat pressure temp curve, Pa/K
    Ref:      Eq. 13, Allen et al. (1998)
    *********************************************************************** */
    double s = exp((tc*17.269)/(tc + 237.3));
    s /= pow((tc + 237.3), 2.0);
    s *= (17.269)*(237.3)*(610.78);
    return s;
}

double EVAP::enthalpy_vap(double tc){
    /* ***********************************************************************
    Name:     EVAP.enthalpy_vap
    Input:    double, air temperature (tc), degrees C
    Output:   double, latent heat of vaporization
    Features: Calculates the enthalpy of vaporization, J/kg
    Ref:      Eq. 8, Henderson-Sellers (1984)
    *********************************************************************** */
    double lv = (tc + 273.15)/(tc + 273.15 - 33.91);
    lv = pow(lv, 2.0);
    lv *= 1.91846e6;
    return lv;
}

double EVAP::elv2pres(double z){
    /* ***********************************************************************
    Name:     EVAP.elv2pres
    Input:    double, elevation above sea level (z), m
    Output:   double, atmospheric pressure, Pa
    Features: Calculates atm. pressure for a given elevation
    Depends:  Global constants
              - Po .... base pressure     - To .... base temperature
              - L ..... temp. lapse rate  - Ma .... molecular wt. dry air
              - G ..... standard gravity  - R ..... universal gas constant
    Ref:      Allen et al. (1998)
    *********************************************************************** */
    double ep = (Global::G * Global::Ma)/(Global::R * Global::L);
    double pa = (1.0 - z*Global::L/Global::To);
    pa = pow(pa, ep);
    pa *= Global::Po;
    return pa;
}

double EVAP::density_h2o(double tc, double p){
    /* ***********************************************************************
    Name:     EVAP.density_h2o
    Input:    - double, air temperature (tc), degrees C
              - double, atmospheric pressure (p), Pa
    Output:   double, density of water, kg/m^3
    Features: Calculates density of water at a given temperature and
              pressure
    Ref:      Chen et al. (1977)
    *********************************************************************** */
    // Calculate density at 1 atm:
    double po = 0.99983952;
    po += (6.788260e-5)*tc;
    po += -(9.08659e-6)*tc*tc;
    po += (1.022130e-7)*tc*tc*tc;
    po += -(1.35439e-9)*tc*tc*tc*tc;
    po += (1.471150e-11)*tc*tc*tc*tc*tc;
    po += -(1.11663e-13)*tc*tc*tc*tc*tc*tc;
    po += (5.044070e-16)*tc*tc*tc*tc*tc*tc*tc;
    po += -(1.00659e-18)*tc*tc*tc*tc*tc*tc*tc*tc;

    // Calculate bulk modulus at 1 atm:
    double ko = 19652.17;
    ko += 148.1830*tc;
    ko += -2.29995*tc*tc;
    ko += 0.01281*tc*tc*tc;
    ko += -(4.91564e-5)*tc*tc*tc*tc;
    ko += (1.035530e-7)*tc*tc*tc*tc*tc;

    // Calculate temperature dependent coefficients:
    double ca = 3.26138;
    ca += (5.223e-4)*tc;
    ca += (1.324e-4)*tc*tc;
    ca += -(7.655e-7)*tc*tc*tc;
    ca += (8.584e-10)*tc*tc*tc*tc;

    double cb = (7.2061e-5);
    cb += -(5.8948e-6)*tc;
    cb += (8.69900e-8)*tc*tc;
    cb += -(1.0100e-9)*tc*tc*tc;
    cb += (4.3220e-12)*tc*tc*tc*tc;

    // Convert atmospheric pressure to bar (1 bar = 100000 Pa)
    double pbar = (1.0e-5)*p;

    // Calculate water density, kg/m^3
    double pw = (ko + ca*pbar + cb*pow(pbar, 2.0));
    pw /= (ko + ca*pbar + cb*pow(pbar, 2.0) - pbar);
    pw *= (1.0e3)*po;

    return pw;
}

//-Name:     -calc_viscosity_h2o-------------------------------------

// Input:    - float, ambient temperature (tc), degrees C
//           - float, ambient pressure (p), Pa
// Return:   float, viscosity of water (mu), Pa s
// Features: Calculates viscosity of water at a given temperature and 
//           pressure.
// Depends:  density_h2o
// Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
//           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
//           international formulation for the viscosity of H2O, J. Phys. 
//           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
//-----------------------------------------------------------------------

float EVAP::calc_viscosity_h2o(float tc, float p) {

  // Define reference temperature, density, and pressure values:
  float tk_ast  = 647.096;    // Kelvin
  float rho_ast = 322.0;      // kg/m^3
  float mu_ast  = 1e-6;       // Pa s

  // Get the density of water, kg/m^3
  float rho = density_h2o(tc, p);

  // Calculate dimensionless parameters:
  float tbar  = (tc + 273.15)/tk_ast;
  float tbarx = pow(tbar, 0.5);
  float tbar2 = tbar*tbar;
  float tbar3 = tbar*tbar*tbar;
  float rbar  = rho/rho_ast;

  // Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  float mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3;
  mu0 = 1e2*tbarx/mu0;

  // Create Table 3, Huber et al. (2009):
  float h_array[7][6] = {
     {0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0},  // hj0
     {0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573}, // hj1
     {-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0}, // hj2
     {0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0}, // hj3
     {-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0}, // hj4
     {0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0}, // hj5
     {0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264} // hj6
  };

  // Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  float mu1 = 0.0;
  float ctbar = (1.0/tbar) - 1.0;
  // print(paste("ctbar",ctbar))
  // for i in xrange(6):
  for (int i=0; i<6; ++i){ // in 1:6){
    float coef1 = pow(ctbar, i);
    // print(paste("i, coef1", i, coef1))
    float coef2 = 0.0;
    for (int j=0; j<7; ++j){ // in 1:7){
      coef2 = coef2 + h_array[j][i] * pow(rbar - 1.0, j);
    }
    mu1 = mu1 + coef1 * coef2;    
  }
  mu1 = exp( rbar * mu1 );
  // print(paste("mu1",mu1))

  // Calculate mu_bar (Eq. 2, Huber et al., 2009)
  //   assumes mu2 = 1
  float mu_bar = mu0 * mu1;

  // Calculate mu (Eq. 1, Huber et al., 2009)
  float mu = mu_bar * mu_ast;    // Pa s

  return mu; 
}


double EVAP::psychro(double tc, double p){
    /* ***********************************************************************
    Name:     EVAP.psychro
    Input:    - double, air temperature (tc), degrees C
              - double, atm. pressure (p), Pa
    Output:   double, psychrometric constant, Pa/K
    Features: Calculates the psychrometric constant for a given temperature
              and pressure
    Depends:  Global constants:
              - kMa
              - kMv
    Refs:     Allen et al. (1998); Tsilingiris (2008)
    *********************************************************************** */
    // Calculate the specific heat capacity of water, J/kg/K
    double cp = specific_heat(tc);

    // Calculate latent heat of vaporization, J/kg
    double lv = enthalpy_vap(tc);

    // Calculate psychrometric constant, Pa/K
    //   Eq. 8, Allen et al. (1998)
    double ps = (Global::Ma*cp*p)/(Global::Mv*lv);

    return ps;
}

double EVAP::specific_heat(double tc){
    /* ***********************************************************************
    Name:     EVAP.specific_heat
    Input:    double, air temperature (tc), degrees C
    Output:   double, specific heat of moist air, J/kg/K
    Features: Calculates the specific heat of moist air for a given temperature
    Refs:     Eq. 47, Tsilingiris (2008); valid only for air temp 0--100 deg C
    *********************************************************************** */
    double cp;
    if (tc < 0) {
        cp = 1004.5714270;
    } else if (tc > 100) {
        cp = 2031.2260590;
    } else {
        cp = 1.0045714270;
        cp += (2.050632750e-3)*tc;
        cp += -(1.631537093e-4)*tc*tc;
        cp += (6.212300300e-6)*tc*tc*tc;
        cp += -(8.830478888e-8)*tc*tc*tc*tc;
        cp += (5.071307038e-10)*tc*tc*tc*tc*tc;
        cp *= (1.0e3);
    }

    return cp;
}

etr EVAP::get_vals(){
    /* ***********************************************************************
    Name:     EVAP.get_vals
    Input:    None
    Output:   etr struct, daily ET and radiation values
    Features: Returns the current set of daily ET and radiation values
    *********************************************************************** */
    d_etr.cond = cn;
    d_etr.eet = eet_d;
    d_etr.pet = pet_d;
    d_etr.aet = aet_d;
    d_etr.snowmelt = snowmelt;
    d_etr.econ = econ;
    d_etr.pw = pw;
    d_etr.rn_d = rn_d;
    d_etr.visc = visc;
    
    return d_etr;
}

void EVAP::display(){
    /* ***********************************************************************
    Name:     EVAP.display
    Input:    None
    Output:   None
    Features: Prints the current set of daily STASH variables
    *********************************************************************** */
    solar.display();

    printf("EVAP variable list:\n");
    printf("  s: %0.6f Pa/K\n", s);
    printf("  lv: %0.6f MJ/kg\n", (1.0e-6)*lv);
    printf("  Patm: %0.6f bar\n", (1.0e-5)*patm);
    printf("  pw: %0.6f kg/m^3\n", pw);
    printf("  gamma: %0.6f Pa/K\n", g);
    printf("  Econ: %0.6f mm^3/J\n", (1.0e9)*econ);
    printf("  Cn: %0.6f mm\n", cn);
    printf("  rx: %0.6f\n", rx);
    printf("  hi: %0.6f degrees\n", hi);
    printf("  EET: %0.6f mm\n", eet_d);
    printf("  PET: %0.6f mm\n", pet_d);
    printf("  AET: %0.6f mm\n", aet_d);
}
