#include <cstdlib>
#include <stdio.h>
#include <vector>

#include "global.h"
#include "SPLASH.h"

using namespace std;

#include <Rcpp.h>
using namespace Rcpp;
RCPP_MODULE(splash_module) {
	class_<SPLASH>( "SPLASH" )
	.constructor<double,double>()
	.method( "run_one_year", &SPLASH::run_one_year )
    .method( "moist_surf", &SPLASH::moist_surf )
    .method( "inf_GA", &SPLASH::inf_GA )
	.method( "spin_up", &SPLASH::spin_up )
	;
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 * SPLASH.cpp
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
 * This class updates daily quantities of radiation, evapotranspiration, soil
 * moisture and runoff based on the SPLASH methodology.
 *
 * ~~~~~~~~~~
 * changelog:
 * ~~~~~~~~~~
 * 01. added header to include list [15.02.19]
 * 02. added iostream to include list [15.02.19]
 * 03. created print_vals function [15.02.19]
 * 04. added smr struct (dsoil) [15.02.19]
 * 05. added vector to include list [15.02.19]
 * 06. created quick_run & spin_up functions [15.02.19]
 * 07. included global.h [16.01.22]
 * 08. wn_vec is now a class variable [16.02.06]
 *
 * //////////////////////////////////////////////////////////////////////// */

//RCPP_MODULE(unif_module) {
//class_<SPLASH>("SPLASH")
//.constructor<double,double>()
//.method("quick_run", &SPLASH::quick_run)
//.method("run_one_day", &SPLASH::run_one_day)
//;
//}


/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Constructors:
   ///////////////////////////////////////////////////////////////////// */
SPLASH::SPLASH(double a, double b)
    : lat(a), elv(b), precip(0.0)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialize class & structure values:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap = EVAP(a, b);

    dvap.cond = 0.0;
    dvap.eet = 0.0;
    dvap.pet = 0.0;
    dvap.aet = 0.0;
    dvap.snowmelt = 0.0;
    dvap.econ = 0.0;
    dvap.pw = 0.0;
    dvap.rn_d = 0.0;

    dsoil.sm = 0.0;
    dsoil.ro = 0.0;
    dsoil.swe = 0.0;
    dsoil.sqout = 0.0;
    dsoil.tdr = 0.0;
}

/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   Class Function Definitions
   ///////////////////////////////////////////////////////////////////// */
double SPLASH::dcos(double x){
    /* ***********************************************************************
    Name:     EVAP.dcos
    Input:    double, angle, degrees (x)
    Output:   double, cos(x*pi/180)
    Features: Calculates the cosine of an angle given in degrees
    *********************************************************************** */
    double v = cos(x*Global::pir);
    return v;
}

double SPLASH::dsin(double x){
    /* ***********************************************************************
    Name:     EVAP.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = sin(x*Global::pir);
    return v;
}

double SPLASH::dtan(double x){
    /* ***********************************************************************
    Name:     EVAP.dsin
    Input:    double, angle, degrees (x)
    Output:   double, sin(x*pi/180)
    Features: Calculates the sine of an angle given in degrees
    *********************************************************************** */
    double v = tan(x*Global::pir);
    return v;
}

void SPLASH::quick_run(int n, int y, double wn, double sw_in, double tc,
                       double pn, smr &dsm, double slop, double asp, double snow, double snowfall, vector <double> &soil_info,
                       double qin, double td) {
    /* ***********************************************************************
    Name:     SPLASH.quick_run
    Input:    - int, day of year (n)
              - int, year (y)
              - double, previous day's soil moisture, mm (wn)
              - double, daily sunshine fraction (sf)
              - double, daily air temperature, deg C (tc)
              - double, daily precipitation, mm (pn)
              - smr, daily soil moisture & runoff
    Output:   None.
    Features: Calculates daily soil moisture and runoff based on STASH
              methods.
    *********************************************************************** */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // . Preprocess
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //read soil hydro points in mm
    //saturation (water  content mm at 0 KPa)
    double SAT = soil_info[0];
    //wilting point (water content mm at 1500 KPa)
    double WP = soil_info[1];
    //field capacity (water content mm at 33KPa)
    double FC = soil_info[2];
    //saturated hydraulic conductivity(mm/hr)
    double Ksat = soil_info[3];
    //slope of the log-log curve soil moisture-matric potetial from brooks and Corey (1964)
    double lambda = soil_info[4];
    double depth = soil_info[5];
    //bubbling pressure/capillarity fringe (mm)
    double bub_press = soil_info[6];
    //residual water content, test as WP?
    double RES = soil_info[1];
    //double RES = soil_info[7];
    //upslope area
    double Au = soil_info[8];
    //Pixel Area
    double Ai = soil_info[9];
    //neighbourhood  cells
    double cellin = soil_info[10];
    double cellout = soil_info[11];
    double theta_s = SAT/(depth *1000.0);
    double theta_r = RES/(depth *1000.0);
    double theta_fc = FC/(depth *1000.0);
    double theta_wp = WP/(depth *1000.0);
    double theta_q0 = theta_wp + 0.001;
    double theta_qs = theta_s;
    double theta_i = (wn)/(depth*1000.0);
    // correct theta_i for NA error reaching boundary conditions
    if (theta_i>=theta_s){
        theta_i = theta_s - 0.001;
    } else if (theta_i<=theta_r){
        theta_i = theta_wp + 0.001;
    }
    double alph = 4.0 + 2.0*lambda;
    // 7.1. assume hydraulic gradient as the slope
    double hyd_grad = (dtan(slop));
    // assuming octogonal cells, ncells draining will have n octagon sides of length[m]:
    double sid_oct = sqrt(Ai/(2.0*(1+sqrt(2.0))));
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 00. Compute real FC
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //double coeff_A = exp(log(33.0) + (1.0/lambda)*log(FC/(depth *1000.0)));
    //double Wmax = pow((coeff_A/(depth*1000)), (1.0/((1/lambda)+1.0))) * (depth *1000.0) ;
    
    double psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // 1. Calculate evaporative supply rate (sw), mm/h using calculate water in the first 50cm of depth
    // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1.0 get wtd
    double wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // 2.0 get unsaturated/saturated depth within root zone if wtd is shallow, MAX ROOTING DEPTH = 2.0 m
    /////////////////////////////////////////////////////////////////////////////////////
    double z_uns =0.0;
    double sat_swc = 0.0;
    if (wtd<=2.0){
        z_uns = wtd*1000.0;
        sat_swc = theta_s*(2.0-wtd)*1000.0;
    }else{
        z_uns = 2.0*1000.0;
        sat_swc = 0.0;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // 4.0 calculate water in the within root zone
    /////////////////////////////////////////////////////////////////////////////////////
    double w_uns_z = (theta_r*z_uns) + (((psi_m+z_uns)*(theta_r-theta_s)*pow((bub_press/(psi_m+z_uns)),lambda))/(lambda-1));
    double w_uns_0 = (theta_r*0.0) + (((psi_m+0.0)*(theta_r-theta_s)*pow((bub_press/(psi_m+0.0)),lambda))/(lambda-1));
	double w_uns = w_uns_z-w_uns_0;
    double w_z = w_uns + sat_swc ;
    /////////////////////////////////////////////////////////////////////////////////////
    // 3.0 calculate real FC
    /////////////////////////////////////////////////////////////////////////////////////
    double coeff_A = exp(log(33.0) + (1.0/lambda)*log(theta_fc));
    double Wmax = pow((coeff_A/(z_uns)), (1.0/((1/lambda)+1.0))) * (z_uns) ;
    /////////////////////////////////////////////////////////////////////////////////////
    // 5.0 calculate supply rate (sw)
    /////////////////////////////////////////////////////////////////////////////////////
    double sw = 0.0;
    if (depth>2.0){
        double RES_z = theta_r * z_uns;
        sw = Global::Cw*((w_z-RES_z)/(Wmax-RES_z));
    } else{
        Wmax = pow((coeff_A/(depth*1000)), (1.0/((1/lambda)+1.0))) * (depth *1000.0) ;
        sw = Global::Cw*((wn-RES)/(Wmax-RES));
    }
   
    //double sw = Global::Cw*((wn)/(SAT));
    if (sw < 0.0 || isnan(sw)==1) {
         sw = 0.0;
     } else if (sw > Global::Cw){
         sw = Global::Cw;
     }


    //double theta_uns = w_uns/( dept_rz*1000.0);
    //psi_m = bub_press/pow((((theta_uns-theta_r)/(theta_s-theta_r))),(1/lambda));

    // // 7.2.3 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    
    // double w_sat = theta_s * (depth-wtd) * 1000.0;
    // //double theta_u = (w_uns)/((wtd)*1000.0);
    // double theta_u = moist_surf(depth,10.0,bub_press,wn,SAT,RES,lambda);
    
    // //
    // double dept_rz = 0.0;

    // if (depth > 0.5){
    //     dept_rz = 0.5;
    // }else{
    //     dept_rz = depth;
    // }
    // double z_uns = dept_rz*1000.0;
    
    // double w_uns_z = (theta_r*z_uns) + (((psi_m+z_uns)*(theta_r-theta_s)*pow((bub_press/(psi_m+z_uns)),lambda))/(lambda-1));
    // double w_uns_0 = (theta_r*0.0) + (((psi_m+0.0)*(theta_r-theta_s)*pow((bub_press/(psi_m+0.0)),lambda))/(lambda-1));
	// double w_uns = w_uns_z-w_uns_0;
    // double theta_uns = w_uns/( dept_rz*1000.0);
    // psi_m = bub_press/pow((((theta_uns-theta_r)/(theta_s-theta_r))),(1/lambda));
    // get wp and fc
    // double wp_uns = theta_wp * z_uns;
    // double fc_uns = theta_fc * z_uns;
    // double sw = Global::Cw*((w_uns-wp_uns)/(fc_uns-wp_uns));

    // if (sw < 0.0 || isnan(sw)==1) {
    //      sw = 0.0;
    //  } else if (sw > Global::Cw){
    //      sw = Global::Cw;
    //  }


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1. Calculate evaporative supply rate (sw), mm/h Original
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // org splash
    // double sw = Global::Cw*((wn-RES)/(Wmax-RES));
    // //double sw = Global::Cw*((wn)/(SAT));
    // if (sw < 0.0 || isnan(sw)==1) {
    //      sw = 0.0;
    //  } else if (sw > Global::Cw){
    //      sw = Global::Cw;
    //  }
   

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //First 10 cm
    //double theta_aet = moist_surf(depth,10,bub_press,wn,SAT,RES,lambda);
    //double sw = Global::Cw*((theta_aet-theta_r)/theta_fc);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Maatselar approach
    // double theta_l = (theta_wp) + 0.9*(theta_s-theta_wp);
    // double TH_l = (theta_l-theta_r)/(theta_s-theta_wp);
    // double TH_w = (theta_wp-theta_r)/(theta_s-theta_wp);
    // double TH = (theta_i-theta_r)/(theta_s-theta_wp);
    // double f_aet = (pow(TH,alph)-pow(TH_w,alph)) /(pow(TH_l,alph)-pow(TH_w,alph));
    //double sw = Global::Cw*(f_aet);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Feddes approach
    // double h1 = 100.0;
    // double h2 = -250.0;
    // double h3 = -2000.0;
    // double h4 = -150000.0;
    // double sw = 0.0;
    // if((psi_m >= h1)||(psi_m <= h4)){
    //     sw = 0.0;
    // }else if ((psi_m >= h2) && (psi_m < h1)) {
    //     sw = (psi_m-h1)/(h2-h1);
    // }else if ((psi_m > h3) && (psi_m <= h2 )) {
    //     sw = 1.0;
    // }else if ((psi_m > h4) && (psi_m <= h3 )) {
    //     sw = (psi_m-h4)/(h3 - h4);
    // }
    // sw *= Global::Cw;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //Campbel approach
    // double awc = (wn-RES)/(Wmax-RES);
    // double sw = 1.0 - pow((1.0 + 1.3 * awc ),(-1.0/lambda));
    
    // sw *= Global::Cw;
    // if (sw < 0.0 || isnan(sw)==1) {
    //      sw = 0.0;
    //  } else if (sw > Global::Cw){
    //      sw = Global::Cw;
    //  }

    //Bonan, 2014 gradient assuming 2mmol root conductance
    // double psi_m_mpa = psi_m * 0.00000980665;
    // double sw = ((psi_m_mpa+1.5)*2.0*(18.0*86400.0)/(1e6)); //* pow((theta_i/theta_s),(3+(2/lambda))) ;
    //if (sw < 0.0 || isnan(sw)==1) {
    //    sw = 0.0;
    //}

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Accumulate snow (mm) first so feedback on albedo happens
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    snow += snowfall;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 3. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap.calculate_daily_fluxes(sw, n, y, sw_in, tc, slop, asp,snow);
    etr dn = evap.get_vals();
    econ = dn.econ;
    pw = dn.pw;
    rn_d = dn.rn_d;
    snow -= dn.snowmelt;
    visc = dn.visc;
    
    double melt_enrg = (dn.snowmelt/1000.0)*pw*Global::kfus;
	// energy available after snowmelt
    double AE = rn_d-melt_enrg;
    double sublimation = min(dn.snowmelt,(AE*econ)*1000.0);
    double snowmelt = dn.snowmelt - sublimation;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate daily drainage, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     // correct Ksat for fuid viscosity
	double int_perm = Ksat/Global::fluidity;

	double Ksat_visc = int_perm*((pw*Global::G)/visc)*3.6;
    //double Ksat_visc = Ksat;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4.1 calculate inflow
    double inflow = pn + 0.1*dn.cond + snowmelt;
    // 4.2 calculate skin moisture at 5 cm
    
    double surf_moist = moist_surf(depth,5.0,bub_press,wn,SAT,RES,lambda);
    // 4.2 calculate infiltration assuming storm duration 6hrs
    double infi = inf_GA(bub_press,surf_moist,Ksat,theta_s,lambda,inflow,6.0,slop);
    // 4.2 calculate Dunne runoff
    double ro_d = max(inflow-infi,0.0);
    // 4.2 calculate recharge
    double R = infi - dn.aet;
    // 4.2 update soil moisture
    //double sm = wn + pn + dn.cond - dn.aet + snowmelt;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.4. Estimate recession constant Kb
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.4.1 estimate baseflow at ~wilting point
    // 7.4.1.1 calculate matric potential (mmH2O)
    double psi_q0 = bub_press/pow((((theta_q0-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.4.1.2 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    double wtd_q0 = ((bub_press-psi_q0)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd_q0 < 0.0 || isnan(wtd_q0)==1){
        wtd_q0 = 0.0;
    }else if (wtd_q0 > depth){
        wtd_q0 = depth;
    }
    //7.4.1.3 calculate pixel's cross sectional area for baseflow output
    double Acs_out_q0 = (depth - wtd_q0) * sid_oct * cellout;
    //7.4.1.4 calculate transmitance over the unsaturated part of th profile [mm^2/h]
    double T_q0 =  (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_q0), (3.0*lambda+1.0))-pow((bub_press/(psi_q0 +(wtd_q0*1000.0))), (3.0*lambda+1.0)) );
    //7.4.1.5 adjust transmitance to m3/day
    double Q_q0 = T_q0 * hyd_grad* ((24.0 * sid_oct * cellout)/(1.0e6));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 7.4.2 estimate baseflow at saturation Q_qs
    // 7.4.2.1 calculate matric potential (mmH2O)
    double psi_qs = bub_press/pow((((theta_s-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.4.2.2 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    double wtd_qs = ((bub_press-psi_qs)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd_qs < 0.0 || isnan(wtd_qs)==1){
        wtd_qs = 0.0;
    }else if (wtd_qs > depth){
        wtd_qs = depth;
    }
    //7.4.2.3 calculate pixel's cross sectional area for baseflow output [m^2]
    double Acs_out_qs = (depth) * sid_oct * cellout;
    //7.4.2.4 calculate transmitance over the unsaturated part of th profile [mm^2/h]
    double T_qs = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_qs), (3.0*lambda+1.0))-pow((bub_press/(psi_qs +(wtd_qs*1000.0))), (3.0*lambda+1.0)) );
    //7.4.2.5 adjust transmitance to m3/day
    //double Q_qs = (T_qs *hyd_grad *((24.0 * sqrt(Ai))/(1.0e6)))+(hyd_grad * Ksat_visc*24.0* (Acs_out_qs)/1000);
    double Q_qs = (hyd_grad * Ksat_visc*24.0* (Acs_out_qs)/1000.0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // 7.4.2 compute Kb
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    double Kb = exp((Q_q0-Q_qs)/((SAT-WP)*(Ai/1000.0)));
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.2. Estimate Transmitance at the beggining of the day
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.2.1 calculate volumetric water content for the whole soil column
    //theta_i = (wn)/(depth*1000.0);
    // 7.2.2 calculate matric potential (mmH2O)
    psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.2.3 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
    //7.2.4 calculate pixel's cross sectional area saturated part
    double Acs_out = (depth - wtd) * cellin * sid_oct;
    // transmitance over the unsaturated part of th profile [mm^2/h]
    double To_uns = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_m), (3.0*lambda+1.0))-pow((bub_press/(psi_m +(wtd*1000.0))), (3.0*lambda+1.0)) );
    // adjust transmitance to m3/day
    double Qo_uns = To_uns * ((24.0 * cellin * sid_oct)/(1.0e6));

    // adjust to flux density  and daily timestep [mm/day]
    To_uns *= ((24.0 * cellin * sid_oct)/(1000.0 * Ai));
    
    // transmitance over the saturated part of the profile adjusted to flux density [mm/day]
    double To_sat = Ksat_visc*24.0* (Acs_out/Ai);
     // adjust transmitance to m3/day saturated part
    double Qo_sat = Ksat_visc*24.0* (Acs_out)/1000.0;
    // total transmitance at initial soil moisture wn
    double To = (To_sat+To_uns)*hyd_grad;
    // total Qo
     double Qo = (Qo_sat+Qo_uns)*hyd_grad;  
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.6. Estimate q_in from upslope (flux density) [mm/day]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double q_in = 0.0; 
      

    if( (td == 0.0) ){ 
      q_in = 0.0;
    } else {
        q_in = qin*Kb;
    }

    double sm = wn + R + q_in;
        

   // ~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Calculate hortonian runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro_h = 0.0;
    //declare some variables
   double dummy = wn + R;
     if (dummy > SAT){
        // Bucket is too full:
        //   allocate excess water to runoff
        
        ro_h = (dummy - SAT);
        dummy = SAT;
        R -= ro_h;
   
    } else if (dummy < RES){
        // Bucket is too empty:
        //   set soil moisture & runoff to zero
        ro_h = 0.0;
        dummy = RES;
    } 
    if (sm > SAT){
        // Bucket is too full:
        //   allocate excess water to runoff
        
        //ro_h = (sm - SAT);
        sm = SAT;
        //R -= ro_h;
   
    } else if (sm < RES){
        // Bucket is too empty:
        //   set soil moisture & runoff to zero
        //ro_h = 0.0;
        sm = RES;
    } 
    // 5.1 update runoff
    double ro  = ro_d + ro_h;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6. Calculate drainage vadose sone (unsaturated) (used only to compare against field measurements at low depths)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // double subro = Ksat_visc*pow((sm/SAT), (3+(2/lambda)))*24*(dsin(slop));
    // sm -= subro;
    // double q_drain = Ksat_visc*pow((sm/SAT), (3+(2/lambda)))*24*(dcos(slop));
    // sm -=q_drain;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate transmitance previous day To
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.3. Estimate Transmitance after recharge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.3.1 calculate volumetric water content for the whole soil column
    theta_i = (sm)/(depth*1000.0);
    // 7.3.2 calculate matric potential (mmH2O)
    psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.3.3 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
       //7.3.4 calculate pixel's cross sectional area saturated part
    Acs_out = (depth - wtd) * sid_oct * cellout;
    // transmitance over the unsaturated part of the profile [mm^2/h]
    double T_uns = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_m), (3.0*lambda+1.0))-pow((bub_press/(psi_m +(wtd*1000.0))), (3.0*lambda+1.0)) );
    // adjust transmitance to m3/day
    double Q_uns = T_uns * ((24.0 * cellout * sid_oct)/(1.0e6));
    // adjust to flux density  and daily timestep [mm/day]
         
    T_uns *= ((24.0 * cellout*sid_oct)/(1000.0* Ai));
     // failsafe for low water contents where there is no saturated section or big storms
    if (T_uns  < 0.0 || isnan(T_uns)==1){
        T_uns  = 0.0;
        Q_uns  = 0.0;
    }
    // transmitance over the saturated part of the profile adjusted to flux density [mm/day]
    double T_sat = Ksat_visc*24.0* (Acs_out/Ai);
     // adjust transmitance to m3/day
    double Q_sat = Ksat_visc*24.0* (Acs_out)/1000;
    // total transmitance at initial soil moisture wn
    double T = (T_sat+T_uns)*hyd_grad;
     // total Q m3/day
    //double Q = (Q_sat+Q_uns)*hyd_grad;
    double Q = (T*Ai)/1000;
   
     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.5. recalculate volume Au from ti-1 if R<0.0 and correct td
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // double vol_prev = ((qin*Ai*pow(Kb,td))/log(Kb))-((qin*Ai)/log(Kb)) ;
    // if (R < 0.0){
    //     td = -1.0*log(1.0-(log(Kb)*((vol_prev-Au*R)/Qo)))/log(Kb);
    // }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.5. Estimate time required to drain upslope recharge [days]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // TO DO: recalculate volume Au from ti-1 if R<0.0 and correct tdrain
    double t_drain = 0.0;

    if((R > 0.0) && (sm > Wmax)){
        t_drain = -1.0*log(1.0-(log(Kb)*(Au*R/Q)))/log(Kb);
        q_in = Q/(Ai*pow(Kb,t_drain));
    } else {
        t_drain = (td -1.0);
        q_in = 0.0;
    }

    if( isnan(t_drain)==1 || (t_drain < 0.0)){
       t_drain = 0.0;
    } else if (t_drain >365.0){
        t_drain = 365.0;
    }
        

    sm += (q_in-T);


   /*  double t_drain = 0.0;
    if( sm >= Wmax){
       t_drain = log(((Au * R * log(Kb))/(Q*1000.0))+1)/log(Kb); 
    } else{
       t_drain = 0.0; 
    }
   

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.7. Update soil moisture
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // time to drain the upslope volume takes less than a day, so q_in happens the same day, only posible if R>0
    if((t_drain < 1.0) && (R > 0.0) && (t_drain > 0.0)){
       
       sm += (R*(Au/Ai)-T);
    } else if ((R > 0.0) && (t_drain > 1.0)) {
       sm += (To-T); 

    }
    // reset the counter if q_in happens the same day
   if((t_drain < 1.0) || isnan(t_drain)==1 || (R < 0.0)){
       t_drain = 0.0;
    }
   
    t_drain += (td -1.0);

    if((t_drain < 0.0) ){
       t_drain = 0.0;
    } */
       
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.7. Update soil moisture
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //sm -= (T);
    

    // failsafe for low water contents or big storms
    if (sm > SAT){
        sm = SAT;
    } else if (sm < RES) {
        sm = RES;
    }       
    
    
   
    
    
    dsm.sm = sm;
    dsm.ro = ro;
    dsm.swe = snow;
    dsm.sqout = max(qin*Kb,q_in);
    dsm.tdr = t_drain;
    
   
}

void SPLASH::run_one_day(int n, int y, double wn, double sw_in, double tc,
                         double pn,smr &dsoil, double slop, double asp,double snow, double snowfall, vector <double> &soil_info,
                         double qin, double td) {
    /* ***********************************************************************
    Name:     SPLASH.run_one_day
    Input:    - int, day of year (n)
              - int, year (y)
              - double, previous day's soil moisture, mm (wn)
              - double, daily sunshine fraction (sf)
              - double, daily air temperature, deg C (tc)
              - double, daily precipitation, mm (pn)
    Output:   None.
    Features: Calculates daily soil moisture and runoff based on STASH
              methods and updates class variables accordingly.
    *********************************************************************** */
   
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 0. Set meteorological variable:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    precip = pn;

    //*********************************************************************** */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // . Preprocess
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //read soil hydro points in mm
    //saturation (water  content mm at 0 KPa)
    double SAT = soil_info[0];
    //wilting point (water content mm at 1500 KPa)
    double WP = soil_info[1];
    //field capacity (water content mm at 33KPa)
    double FC = soil_info[2];
    //saturated hydraulic conductivity(mm/hr)
    double Ksat = soil_info[3];
    //slope of the log-log curve soil moisture-matric potetial from brooks and Corey (1964)
    double lambda = soil_info[4];
    double depth = soil_info[5];
    //bubbling pressure/capillarity fringe (mm)
    double bub_press = soil_info[6];
    //residual water content, test as WP?
    double RES = soil_info[1];
    //double RES = soil_info[7];
    //upslope area
    double Au = soil_info[8];
    //Pixel Area
    double Ai = soil_info[9];
    //neighbourhood  cells
    double cellin = soil_info[10];
    double cellout = soil_info[11];
    double theta_s = SAT/(depth *1000.0);
    double theta_r = RES/(depth *1000.0);
    double theta_fc = FC/(depth *1000.0);
    double theta_wp = WP/(depth *1000.0);
    double theta_q0 = theta_wp + 0.001;
    double theta_qs = theta_s;
    double theta_i = (wn)/(depth*1000.0);
    // correct theta_i for NA error reaching boundary conditions
    if (theta_i>=theta_s){
        theta_i = theta_s - 0.001;
    } else if (theta_i<=theta_r){
        theta_i = theta_wp + 0.001;
    }
    double alph = 4.0 + 2.0*lambda;
    // 7.1. assume hydraulic gradient as the slope
    double hyd_grad = (dtan(slop));
    // assuming octogonal cells, ncells draining will have n octagon sides of length[m]:
    double sid_oct = sqrt(Ai/(2.0*(1+sqrt(2.0))));
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 00. Compute real FC
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //double coeff_A = exp(log(33.0) + (1.0/lambda)*log(FC/(depth *1000.0)));
    //double Wmax = pow((coeff_A/(depth*1000)), (1.0/((1/lambda)+1.0))) * (depth *1000.0) ;
    
    double psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // // 1. Calculate evaporative supply rate (sw), mm/h using calculate water in the first 50cm of depth
    // // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 1.0 get wtd
    double wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // 2.0 get unsaturated/saturated depth within root zone if wtd is shallow, MAX ROOTING DEPTH = 2.0 m
    /////////////////////////////////////////////////////////////////////////////////////
    double z_uns =0.0;
    double sat_swc = 0.0;
    if (wtd<=2.0){
        z_uns = wtd*1000.0;
        sat_swc = theta_s*(2.0-wtd)*1000.0;
    }else{
        z_uns = 2.0*1000.0;
        sat_swc = 0.0;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // 4.0 calculate water in the within root zone
    /////////////////////////////////////////////////////////////////////////////////////
    double w_uns_z = (theta_r*z_uns) + (((psi_m+z_uns)*(theta_r-theta_s)*pow((bub_press/(psi_m+z_uns)),lambda))/(lambda-1));
    double w_uns_0 = (theta_r*0.0) + (((psi_m+0.0)*(theta_r-theta_s)*pow((bub_press/(psi_m+0.0)),lambda))/(lambda-1));
	double w_uns = w_uns_z-w_uns_0;
    double w_z = w_uns + sat_swc ;
    /////////////////////////////////////////////////////////////////////////////////////
    // 3.0 calculate real FC
    /////////////////////////////////////////////////////////////////////////////////////
    double coeff_A = exp(log(33.0) + (1.0/lambda)*log(theta_fc));
    double Wmax = pow((coeff_A/(z_uns)), (1.0/((1/lambda)+1.0))) * (z_uns) ;
    /////////////////////////////////////////////////////////////////////////////////////
    // 5.0 calculate supply rate (sw)
    /////////////////////////////////////////////////////////////////////////////////////
    double sw = 0.0;
     if (depth>2.0){
        double RES_z = theta_r * z_uns;
        sw = Global::Cw*((w_z-RES_z)/(Wmax-RES_z));
    } else{
        Wmax = pow((coeff_A/(depth*1000)), (1.0/((1/lambda)+1.0))) * (depth *1000.0) ;
        sw = Global::Cw*((wn-RES)/(Wmax-RES));
    }
   
    //double sw = Global::Cw*((wn)/(SAT));
    if (sw < 0.0 || isnan(sw)==1) {
         sw = 0.0;
     } else if (sw > Global::Cw){
         sw = Global::Cw;
     }
      
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Accumulate snow (mm) first so feedback on albedo happens
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    snow += snowfall;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 2. Calculate radiation and evaporation quantities
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    evap.calculate_daily_fluxes(sw, n, y, sw_in, tc, slop, asp, snow);
    dvap = evap.get_vals();
    econ = dvap.econ;
    pw = dvap.pw;
    rn_d = dvap.rn_d;
    snow -= dvap.snowmelt;
    visc = dvap.visc;

    double melt_enrg = (dvap.snowmelt/1000.0)*pw*Global::kfus;
	// energy available after snowmelt
    double AE = rn_d - melt_enrg;
    double sublimation = min(dvap.snowmelt,(AE*econ)*1000.0);
    double snowmelt = dvap.snowmelt - sublimation;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate daily drainage, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // correct Ksat for fuid viscosity
	double int_perm = Ksat/Global::fluidity;

	double Ksat_visc = int_perm*((pw*Global::G)/visc)*3.6;
    //double Ksat_visc = Ksat;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4. Calculate today's soil moisture (sm), mm
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 4.1 calculate inflow
    double inflow = pn + 0.1*dvap.cond + snowmelt;
    // 4.2 calculate skin moisture at 5 cm
    
    double surf_moist = moist_surf(depth,5.0,bub_press,wn,SAT,RES,lambda);
    // 4.2 calculate infiltration assuming storm duration 6hrs
    double infi = inf_GA(bub_press,surf_moist,Ksat_visc,theta_s,lambda,inflow,6.0,slop);
    // 4.2 calculate Dunne runoff
    double ro_d = max(inflow-infi,0.0);
    // 4.2 calculate recharge
    double R = infi - dvap.aet;
   // 4.2 update soil moisture
    //double sm = wn + pn + dn.cond - dn.aet + snowmelt;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.4. Estimate recession constant Kb
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.4.1 estimate baseflow at ~wilting point
    // 7.4.1.1 calculate matric potential (mmH2O)
    double psi_q0 = bub_press/pow((((theta_q0-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.4.1.2 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    double wtd_q0 = ((bub_press-psi_q0)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd_q0 < 0.0 || isnan(wtd_q0)==1){
        wtd_q0 = 0.0;
    }else if (wtd_q0 > depth){
        wtd_q0 = depth;
    }
    //7.4.1.3 calculate pixel's cross sectional area for baseflow output
    double Acs_out_q0 = (depth - wtd_q0) * sid_oct * cellout;
    //7.4.1.4 calculate transmitance over the unsaturated part of th profile [mm^2/h]
    double T_q0 =  (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_q0), (3.0*lambda+1.0))-pow((bub_press/(psi_q0 +(wtd_q0*1000.0))), (3.0*lambda+1.0)) );
    //7.4.1.5 adjust transmitance to m3/day
    double Q_q0 = T_q0 * hyd_grad* ((24.0 * sid_oct * cellout)/(1.0e6));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 7.4.2 estimate baseflow at saturation Q_qs
    // 7.4.2.1 calculate matric potential (mmH2O)
    double psi_qs = bub_press/pow((((theta_qs-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.4.2.2 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    double wtd_qs = ((bub_press-psi_qs)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd_qs < 0.0 || isnan(wtd_qs)==1){
        wtd_qs = 0.0;
    }else if (wtd_qs > depth){
        wtd_qs = depth;
    }
    //7.4.2.3 calculate pixel's cross sectional area for baseflow output [m^2]
    double Acs_out_qs = (depth) * sid_oct * cellout;
    //7.4.2.4 calculate transmitance over the unsaturated part of th profile [mm^2/h]
    double T_qs = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_qs), (3.0*lambda+1.0))-pow((bub_press/(psi_qs +(wtd_qs*1000.0))), (3.0*lambda+1.0)) );
    //7.4.2.5 adjust transmitance to m3/day
   // double Q_qs = (T_qs *hyd_grad *((24.0 * sqrt(Ai))/(1.0e6)))+(hyd_grad * Ksat_visc*24.0* (Acs_out_qs)/1000);
    double Q_qs = (hyd_grad * Ksat_visc*24.0* (Acs_out_qs)/1000.0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // 7.4.2 compute Kb
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    double Kb = exp((Q_q0-Q_qs)/((SAT-WP)*(Ai/1000.0)));
   // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.2. Estimate Transmitance at the beggining of the day
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.2.1 calculate volumetric water content for the whole soil column
    //theta_i = (wn)/(depth*1000.0);
    // 7.2.2 calculate matric potential (mmH2O)
    psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.2.3 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
    //7.2.4 calculate pixel's cross sectional area saturated part
    double Acs_out = (depth - wtd) * cellin * sid_oct;
    // transmitance over the unsaturated part of th profile [mm^2/h]
    double To_uns = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_m), (3.0*lambda+1.0))-pow((bub_press/(psi_m +(wtd*1000.0))), (3.0*lambda+1.0)) );
    // adjust transmitance to m3/day
    double Qo_uns = To_uns * ((24.0 * cellin * sid_oct)/(1.0e6));

    // adjust to flux density  and daily timestep [mm/day]
    To_uns *= ((24.0 * cellin * sid_oct)/(1000.0 * Ai));
    
    // transmitance over the saturated part of the profile adjusted to flux density [mm/day]
    double To_sat = Ksat_visc*24.0* (Acs_out/Ai);
     // adjust transmitance to m3/day saturated part
    double Qo_sat = Ksat_visc*24.0* (Acs_out)/1000.0;
    // total transmitance at initial soil moisture wn
    double To = (To_sat+To_uns)*hyd_grad;
    // total Qo
     double Qo = (Qo_sat+Qo_uns)*hyd_grad;  
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.6. Estimate q_in from upslope (flux density) [mm/day]
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    double q_in = 0.0; 
      

    if( (td == 0.0) ){ 
      q_in = 0.0;
    } else {
        q_in = qin*Kb;
    }

    double sm = wn + R + q_in;
        

   // ~~~~~~~~~~~~~~~~~~~~~~~
    // 5. Calculate hortonian runoff, mm
    // ~~~~~~~~~~~~~~~~~~~~~~~
    double ro_h = 0.0;
    //declare some variables
    //drainage efficiency
    double e_dr =0.0;
    // Baseflow
    double q_out = 0.0;
    //initial water content
    //double theta_i = 0.0;
    //matric potential
    //double psi_m = 0.0;
    //thicknes of the saturated part of the soil column
    double thickness = 0.0;
    //cross- sectional area of the saturated part
    //double Acs_out = 0.0;
    //water table depth
    double dummy = wn + R;
     if (dummy > SAT){
        // Bucket is too full:
        //   allocate excess water to runoff
        
        ro_h = (dummy - SAT);
        dummy = SAT;
        R -= ro_h;
   
    } else if (dummy < RES){
        // Bucket is too empty:
        //   set soil moisture & runoff to zero
        ro_h = 0.0;
        dummy = RES;
    } 
    if (sm > SAT){
        // Bucket is too full:
        //   allocate excess water to runoff
        
        //ro_h = (sm - SAT);
        sm = SAT;
        //R -= ro_h;
   
    } else if (sm < RES){
        // Bucket is too empty:
        //   set soil moisture & runoff to zero
        //ro_h = 0.0;
        sm = RES;
    } 
    // 5.1 update runoff
    double ro  = ro_d + ro_h;
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 6. Calculate drainage vadose sone (unsaturated) (used only to compare against field measurements at low depths)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // double subro = Ksat_visc*pow((sm/SAT), (3+(2/lambda)))*24*(dsin(slop));
    // sm -= subro;
    // double q_drain = Ksat_visc*pow((sm/SAT), (3+(2/lambda)))*24*(dcos(slop));
    // sm -=q_drain;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7. Calculate transmitance previous day To
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.3. Estimate Transmitance after recharge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.3.1 calculate volumetric water content for the whole soil column
    theta_i = (sm)/(depth*1000.0);
    // 7.3.2 calculate matric potential (mmH2O)
    psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 7.3.3 calculate the thickness ot the saturated section of the soil column as: depth - wtd (m)
    wtd = ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
       //7.3.4 calculate pixel's cross sectional area saturated part
    Acs_out = (depth - wtd) * sid_oct * cellout;
    // transmitance over the unsaturated part of the profile [mm^2/h]
    double T_uns = (Ksat_visc*bub_press/(3.0*lambda+1.0)) * ( pow((bub_press/psi_m), (3.0*lambda+1.0))-pow((bub_press/(psi_m +(wtd*1000.0))), (3.0*lambda+1.0)) );
    // adjust transmitance to m3/day
    double Q_uns = T_uns * ((24.0 * cellout * sid_oct)/(1.0e6));
    // adjust to flux density  and daily timestep [mm/day]
         
    T_uns *= ((24.0 * cellout*sid_oct)/(1000.0 * Ai));
    // failsafe for big storms
    if (T_uns  < 0.0 || isnan(T_uns)==1){
        T_uns  = 0.0;
        Q_uns  = 0.0;
    }
    // transmitance over the saturated part of the profile adjusted to flux density [mm/day]
    double T_sat = Ksat_visc*24.0* (Acs_out/Ai);
     // adjust transmitance to m3/day
    double Q_sat = Ksat_visc*24.0* (Acs_out)/1000;
    // total transmitance at initial soil moisture wn
    double T = (T_sat+T_uns)*hyd_grad;
     // total Q
    //double Q = (Q_sat+Q_uns)*hyd_grad;
    //double Q = (Q_sat+Q_uns)*0.3;
    double Q = (T*Ai)/1000;
    
    
   // TO DO: recalculate volume Au from ti-1 if R<0.0 and correct tdrain
    double t_drain = 0.0;

    if((R > 0.0) && (sm > Wmax)){
        t_drain = -1.0*log(1.0-(log(Kb)*(Au*R/Q)))/log(Kb);
        q_in = Q/(Ai*pow(Kb,t_drain));
    } else {
        t_drain = (td -1.0);
        q_in = 0.0;
    }

    if( isnan(t_drain)==1 || (t_drain < 0.0)){
       t_drain = 0.0;
    } else if (t_drain >365.0){
        t_drain = 365.0;
    }
        

    sm += (q_in-T);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 7.7. Update soil moisture
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //sm -= (T);
    

    // failsafe for low water contents or big storms
    if (sm > SAT){
        sm = SAT;
    } else if (sm < RES) {
        sm = RES;
    }       
    
    
            

    dsoil.sm = sm;
    dsoil.ro = ro;
    dsoil.swe = snow;
    dsoil.bflow = T;
    dsoil.sqout = max(qin*Kb,q_in); 
    dsoil.tdr = t_drain;

       
}

// How to return Lists to R
// In .h and .cpp -- Make return type of the function as List (in each file, at the begining, put #include <Rcpp.h>)
// at the end of the function, return List::create(values...)

List SPLASH::spin_up(int n, int y, vector<double> &sw_in, vector <double> &tair, vector <double> &pn, 
                        double slop, double asp, vector <double> &snowfall,vector <double> &soil_info){
    /* ***********************************************************************
    Name:     SPLASH.spin_up
    Input:    DATA class (d)
    Output:   None
    Features: Spins up the daily soil moisture
    Depends:  quick_run
    *********************************************************************** */
    double wn;  // previous day's soil moisture, mm
    double snow;  // previous day's swe, mm
    double qin;  // previous day's drainage
    double td;  // days left to drain the area upslope
    smr dsm;    // daily soil moisture, swe and runoff
    //residual water content
    //read soil hydro points in mm
    //saturation (water  content mm at 0 KPa)
    double SAT = soil_info[0];
    //wilting point (water content mm at 1500 KPa)
    double WP = soil_info[1];
    //slope of the log-log curve soil moisture-matric potetial from brooks and Corey (1964)
    double lambda = soil_info[4];
    double depth = soil_info[5];
    //bubbling pressure/capillarity fringe (mm)
    double bub_press = soil_info[6];
    //residual water content, test as WP?
    double RES = soil_info[1];
    //double RES = soil_info[7];

    double theta_s = SAT/(depth *1000.0);
    double theta_r = RES/(depth *1000.0);
    

    // Prepare daily soil moisture vector
    //    int n = d.nlines();
    //    int y = d.get_year();
    // wn_vec.resize(n, RES);
    vector <double> wn_vec(n,RES);
    vector <double> ro_vec(n);
	vector <double> snow_vec(n,0.0);
    vector <double> tdrain_vec(n,0.0);
    vector <double> qin_prev_vec(n,0.0);
    // Run one year:
    for (int i=0; i<n; i++){
        // Get preceeding soil moisture status:
        if (i == 0){
            wn = wn_vec[n-1];
            snow = snow_vec[n-1];
            qin = qin_prev_vec[n-1];
            td = tdrain_vec[n-1];
        } else {
            wn = wn_vec[(i-1)];
            snow = snow_vec[(i-1)];
            qin = qin_prev_vec[(i-1)];
            td = tdrain_vec[(i-1)];
        }

        // Calculate soil moisture and runoff
        quick_run((i+1), y, wn, sw_in[i], tair[i],
                  pn[i], dsm, slop, asp, snow ,snowfall[i],soil_info, qin, td);

        wn_vec[i] = dsm.sm;
        ro_vec[i] = dsm.ro;
        snow_vec[i] = dsm.swe;
        tdrain_vec[i] = dsm.tdr;
        qin_prev_vec[i] = dsm.sqout;
        
    }

    // Calculate change in starting soil moisture:
    double start_sm = wn_vec[0];
    double start_wtd = get_wtd( start_sm,  depth,  bub_press, theta_s, theta_r, lambda);
    quick_run(1, y, wn_vec[n-1], sw_in[0], tair[0],
              pn[0], dsm, slop, asp,snow_vec[n-1],snowfall[0],soil_info,qin_prev_vec[n-1],tdrain_vec[n-1]);
    double end_sm = dsm.sm;
    double end_wtd = get_wtd(end_sm,  depth,  bub_press, theta_s, theta_r, lambda);
    double diff_sm = (end_sm - start_sm);
    double diff_swe = (dsm.swe - snow_vec[0]);
    double diff_wtd = (end_wtd - start_wtd);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
   
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //abs() works weird inside the while loop
    if (diff_sm < 0){
        diff_sm = (start_sm - end_sm);
    }
    if (diff_swe < 0){
        diff_swe = (snow_vec[0] - dsm.swe);
    }
    if (diff_wtd < 0){
        diff_wtd = (start_wtd - end_wtd);
    }
    // Equilibrate (diff_sm > 1.0) && (diff_swe > 1.0) && (spin_count < 100)
    int spin_count = 1;
    while ((diff_sm > 1.0) && (diff_swe > 1.0) ){
        for (int i=0; i<n; i++){
            // Get preceeding soil moisture status:
            if (i == 0){
                wn = wn_vec[n-1];
                snow = snow_vec[n-1];
                qin = qin_prev_vec[n-1];
                td = tdrain_vec[n-1];
            } else {
                wn = wn_vec[(i-1)];
                snow = snow_vec[(i-1)];
                td = tdrain_vec[(i-1)];
                qin = qin_prev_vec[(i-1)];
            }

            // Calculate soil moisture and runoff
            quick_run((i+1), y, wn, sw_in[i], tair[i],
                      pn[i], dsm, slop, asp, snow ,snowfall[i],soil_info,qin,td);

            wn_vec[i] = dsm.sm;
            ro_vec[i] = dsm.ro;
            snow_vec[i] = dsm.swe;
            tdrain_vec[i] = dsm.tdr;
            qin_prev_vec[i] = dsm.sqout;
        }

        // Calculate difference
        start_sm = wn_vec[0];
        start_wtd = get_wtd( start_sm,  depth,  bub_press, theta_s, theta_r, lambda);
        quick_run(1, y, wn_vec[n-1], sw_in[0], tair[0],
                  pn[0], dsm, slop, asp,snow_vec[n-1],snowfall[0],soil_info,qin_prev_vec[n-1],tdrain_vec[n-1]);
        end_sm = dsm.sm;
        end_wtd = get_wtd(end_sm,  depth,  bub_press, theta_s, theta_r, lambda);
        diff_sm = (end_sm - start_sm);
        diff_swe = (dsm.swe - snow_vec[0]);
        diff_wtd = (end_wtd - start_wtd);
        if (diff_sm < 0){
            diff_sm = (start_sm - end_sm);
        }
        if (diff_wtd < 0){
        diff_wtd = (start_wtd - end_wtd);
        }
        spin_count++;
    }

    // Save initial soil moisture condition:
    dsoil.sm = wn_vec[n-1];
    
    return List::create(Named("sm") = wn_vec, Named("snow") = snow_vec,Named("qin") = qin_prev_vec,Named("tdrain") = tdrain_vec, Named("ro") = ro_vec);
}

List SPLASH::run_one_year(int n, int y, vector<double> &sw_in, vector <double> &tair, vector <double> &pn, vector <double> &wn_vec, 
                            double slop, double asp, vector <double> &snow, vector <double> &snowfall, vector <double> &soil_info,
                            vector <double> &qin_vec,vector <double> &td_vec){
    /* ***********************************************************************
    Name:     SPLASH.spin_up
    Input:    DATA class (d)
    Output:   None
    Features: Spins up the daily soil moisture
    Depends:  quick_run
    *********************************************************************** */
    double wn;  // previous day's soil moisture, mm
    double swe;  // previous day's swe, mm
    double qin;  // previous day's drainage
    double td;  // days left to drain the area upslope
    smr dsoil;    // daily soil moisture and runoff
    
    // Prepare daily outputs vector
    int n_end = wn_vec.size();
    //    int y = d.get_year();
    vector <double> pet_vec(n,0.0);
	vector <double> aet_vec(n,0.0);
    vector <double> cond_vec(n,0.0);
    vector <double> ro_vec(n,0.0);
    vector <double> snow_vec(n,0.0);
    vector <double> bflow_vec(n,0.0);
    vector <double> sqout_vec(n,0.0);
    vector <double> netr_vec(n,0.0);
    vector <double> tdrain_vec(n,0.0);
    //vector <double> sm_vec(n,0.0);
    // Run one year:
    for (int i=0; i<n; i++){
        // Get preceeding soil moisture status:
        if (i == 0){
            wn = wn_vec[n_end-1];
            swe = snow[n_end-1];
            qin = qin_vec[n_end-1];
            td = td_vec[n_end-1];
        } else {
            wn = wn_vec[(i-1)];
            swe = snow_vec[(i-1)];
            qin = sqout_vec[(i-1)];
            td = tdrain_vec[(i-1)];
        }

        // Calculate soil moisture and runoff
        run_one_day((i+1), y, wn, sw_in[i], tair[i],pn[i], dsoil, slop, asp, swe, snowfall[i], soil_info, qin, td);

        wn_vec[i] = dsoil.sm;
        ro_vec[i] = dsoil.ro;
        snow_vec[i]=dsoil.swe;
        bflow_vec[i]=dsoil.bflow;
        sqout_vec[i]=dsoil.sqout;
        pet_vec[i] =dvap.pet;
        aet_vec[i] =dvap.aet;
        cond_vec[i] =dvap.cond;
        netr_vec[i] =dvap.rn_d/1e6;
        tdrain_vec[i] = dsoil.tdr;
    }
 
    dsoil.sm = wn_vec[n-1];
    
return List::create(Named("wn") = wn_vec, Named("ro") = ro_vec, Named("pet") = pet_vec, Named("aet") = aet_vec, 
                    Named("snow") = snow_vec, Named("cond") = cond_vec, Named("bflow") = bflow_vec,Named("netr") = netr_vec,
                    Named("tdrain") = tdrain_vec, Named("qin_prev") = sqout_vec);
}



double SPLASH::moist_surf(double depth,double z,double bub_p,double wn,double SAT,double RES,double lambda){
    /* ***********************************************************************
    Name:     SPLASH.moist_surf
    Input:    - double depth, depth of the soil column (m)
              - double z, depth of the superficial layer (cm)
              - double bub_p, bubbling pressure (capillarity fringe) (mmH2O)
              - double wn, water content (mm)
              - double SAT, water content at saturation(mm)
              - double RES, residual water content(mm)
              - double lambda, slope of the log-log curve soil moisture-Matric potential (pore distribution index)
    Output:   tetha_i
    Features: Estimates the volumetric water content in the soil most superficial layer, following Brooks and Corey (1964)
    *********************************************************************** */
    //1. get hydro points vol/vol,
    double theta_r = RES/(depth*1000.0);
    double theta_s = SAT/(depth*1000.0);
    //2. convert lengths to cm,
	//depth *= 100;
	//wn /= 10;
    double bubbling_pr = bub_p/10;
    //2. get average volumetric content,
	double theta_mean = (wn)/(depth*1000.0);
    //2. get matric potential (cm),
	double water_pot_BC = bubbling_pr/pow((((theta_mean-theta_r)/(theta_s-theta_r))),(1/lambda));
    //2. get total hydraulic head (cm),
    double total_head_BC = water_pot_BC + z;
	//2. volumetric w content at head x (cm), BC overstimate theta, so theta r is not added to theta
	double theta_BC = (theta_s-theta_r)* pow((total_head_BC/bubbling_pr),(-1*lambda)) + theta_r;
	//wtd<-bubbling_pr-water_pot_BC
    //failsafes, ptfs are sometimes rubbish
	if(theta_mean<theta_r){
		theta_BC = theta_r;
        } else if (isnan(theta_BC)==1){
		theta_BC = theta_s;	
		}
	
	return theta_BC;



}

double SPLASH::inf_GA(double bub_press,double theta_i,double Ksat,double theta_s,double lambda,double P,double td,double slop){
	// ************************************************************************
	// # Name:     inf_GA
	// # Input:    - double, bub_press, air entry or bubbling pressure (mm)
	// #           - double, theta_s, saturation(volumetric fraction)
	// #           - double, theta_i, initial volumetric moisture
	// #           - double, td, duration of the event (hr)
	// #           - double,Ksat, sarurated hydraulic conductivity (mm/hr)
    // #           - double,slop, slope
	
	// # Output:   double:
	// #           - infiltration accumulated (mm)
	// #           
	// # Features: calculate infiltration accumulated using greem-Amp, modified to consider slopes, and non iterative
	// # Ref:      Assouline, S., 2013. Infiltration into soils: Conceptual approaches and solutions. 
	// #           Water Resour. Res. 49, 1755–1772. doi:10.1002/wrcr.20155
	// # ************************************************************************
	// # testing
	
	// # # end testing
	// # rainfall rate r
	double r = P/td;
	
	//# efective suction at wetting front
	double h_f = ((2+3*lambda)/(1+3*lambda))*(bub_press/2);
	double delta_head = h_f;
	double delta_theta = (theta_s-theta_i);
    double I = 0.0;
    double tp= 0.1;
    double tp_s =0.1;
	if(r <= Ksat){
        
		I = P;
	}else{
		if(delta_theta<=0.0){
			// if soil is at saturation, infiltration rate = ksat
            tp= 0.0;
			I = Ksat*td;
		}else{
            //time ponding
			tp = delta_theta*delta_head*((1/(r-Ksat))+(log((Ksat/(Ksat-r))+1)/Ksat));
            //tp = (Ksat*-1.0*delta_theta*delta_head)/(r*(r-Ksat));
			if(isnan(tp)==1){
				// massive storms, (Ksat/(Ksat-r))+1) result in negative
                tp=0.01;
			}
			// correction of the ponding time for slopes
			tp_s = tp/pow(dcos(slop),2.0);
			I = r*tp_s+(Ksat*(td-tp_s)-(delta_head*delta_theta*log(1-(r*tp_s/(delta_head*delta_theta)))));
		}
		
	}

	if(I>P){
        // errors with the log on big storms
		I = P;
	}
		
	return I;
}

double SPLASH::get_wtd(double wn, double depth, double bub_press,double theta_s,double theta_r,double lambda){
    /* ***********************************************************************
    Name:     SPLASH::get_wtd
    Input:    double wn, double depth, double bub_press,double theta_s,double theta_r,double lambda
    Output:   double
    Features: water table depth, meters.
    *********************************************************************** */
    double theta_i = (wn)/(depth*1000.0);
    // correct theta_i for NA error reaching boundary conditions
    if (theta_i>=theta_s){
        theta_i = theta_s - 0.001;
    } else if (theta_i<=theta_r){
        theta_i = theta_r + 0.001;
    }
    double psi_m = bub_press/pow((((theta_i-theta_r)/(theta_s-theta_r))),(1/lambda));
    // 1.0 get wtd
    double wtd= ((bub_press-psi_m)/1000.0);
    // failsafe for low water contents where there is no saturated section or big storms
    if (wtd < 0.0 || isnan(wtd)==1){
        wtd = 0.0;
    }else if (wtd > depth){
        wtd = depth;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return wtd;
}
double SPLASH::get_elv(){
    /* ***********************************************************************
    Name:     SPLASH.get_elv
    Input:    None
    Output:   double
    Features: Returns the elevation, meters.
    *********************************************************************** */
    return elv;
}

double SPLASH::get_lat(){
    /* ***********************************************************************
    Name:     SPLASH.get_lat
    Input:    None
    Output:   double
    Features: Returns the latitude, degrees.
    *********************************************************************** */
    return lat;
}

void SPLASH::print_vals(){
    /* ***********************************************************************
    Name:     SPLASH.print_vals
    Input:    None
    Output:   None
    Features: Prints the current dvap values.
    *********************************************************************** */
    evap.display();
    printf("Daily SPLASH values:\n");
    printf("  EET: %0.6f mm\n", dvap.eet);
    printf("  PET: %0.6f mm\n", dvap.pet);
    printf("  AET: %0.6f mm\n", dvap.aet);
    printf("  Cn: %0.6f mm\n", dvap.cond);
    printf("  Wn: %0.6f mm\n", dsoil.sm);
    printf("  RO: %0.6f mm\n", dsoil.ro);
}

void SPLASH::print_daily_wn(){
    /* ***********************************************************************
    Name:     SPLASH.print_daily_wn
    Input:    None
    Output:   None
    Features: Prints the current wn_vec values.
    *********************************************************************** */
    printf("%s\n", "Day,Wn (mm)");
    for (unsigned int i=0; i<wn_vec.size(); i++){
        printf("%d,%0.6f\n", i, wn_vec[i]);
    }
}
