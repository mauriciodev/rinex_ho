/*-------------------------------------------------------------------
    Purpose: Constants to be used in RINEX_HO
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1
    Date: July of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

#ifndef Constante
#define Constante
#include <math.h>

//CONSTANTS:
const double u = 3.986004418E+14;                                          //Gravitational cosntant in m/s2
const double c = 299792458.0;                                              //Velocity of light
const double we = 7.2921151467E-5;                                         //Earth rotation velocity
//const double we = 7.292115e-5;                                           //Earth rotation velocity
const double RE = 6371000.0;                                               //Earth Equatorial radius
const double hion = 450000.0;                                              //height of ionospheric layer
const double PI = 3.1415926535898;                                         //pi
const double flattwgs84 = 1.0/298.257222101;                               //WGS84 flattening
const double semjaxwgs84 = 6378137.0;                                      //WGS84 semi-major axis
const double exc_wgs84 = sqrt( (2.0*flattwgs84) - (pow(flattwgs84,2)) );   //excentricity
const double freq1 = 1575420000.0;                                         //L1 frequency in Hz
const double freq2 = 1227600000.0;                                         //L2 frequency in Hz

const double wlwl = c/(freq1-freq2);                                       //wide lane wavelenght

//iono free coefficients:
const double m1 = pow(freq1,2)/(pow(freq1,2)-pow(freq2,2));
const double m2_code = pow(freq2,2)/(pow(freq1,2)-pow(freq2,2));
const double m2_phase = (freq1*freq2)/( pow(freq1,2)-pow(freq2,2) );
const double m2_phase_tgd = (freq2/(freq1-freq2));

const double beta = (pow(freq1,2)/pow(freq2,2));

const double wlf1 = c/freq1;         //L1 wavelengh ~ 19 cm
const double wlf2 = c/freq2;         //L2 wavelengh ~ 24.4 cm

//for widelane Code and Phase
const double wl1r = freq1/(freq1+freq2);
const double wl2r = freq2/(freq1+freq2);
const double wl1p = wlf1*freq1/(freq1-freq2);
const double wl2p = -wlf2*freq2/(freq1-freq2);
const unsigned int NSAT = 33;
const int max_data = 300;

#endif