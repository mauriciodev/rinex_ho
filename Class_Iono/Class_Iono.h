/*-------------------------------------------------------------------
    Purpose:  C++ Class to compute second and third order ionospheric effects
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation: 
    -------------------------------------------------------------------*/

//------------------------------------------------------------------------------------

#ifndef CLASSIONO_H
#define CLASSIONO_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <ios>


//#include ".\IGRF11\IGRF11.h"
//#include "fortran.h"  //Class to mix C/C++ and Fortran (available at: http://arnholm.org/software/index.htm)

#include "../IGRF11/Class_IGRF11.h"

using namespace std;
using namespace IGRF11;

//for linux compiling C/C++ and fortran together
extern "C" void geodcgeomag_ (int *day, float *ut, float *lat, float *lon,
		                        float *mlat, float *mlon, float *mlt, int *ler);

//for Windows DLL
//extern "C" __declspec(dllimport) void _stdcall GEODCGEOMAG( int *DAY,float *UT,float *LAT,  float *LON,
//                                                            float *MLAT,float *MLON,float *MLT,int *ler);
namespace CIONO{
class Class_Iono : public Class_IGRF11{


    public:

           double zl,Elevm,Azm, //zenith, elevation and azimuth angle (geomagnetic system)
                  Az,El,Zen,    //zenith, elevation and azimuth angle(Terrestrial system)
                  var_I2_L1,
                  var_I2_L2;

           //Class iono constructor
           Class_Iono( double phi_pn_mag = 79.8,    //magnetic north pole latitude
                 double long_pn_mag = 288.1,  //magnetic north pole longitude
                 double h_ion = 450000.0,     //ionospheric layer height (m)
                 double R_e = 6371000.0,      //Equatorial Earth axis
                 double B_eq = 3.12E-05,      //Geomagnetic induction magnitude at equator (Tesla)
                 double Ne_max_in = 3.0e12);  //Electrons Density Maximum

           ~Class_Iono();//destructor

           //transform to geomagnetic e compute ionospheric effects
           void Calc_Iono( double TEC,int read_tec, int MJD,double fracOfDay,
                           double Xs,double Ys,double Zs,
                           double Xest,double Yest,double Zest,
                           double lat,double lamb,double h,
                           const int Geomag_Type,char Name_IGRF_Model[]);

           //Cartesianas para Curvilinear coordinates
           void Cart_Curv( double x,double y,double z,double &fi,double &lamb,
                           double &h);

           //Curvilinear coordinates to cartesians
           void Curv_Cart( double lat,double longi,double h,double &X,
                           double &Y, double &Z);

           //geodetic coordinates of the pierce point
           void Pierce_Point_Coord( double Xest,double Yest,double Zest,
                                    double Xsat,double Ysat,double Zsat,
                                    double lat,double lamb, double h);
           //Get Methods
           double Get_I1_L1(); //first order iono
           double Get_I1_L2();
           double Get_I2_L1(); //Second order iono
           double Get_I2_L2();
           double Get_I3_L1(); //Third order iono
           double Get_I3_L2();

           int Get_year();
           int Get_month();
           int Get_day();

           int Get_dayofyear();
           int Get_hour();
           int Get_min();
           double Get_sec();

           //ellipsoide parameters
           double Get_Major_Semi_Axis(); //Semi-major axis
           double Get_Minor_Semi_Axis(); //semi-minor axis
           double Get_f(); //flattenig
           double Get_exc(); //excentricity

           //magnetic north pole coordenates
           double Get_lat_pn_mag();
           double Get_lamb_pn_mag();
           double Get_lat_ps_mag();
           double Get_lamb_ps_mag();

           double Get_Re();   
           double Get_hion(); 
           double Get_Beq(); 
           double Get_f1();
           double Get_f2();
           double Get_lamb_l1();
           double Get_lamb_l2();   
           double Get_c();
           double Get_Ne_max();

           //Pierce point coordinates
           double Get_lat_ion();
           double Get_lamb_ion(); 

           double Get_PI();

           //Set Methods
           void Set_year(int value);
           void Set_month(int value);
           void Set_day(int value);
           void Set_dayofyear(int value);
           void Set_hour(int value);
           void Set_min(int value);
           void Set_sec(double value);

           void Set_I1_L1(double value); //first order iono
           void Set_I1_L2(double value);
           void Set_I2_L1(double value); //Second order iono
           void Set_I2_L2(double value);
           void Set_I3_L1(double value); //Third order iono
           void Set_I3_L2(double value);

           //ellipsoid parameters
           void Set_Major_Semi_Axis(double value); 
           void Set_Minor_Semi_Axis(double value); 
           void Set_f(double value); 
           void Set_exc(double value); 

           void Set_lat_pn_mag(double value);
           void Set_lamb_pn_mag(double value);
           void Set_lat_ps_mag(double value);
           void Set_lamb_ps_mag(double value);

           void Set_Re(double value);   
           void Set_hion(double value); 
           void Set_Beq(double value); 
           void Set_f1(double value);
           void Set_f2(double value);
           void Set_lamb_l1(double value);
           void Set_lamb_l2(double value);
           void Set_c(double value);
           void Set_Ne_max(double value);


    private:

            void Local_Terrestrial_Vector(double Xest,double Yest,double Zest);

            void Project_Vector( double Xs,double Ys,double Zs,
                                 double Xest,double Yest,double Zest,
                                 double lat,double lamb,double h,
                                 double Nm,double Em,double Um,
                                 double &V_proj);

           void Ang_Local( double E, double N, double U,
                           double &Am, double &Em, double &Zm);
 
           
           double BtJ(double lat_ion_m,double z_m,double a_m);

           
           void Calcula_N(double lat,double &N);
          

           
           void Cart_Geomag(double Xest,double Yest, double Zest,
                            double &Xm,double &Ym,double &Zm);

           void cart_local (double xestacao,double yestacao,double zestacao,
                            double xorigem, double yorigem,double zorigem,
                            double &e,double &n,double &u,double a,double f);

           
           void CartEsfeGeoc_CartGeod(double XEsf,double YEsf,double ZEsf,double quisi,
                                      double &XGeod,double &YGeod,double &ZGeod);                     

          
           void Curv_EsfeGeoc_Cart(double lat_geoc,double long_geoc,double R,
                                   double &XEsf,double &YEsf,double &ZEsf);

           void Geod_Dipolo(double lat,double lamb,double &latdip,double &lambdip);

           //First order ionospheric effect
           void Iono_First_Order(double TEC);  //Not in use

           //Second order ionospheric effect
           void Iono_Second_Order(double TEC,double BtJ);

           //Third order ionospheric effect
           void Iono_Third_Order(double TEC);

           //Geodetic latitude in geocentric and vice-versa
           void Lat_Geod_Geoc(double &lat_geod,double &lat_geoc,int op);

           //transform to local system
           void Sist_Local(double xest,double yest,double zest,
                           double xsat,double ysat,double zsat,
                           double fi, double lambda,
                           double &e,double &n,double &u);

          double I1_L1,I1_L2, //first order iono
                 I2_L1,I2_L2, //Second order iono
                 I3_L1,I3_L2; //Third order iono

          double lat_ion,lamb_ion; //Pierce point coordinates

          //ellipsoid parameters
          double a, //semi-major axis
                 b, //semi-minor axis
                 f, //flattening
                 exc, //excentricity

                 //north pole magnetic
                 lat_pn_mag,lamb_pn_mag,
                 lat_ps_mag,lamb_ps_mag,
                 Re,    //Equatorial Earth axis
                 hion,  //ionospheric layer height (m)
                 Beq,   //Geomagnetic induction magnitude at equator (Tesla)
                 f1,f2,lamb_l1,lamb_l2,   //frequencies and wavelenghts
                 c,Ne_max;

           double pi;

           int  year,
                month,
                day,
                dayofyear,
                hour,
                min;

           double sec;

           double Et[3],Nt[3],Ut[3];
           bool Compute_Terr_Local;

};

}//namespace CIONO

#endif
//------------------------------------------------------------------------------------

