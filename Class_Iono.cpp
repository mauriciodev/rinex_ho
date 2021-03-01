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

//c::1304.13, SAH, comment out geodcgeomag() everywhere.



#include "Class_Iono.h"

namespace CIONO {
//------------------------------------------------------------------------------------
Class_Iono::Class_Iono(double phi_pn_mag,double long_pn_mag, double h_ion,double R_e,
           double B_eq,double Ne_max_in)
{ /*-------------------------------------------------------------------
    Purpose: Class constructor
    -------------------------------------------------------------------
    Input: phi_pn_mag - magnetic north pole latitude
           long_pn_mag -magnetic north pole longitude
           h_ion - ionospheric layer height (m)
           R_e - Equatorial Earth axis
           B_eq - Geomagnetic induction magnitude at equator
           Ne_max_in - Electrons Density Maximum

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    pi = 3.1415926535897932384626433832795;
    Set_f1(1575420000.0);                     //L1 - frequency (Hz)
    Set_f2(1227600000.0);                     //L2 - frequency (Hz)
    Set_f5(1176450000.0);                     //L2 - frequency (Hz)

    Set_c(299792458.0);                       //Light velocity
    Set_lamb_l1(c/f1);                        //wavelanght - L1
    Set_lamb_l2(c/f2);                        //wavelenght - L2
    Set_lamb_l5(c/f5);                        //wavelenght - L2

    //GRS80 parameters
    Set_Major_Semi_Axis(6378137.0);           //semi-major axis
    Set_f(1.0/298.257222101);                 //flattening
    Set_Minor_Semi_Axis(a - a*f);             //semi-minor axis
    Set_exc(sqrt( (2.0*f) - (pow(f,2)) ));    //excentricity

    //initial parameters - To be updated later
    Set_lat_pn_mag(phi_pn_mag*pi/180.0);      //magnetic north pole latitude
    Set_lamb_pn_mag(long_pn_mag*pi/180.0);    //magnetic north pole longitude

    Set_hion(h_ion);                          //ionospheric layer height (m)
    Set_Re(R_e);                              //Equatorial Earth axis

    Set_Ne_max(Ne_max_in);                    //Electrons Density Maximum

    Set_Beq(B_eq);                            //Geomagnetic induction magnitude at equator

    Set_year(0);
    Set_dayofyear(0);
    Set_hour(0);
    Set_min(0);
    Set_sec(0);


    for(int i=0;i<3;i++)

    {
        Et[i] = 0.0;
        Nt[i] = 0.0;
        Ut[i] = 0.0;
    }//for

    Compute_Terr_Local = false;

}
//------------------------------------------------------------------------------------
Class_Iono::~Class_Iono()//destrutor
{ /*-------------------------------------------------------------------
    Purpose: Class Destrucotr
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
   //Insert code here

}



//---- added by SH

void Class_Iono::azelIPP(double staX, double staY, double staZ,
     double svX, double svY, double svZ,
     double Re, double htIono,
     double *az, double *elv, double *dist,
     double *pplat, double *pplon, double *zp)
{
//   Input:
//   staX,staY,staZ -- site coordinates in meters
//   svX,svY,svZ -- satellite coordinates in meters
//   Re, htIono -- Mean Radius of the Earth (6371 km) & height on Iono shell (450 km).
//   Output:
//   az, elv, dist -- azimuth, elevation in radians (dist in meters)
//   pplat, pplon -- geocentric lat.long of IPP (in degrees)
//   zp -- z' = zenith angle at the IPP (in radians).
//
//   Author: S. Hilla,  13 Feb 2013
//

// const double DTR  = 0.0174532925199433;
   double rlat, rlon, srlon, crlon, srlat, crlat, dx, dy, dz, du, dv, dw, ppii;
   double dlat, dlon, elliph, angleA, c, b, dipp, up, vp, wp;
   double xipp, yipp, zipp, dequa, DTR;
   double geod_pp_lat, geod_pp_lon, geod_pp_hgt, mapfunc;
   int iJMJ;

   ppii = (atan(1.0)) * 4.0; // compute a double precision PI
   DTR = ppii/180.0;

   xyz_to_geod( 6378137.0, 298.257222101, staX, staY, staZ,
                &dlat, &dlon, &elliph);

   rlat = dlat * DTR;
   rlon = dlon * DTR;

//   rlat = rlat - ppii/2.0;
//   rlon = rlon - ppii;

   srlat = sin(rlat);
   crlat = cos(rlat);
   srlon = sin(rlon);
   crlon = cos(rlon);

   dx = svX - staX;
   dy = svY - staY;
   dz = svZ - staZ;

// du = crlat * crlon * dx + crlat * srlon * dy - srlat * dz;
// dv = srlon * dx - crlon * dy;
// dw = srlat * crlon * dx + srlat * srlon * dy + crlat * dz;

   du = -1.0*srlat*crlon*dx - srlat*srlon*dy + crlat*dz;
   dv = -1.0*srlon*dx + crlon*dy;
   dw = crlat*crlon*dx + crlat*srlon*dy + srlat*dz;

   *dist = sqrt( du*du + dv*dv + dw*dw );

//cout << "du,dv,dw,dist: "
//<< setw(25) << setprecision(5) << du << " "
//<< setw(25) << setprecision(5) << dv << " "
//<< setw(25) << setprecision(5) << dw << " "
//<< setw(25) << setprecision(5) << *dist << endl;

   *az = atan2(dv,du);
   *elv = asin( (dw)/(*dist) );

   *zp = asin( (Re/(Re+htIono))*sin(ppii/2.0 - *elv) );

   // compute the distance from receiver to IPP
   angleA = ppii/2.0 - *zp - *elv;
   c = Re;
   b = Re + htIono;
   dipp = sqrt( b*b + c*c - 2.0*b*c*cos(angleA) );

   // compute uvw for IPP
   up = dipp * cos(*elv) * cos(*az);
   vp = dipp * cos(*elv) * sin(*az);
   wp = dipp * sin(*elv);

   // convert uvw for IPP to dx,dy,dz (reverse of 4.68-4.70, Geometric Geodesy).
   dx = -1.0*srlat*crlon*up - srlon*vp + crlat*crlon*wp;
   dy = -1.0*srlat*srlon*up + crlon*vp + crlat*srlon*wp;
   dz = crlat*up + srlat*wp;

   // compute XYZ for IPP
   xipp = staX + dx;
   yipp = staY + dy;
   zipp = staZ + dz;

   // compute GEODETIC latitude and longitude for IPP
   xyz_to_geod( 6378137.0, 298.257222101, xipp, yipp, zipp,
                &geod_pp_lat, &geod_pp_lon, &geod_pp_hgt);


   //cout<<"geod_pp_lat: "<<geod_pp_lat<<" geod_pp_lon: "<<geod_pp_lon<<endl;
   //getchar();

   // compute geocentric latitude and longitude for IPP
   dequa = sqrt( xipp*xipp + yipp*yipp);
   *pplat = atan2( zipp, dequa);
   *pplon = atan2( yipp, xipp);

}


/*---------------------------------------------------------------------*/

int Class_Iono::xyz_to_geod ( double a, double finv, double x, double y, double z,
                   double *dlat, double *dlon, double *h)

   /********************************************************************|
   | void xyz_to_geod() converts geocentric Cartesian coordinates into
   |                    geodetic latitude, longitude, and ellipsoid
   |                    height ( in decimal degrees and ?????? ) given
   |                    the semi-major axis "a" and the inverse of the
   |                    flattening "finv".
   |
   |                    The units of h will be the same as the units of
   |                    x,y,z, and a (m, km, miles, feet, etc.).
   |
   |   Original Fortran subroutine TOGEOD by C. Goad, 1987.
   |   Converted to C function by S. Hilla, July 1992.
   |********************************************************************/
 {
   double rtd, esq, oneesq, psq, p, r;
   double tolsq = 1.0e-10, sinlat, coslat, N, dp, dz;
   int  maxit = 10, iter;

   /* compute radians-to-degrees factor */
       rtd = 180.0 / PI;                     /* PI is in consts.h */
   /* compute square of eccentricity */
       esq = (2.0 - 1.0/finv)/finv;
       oneesq = 1.0 - esq;
   /* direct calculation of longitude */
       *dlon = ( atan2(y,x) ) * rtd;
       if ( *dlon < 0.0 ) *dlon = *dlon + 360.0;
   /* first guesses */
   /* p is distance from the spin axis */
       psq = (x*x) + (y*y);
       p = sqrt(psq);
   /* r is distance from origin (0,0,0) */
       r = sqrt( psq + (z*z) );
       sinlat = z/r;
       *dlat = asin(sinlat);
   /* initial value of height = distance from origin minus approximate */
   /* distance from origin to surface of ellipsoid */
       *h = r - a*(1.0 - ( (sinlat*sinlat)/finv ) );
   /* iterate */
       for ( iter=1; iter <= maxit; ++iter )
       {
         sinlat = sin(*dlat);
         coslat = cos(*dlat);
      /* compute radius of curvature in prime vertical direction */
         N = a / sqrt( 1.0 - esq*sinlat*sinlat );
      /* compute residuals in p and z */
         dp = p - ( (N + *h)*coslat );
         dz = z - (N*oneesq + *h)*sinlat;
      /* update height and latitude */
         *h = *h + (sinlat*dz + coslat*dp);
         *dlat = *dlat + (coslat*dz - sinlat*dp)/(N + *h);
      /* test for convergence */
         if (  ( dp*dp + dz*dz ) < tolsq )
          {
            *dlat = *dlat * rtd;
            return(0);
          }
       }       /* end of for loop */

   /* Not Converged -- Warn user */
       cout << " Problem in xyz_to_geod, did not " <<
       "converge in " << maxit << " iterations." << endl;
       return (1);

 }  /*  end of function xyz_to_geod  */


/*---------------------------------------------------------------------*/

void Class_Iono::xyzIPP(double staX, double staY, double staZ,
     double svX, double svY, double svZ,
     double Re, double htIono,
     double *xipp, double *yipp, double *zipp)
{
//   Input:
//   staX,staY,staZ -- site coordinates in meters
//   svX,svY,svZ -- satellite coordinates in meters
//   Re, htIono -- Mean Radius of the Earth (6371 km) & height on Iono shell (450 km).
//   Output:
//   *xipp,*yipp,*zipp --  X,Y,Z coordinates of the IPP (in meters).
//
//   Author: S. Hilla,  13 Feb 2013
//

// const double DTR  = 0.0174532925199433;
   double rlat, rlon, srlon, crlon, srlat, crlat, dx, dy, dz, du, dv, dw, ppii;
   double dlat, dlon, elliph, angleA, c, b, dipp, up, vp, wp;
   double dist, az, elv, zp, dequa, DTR;
   double geod_pp_lat, geod_pp_lon, geod_pp_hgt, mapfunc;
   int iJMJ;

   ppii = (atan(1.0)) * 4.0; // compute a double precision PI
   DTR = ppii/180.0;

   xyz_to_geod( 6378137.0, 298.257222101, staX, staY, staZ,
                &dlat, &dlon, &elliph);

   rlat = dlat * DTR;
   rlon = dlon * DTR;

//   rlat = rlat - ppii/2.0;
//   rlon = rlon - ppii;

   srlat = sin(rlat);
   crlat = cos(rlat);
   srlon = sin(rlon);
   crlon = cos(rlon);

   dx = svX - staX;
   dy = svY - staY;
   dz = svZ - staZ;

// du = crlat * crlon * dx + crlat * srlon * dy - srlat * dz;
// dv = srlon * dx - crlon * dy;
// dw = srlat * crlon * dx + srlat * srlon * dy + crlat * dz;

   du = -1.0*srlat*crlon*dx - srlat*srlon*dy + crlat*dz;
   dv = -1.0*srlon*dx + crlon*dy;
   dw = crlat*crlon*dx + crlat*srlon*dy + srlat*dz;

   dist = sqrt( du*du + dv*dv + dw*dw );

//cout << "du,dv,dw,dist: "
//<< setw(25) << setprecision(5) << du << " "
//<< setw(25) << setprecision(5) << dv << " "
//<< setw(25) << setprecision(5) << dw << " "
//<< setw(25) << setprecision(5) << dist << endl;

   az = atan2(dv,du);
   elv = asin( (dw)/(dist) );

   zp = asin( (Re/(Re+htIono))*sin(ppii/2.0 - elv) );

   // compute the distance from receiver to IPP
   angleA = ppii/2.0 - zp - elv;
   c = Re;
   b = Re + htIono;
   dipp = sqrt( b*b + c*c - 2.0*b*c*cos(angleA) );

   // compute uvw for IPP
   up = dipp * cos(elv) * cos(az);
   vp = dipp * cos(elv) * sin(az);
   wp = dipp * sin(elv);

   // convert uvw for IPP to dx,dy,dz (reverse of 4.68-4.70, Geometric Geodesy).
   dx = -1.0*srlat*crlon*up - srlon*vp + crlat*crlon*wp;
   dy = -1.0*srlat*srlon*up + crlon*vp + crlat*srlon*wp;
   dz = crlat*up + srlat*wp;

   // compute XYZ for IPP
   *xipp = staX + dx;
   *yipp = staY + dy;
   *zipp = staZ + dz;

} // end of function xyzIPP()


//---- above code added by SH



//------------------------------------------------------------------------------------
int Class_Iono::Get_year()
{/*-------------------------------------------------------------------
    Purpose: Return year
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(year);
}
//------------------------------------------------------------------------------------
int Class_Iono::Get_month()
{/*-------------------------------------------------------------------
    Purpose: Return month
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(month);
}
//------------------------------------------------------------------------------------
int Class_Iono::Get_day()
{/*-------------------------------------------------------------------
    Purpose: Return day
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(day);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I1_L1()
{/*-------------------------------------------------------------------
    Purpose: Return first order ionospheric effect - L1
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
//first order iono - L1
    return (I1_L1);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I1_L2()
{/*-------------------------------------------------------------------
    Purpose: Return first order ionospheric effect - L2
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
//first order iono - L2
    return (I1_L2);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I1_L5()
{/*-------------------------------------------------------------------
    Purpose: Return first order ionospheric effect - L5
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //first order iono - L5
    return (I1_L5);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I2_L1()
{/*-------------------------------------------------------------------
    Purpose: Return second order ionospheric effect - L1
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Second order iono - L1
    return (I2_L1);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I2_L2()
{/*-------------------------------------------------------------------
    Purpose: Return second order ionospheric effect - L2
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Second order iono - L2
    return (I2_L2);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I2_L5()
{/*-------------------------------------------------------------------
    Purpose: Return second order ionospheric effect - L5
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Second order iono - L5
    return (I2_L5);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I3_L1()
{/*-------------------------------------------------------------------
    Purpose: Return third order ionospheric effect - L1
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Third order iono - L1
    return (I3_L1);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I3_L2()
{/*-------------------------------------------------------------------
    Purpose: Return third order ionospheric effect - L2
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Third order iono - L2
    return (I3_L2);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_I3_L5()
{/*-------------------------------------------------------------------
    Purpose: Return third order ionospheric effect - L5
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
     //Third order iono - L5
    return (I3_L5);
}
//------------------------------------------------------------------------------------
int Class_Iono::Get_dayofyear()
{/*-------------------------------------------------------------------
    Purpose: Return day of the year
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(dayofyear);
}
//------------------------------------------------------------------------------------
int Class_Iono::Get_hour()
{/*-------------------------------------------------------------------
    Purpose: Return hours
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(hour);
}
//------------------------------------------------------------------------------------
int Class_Iono::Get_min()
{/*-------------------------------------------------------------------
    Purpose: Return minutes
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(min);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_sec()
{/*-------------------------------------------------------------------
    Purpose: Return seconds
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(sec);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_PI()
{/*-------------------------------------------------------------------
    Purpose: Return Pi
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    return(pi);
}
//------------------------------------------------------------------------------------
//parametros do elipsoide
double Class_Iono::Get_Major_Semi_Axis()
{/*-------------------------------------------------------------------
    Purpose: Return Semi-Major axis
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(a);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_Minor_Semi_Axis()
{/*-------------------------------------------------------------------
    Purpose: Return semi-minor axis
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(b);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_f()
{/*-------------------------------------------------------------------
    Purpose: Return flattening
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(f);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_exc()
{/*-------------------------------------------------------------------
    Purpose: Return Excentricity
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(exc);
}
//------------------------------------------------------------------------------------
//coordenadas dos polos magneticos
double Class_Iono::Get_lat_pn_mag()
{/*-------------------------------------------------------------------
    Purpose: Return latitude of magnetic north pole
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lat_pn_mag);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_pn_mag()
{/*-------------------------------------------------------------------
    Purpose: Return longitude of magnetic north pole
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lamb_pn_mag);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lat_ps_mag()
{/*-------------------------------------------------------------------
    Purpose: Return latitude of magnetic south pole
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lat_ps_mag);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_ps_mag()
{/*-------------------------------------------------------------------
    Purpose: Return longitude of magnetic south pole
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lamb_ps_mag);
}
//------------------------------------------------------------------------------------

double Class_Iono::Get_Re()
{/*-------------------------------------------------------------------
    Purpose: Equatorial earth axis
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(Re);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_hion()
{/*-------------------------------------------------------------------
    Purpose: Ionospheric layer height
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(hion);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_Beq()
{/*-------------------------------------------------------------------
    Purpose: Return Geomagnetic induction magnitude at Equator (Tesla)
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(Beq);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_f1()
{/*-------------------------------------------------------------------
    Purpose: Return L1 frequency
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(f1);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_f2()
{/*-------------------------------------------------------------------
    Purpose: Return L2 frequency
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(f2);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_f5()
{/*-------------------------------------------------------------------
    Purpose: Return L5 frequency
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(f5);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_l1()
{/*-------------------------------------------------------------------
    Purpose: Return L1 wavelenght
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lamb_l1);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_l2()   //frequencias e comprimentos de ondas
{/*-------------------------------------------------------------------
    Purpose: Return L2 wavelenght
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lamb_l2);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_l5()   //frequencias e comprimentos de ondas
{/*-------------------------------------------------------------------
    Purpose: Return L2 wavelenght
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(lamb_l5);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_c()
{/*-------------------------------------------------------------------
    Purpose: Return light velocity
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(c);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_Ne_max()
{/*-------------------------------------------------------------------
    Purpose: Return electrons density maximun
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   return(Ne_max);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lat_ion()
{/*-------------------------------------------------------------------
    Purpose: Return geodetic pierce point latitude
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    return(lat_ion);
}
//------------------------------------------------------------------------------------
double Class_Iono::Get_lamb_ion()
{/*-------------------------------------------------------------------
    Purpose: Return geodetic pierce point longitude
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    return(lamb_ion);
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_year(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    year = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_month(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    month = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_day(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    day = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_dayofyear(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    dayofyear = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_hour(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    hour = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_min(int value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    min = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_sec(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    sec = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I1_L1(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
//first order iono - L1
    I1_L1 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I1_L2(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //first order iono - L2
    I1_L2 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I1_L5(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    //first order iono - L5
    I1_L5 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I2_L1(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    //Second order iono - L1
    I2_L1 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I2_L2(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    //Second order iono - L2
    I2_L2 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I2_L5(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    //Second order iono - L5
    I2_L5 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I3_L1(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Third order iono - L1
    I3_L1 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I3_L2(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Third order iono - L2
    I3_L2 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_I3_L5(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    //Third order iono - L5
    I3_L5 = value;
}
//------------------------------------------------------------------------------------
//parametros do elipsoide
void Class_Iono::Set_Major_Semi_Axis(double value) //Semi-eico maior
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    a = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_Minor_Semi_Axis(double value) //semi-eixo menor
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    b = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_f(double value) //Achatamento
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    f = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_exc(double value) //excentriciddade
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    exc = value;
}
//------------------------------------------------------------------------------------
//coordenadas dos polos magneticos
void Class_Iono::Set_lat_pn_mag(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lat_pn_mag = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lamb_pn_mag(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lamb_pn_mag = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lat_ps_mag(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lat_ps_mag = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lamb_ps_mag(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lamb_ps_mag = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_Re(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    Re = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_hion(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    hion = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_Beq(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    Beq = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_f1(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    f1 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_f2(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    f2 = value;
}
//------------------------------------------------------------------------------------
 void Class_Iono::Set_f5(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    f5 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lamb_l1(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lamb_l1 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lamb_l2(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lamb_l2 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_lamb_l5(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: April of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    lamb_l5 = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_c(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    c = value;
}
//------------------------------------------------------------------------------------
void Class_Iono::Set_Ne_max(double value)
{/*-------------------------------------------------------------------
    Purpose:
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    Ne_max = value;
}
//------------------------------------------------------------------------------------
/*void iono::Lat_Geod_Esphe_Geoc(double lat,double h,double &lat_geoc)
{
    //WGS84
    a2 = 6378137^2;
    b2 = 6356752^2;

    cosfi = cos(lat);
    senfi = sin(lat);

    arg = (a2*cosfi^2 + b2*senfi^2);

    tanfi = ((sqrt(arg)*h + b2)/(sqrt(arg)*h + a2))*tan(lat);
    lat_geoc = atan(tanfi)

return;
} */
void Class_Iono::Lat_Geod_Geoc(double &lat_geod,double &lat_geoc,int op)
{/*-------------------------------------------------------------------
    Purpose: Transform geodetic latitude to geocentric latitude and vice-versa
    -------------------------------------------------------------------
    Input: lat_geod - geodetic latitude(WGS84) - radian
           lat_geoc - geocentric latitude - radianos
           op = 0 -> geodetic to geocentric
           op = 1 -> geocentric to geodetic

    Output: lat_geod - geodetic latitude(WGS84) - radian
            lat_geoc - geocentric latitude - radian
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    if(op==0)
    {   //geodetic to geocentric
        lat_geoc = atan( (1.0 - pow(exc,2))*tan(lat_geod));
    }
    else if(op==1)
    {   //geocentric para geodetic
        lat_geod = atan2( tan(lat_geoc),(1.0 - pow(exc,2) ) );
    }//else if

 return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Local_Terrestrial_Vector(double Xest,double Yest,double Zest)
{/*-------------------------------------------------------------------
     Purpose: To compute terrestrial local vector to be used in the
              geomagnetic propagation vector
     -------------------------------------------------------------------
     Input: Xest, Yest, Zest - Receiver coordinates - WGS84

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

double X_Norm[3]={0.0},Norm;
double Z_Axis[3] = {0, 0, 1.0};

     //Receiver vector Normalized
     Norm = sqrt(pow(Xest,2) + pow(Yest,2) + pow(Zest,2) );
     X_Norm[0] = Xest/Norm;
     X_Norm[1] = Yest/Norm;
     X_Norm[2] = Zest/Norm;

     //Radial Vector
     Ut[0] = X_Norm[0]*-1.0;
     Ut[1] = X_Norm[1]*-1.0;
     Ut[2] = X_Norm[2]*-1.0;

     //East direction vector
     Et[0] = Ut[1]*Z_Axis[2] - Ut[2]*Z_Axis[1];
     Et[1] = Ut[2]*Z_Axis[0] - Ut[0]*Z_Axis[2];
     Et[2] = Ut[0]*Z_Axis[1] - Ut[1]*Z_Axis[0];

     //East Normalized
     Norm = sqrt(pow(Et[0],2) + pow(Et[1],2) + pow(Et[2],2) );
     Et[0] = Et[0]/Norm;
     Et[1] = Et[1]/Norm;
     Et[2] = Et[2]/Norm;

     //North direction vector
     Nt[0] = Et[1]*Ut[2] - Et[2]*Ut[1];
     Nt[1] = Et[2]*Ut[0] - Et[0]*Ut[2];
     Nt[2] = Et[0]*Ut[1] - Et[1]*Ut[0];

     //North Normalized
     Norm = sqrt(pow(Nt[0],2) + pow(Nt[1],2) + pow(Nt[2],2));
     Nt[0] = Nt[0]/Norm;
     Nt[1] = Nt[1]/Norm;
     Nt[2] = Nt[2]/Norm;

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Project_Vector( double Xs,double Ys,double Zs,
                           double Xest,double Yest,double Zest,
                           double lat,double lamb,double h,
                           double Nm,double Em,double Um,
                           double &V_proj)
{/*-------------------------------------------------------------------
     Purpose: To compute the projection vector between geomagnetic and
              satellite-receiver vector in the pierce point
     -------------------------------------------------------------------
     Input: Xs, Ys, Zs - Satellite coordinates - WGS84
            Xest, Yest, Zest - Receiver coordiantes
            lat, lamb, h - geodetic latitude longitude and altitude
            Nm, Em, Um - geomagnetic coordinates from IGRF subroutine (in nT)

     Output: V_proj - projection vector between geomagnetic and
                      satellite-receiver vector in the pierce point
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              FAPESP PROCESS: 05/03522-1

     Date: July of 2010
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

double xyz[3]={0.0},dx[3]={0.0},Norm;

double Norm_dx;

double zxipp,zyipp,zzipp;

     //Vector in the satellite receiver direction
     dx[0] = Xs - Xest;
     dx[1] = Ys - Yest;
     dx[2] = Zs - Zest;

     Norm_dx = sqrt(pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2) );
     dx[0] = dx[0]/Norm_dx;
     dx[1] = dx[1]/Norm_dx;
     dx[2] = dx[2]/Norm_dx;
//cout << "#:PVec dx[0]: " << setw(30) << setprecision(15) << dx[0] << endl;
//cout << "#:PVec dx[1]: " << setw(30) << setprecision(15) << dx[1] << endl;
//cout << "#:PVec dx[2]: " << setw(30) << setprecision(15) << dx[2] << endl;

     xyzIPP(Xest,Yest,Zest,Xs,Ys,Zs,Re,hion,&zxipp,&zyipp,&zzipp);
     //Vector in the satellite receiver direction AT IPP
     dx[0] = zxipp - Xest;
     dx[1] = zyipp - Yest;
     dx[2] = zzipp - Zest;

     dx[0] = zxipp - Xs;
     dx[1] = zyipp - Ys;
     dx[2] = zzipp - Zs;

     Norm_dx = sqrt(pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2) );

     dx[0] = dx[0]/Norm_dx;
     dx[1] = dx[1]/Norm_dx;
     dx[2] = dx[2]/Norm_dx;
//cout << "#:PVec dx[0]ipp: " << setw(30) << setprecision(15) << dx[0] << endl;
//cout << "#:PVec dx[1]ipp: " << setw(30) << setprecision(15) << dx[1] << endl;
//cout << "#:PVec dx[2]ipp: " << setw(30) << setprecision(15) << dx[2] << endl;


      if(!Compute_Terr_Local)
      {
         Local_Terrestrial_Vector(Xest,Yest,Zest);
         Compute_Terr_Local = true;
      }

      
     //Local geomagnetic coordinates aligned with terrestrial system
     xyz[0] = Nm*Nt[0] + Em*Et[0] + Um*Ut[0];
     xyz[1] = Nm*Nt[1] + Em*Et[1] + Um*Ut[1];
     xyz[2] = Nm*Nt[2] + Em*Et[2] + Um*Ut[2];

/*cout << "#:PVec Nm: " << setw(30) << setprecision(15) << Nm << endl;
cout << "#:PVec Em: " << setw(30) << setprecision(15) << Em << endl;
cout << "#:PVec Um: " << setw(30) << setprecision(15) << Um << endl;
cout << "#:PVec xyz[0]: " << setw(30) << setprecision(15) << xyz[0] << endl;
cout << "#:PVec xyz[1]: " << setw(30) << setprecision(15) << xyz[1] << endl;
cout << "#:PVec xyz[2]: " << setw(30) << setprecision(15) << xyz[2] << endl;
*/

     //Projection vector between geomagnetic and satellite-receiver vector in the pierce point
     V_proj = (xyz[0]*dx[0] + xyz[1]*dx[1] + xyz[2]*dx[2])*1.0e-9;

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Calc_Iono(double TEC,int read_tec, int MJD,double fracOfDay,
                     double Xs,double Ys,double Zs,
                     double Xest,double Yest,double Zest,
                     double lat,double lamb,double h,
                     const int Geomag_Type,char Name_IGRF_Model[])
{/*-------------------------------------------------------------------
    Purpose: To compute ionospheric effects
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

 double lat_dip(0.0),lamb_dip(0.0),
        Xm(0.0),Ym(0.0),Zm(0.0),
        Em(0.0),Nm(0.0),Um(0.0),
        lat_ion_dip(0.0),lamb_ion_dip(0.0),h_dip(0.0),
        BTJ(0.0),
        Xsm(0.0),Ysm(0.0),Zsm(0.0),
        hs_dip(0.0),
        Xion,Yion,Zion,
        Xionm,Yionm,Zionm,hion_dip;

  double lats(0.0),lambs(0.0),
         lats_dip(0.0),lambs_dip(0.0);


  int I_ERROR;
  double late,longe;

  int ler;

  float MLT,ut,lats_temp,lambs_temp,lats_dip_temp,lambs_dip_temp;
  float lat_ion_dip_temp,lamb_ion_dip_temp;
  float lat_temp, lamb_temp,lat_dip_temp,lamb_dip_temp;
  FILE *f1;

 if (Geomag_Type == 1)
  if ((f1 = fopen("Geomag_Equator.DAT", "w+"))== NULL)
  {
      fprintf(stderr, "Cannot open input file.\n");
      return;
  }

  try{
         //updating geomagnetic latitude and longitude
         lat_pn_mag = (78.8 + 4.283E-2*(((MJD+fracOfDay)-46066)/365.25))*pi/180.0;
         lamb_pn_mag = (289.1 - 1.413E-2*(((MJD+fracOfDay)-46066)/365.25))*pi/180.0;

         //Compute geodetic pierce point coordinates
         Pierce_Point_Coord(Xest,Yest,Zest,Xs,Ys,Zs,lat,lamb, h);

         if(read_tec==2)
         {
             //mapping VTEC from GIM to slant TEC (STEC - )
             TEC = (TEC/cos(zl))*1.0e16;
         }

         if(Geomag_Type==0){

                //using dipolar model

                //Satellite coordinates in the geomagnetic field
                Cart_Geomag(Xs,Ys,Zs,Xsm,Ysm,Zsm);  //geomagnetic cartesians by rotations

                Cart_Curv(Xsm,Ysm,Zsm,lats_dip,lambs_dip,hs_dip);

                //---------------------------------------------------------------------------------------
                //Station coordinates in the geomagnetic field
                Cart_Geomag(Xest,Yest,Zest,Xm,Ym,Zm);

                Cart_Curv(Xm,Ym,Zm,lat_dip,lamb_dip,h_dip);

                //---------------------------------------------------------------------------------------
                //Geomagnetic pierce point ionosphere coordinates
                Curv_Cart(lat_ion,lamb_ion,hion,Xion,Yion,Zion);

                Cart_Geomag(Xion,Yion,Zion,Xionm,Yionm,Zionm);

                //Geomagnetic latitude and longitude of the pierce point
                Cart_Curv(Xionm,Yionm,Zionm,lat_ion_dip,lamb_ion_dip,hion_dip);

                //Get difference among geomagnetic coordinates (Xsm - Xm) in the Geomagnetic Local system
                Sist_Local(Xm,Ym,Zm,Xsm,Ysm,Zsm,lat_dip,lamb_dip,Em,Nm,Um);

                //Compute elevation and azimuth in the geomagnetic local system
                Ang_Local(Em,Nm,Um,Azm,Elevm,Zm);

                //Geomagnetic vector
                BTJ = BtJ(lat_ion_dip,Zm,Azm);

                BTJ = BTJ*Beq;

                Iono_Second_Order(TEC,BTJ);

                Iono_Third_Order(TEC);

                return;
         }
         else if (Geomag_Type == 1){
                 //Corrected Geomagnetic Model - PIM

                 //-----------------------Satellite coordinates-----------------
                 Cart_Curv(Xs,Ys,Zs,lats,lambs,hs_dip);

                 for(lambs=-55;lambs<-30;lambs++){
                 lats = 0;
                 lambs_temp = lambs;

                 //Geodetic to geocentric
                 double lat_geoc_sat = lats;
                 Lat_Geod_Geoc(lat_geoc_sat,lat_geoc_sat,0);

                 //Decimal degrees
                 lats_temp = float(lat_geoc_sat*180.0/pi);
                // lambs_temp = float(lambs*180.0/pi);

                 //negative for west hemisphere
                 if( lambs_temp > 180.0 )
                    lambs_temp = lambs-360.0;

                 ut = (hour)+(min/60.0)+(sec/3600.0);

                 static int cont=0; //used to read data from CGLATLON.dat

                 if(cont==0)
                   ler = 1;
                 else
                   ler = 0;

                 cont++;
                 //Corrected Geomagnetic Model - CGM

                 //for linux
//               geodcgeomag_(&MJD,&ut,&lats_temp,&lambs_temp,&lats_dip_temp,&lambs_dip_temp,&MLT,&ler);

                 //for windows
                // GEODCGEOMAG(&MJD,&ut,&lats_temp,&lambs_temp,&lats_dip_temp,&lambs_dip_temp,&MLT,&ler);


                 lats_dip = (lats_dip_temp*pi/180.0); //radians

                 //negative for west hemisphere
                 if( lambs_dip_temp > 180 )
                   lambs_dip = ((lambs_dip_temp - 360.0)*pi/180.0);
                 else
                   lambs_dip = (lambs_dip_temp*pi/180.0);

                 //Geocentric to Geodetic
                 Lat_Geod_Geoc(lats_dip,lats_dip,1);

                 //Geomagnetic cartesian
                 Curv_Cart(lats_dip,lambs_dip,hs_dip,Xsm,Ysm,Zsm);


                  fprintf(f1,"%lf %lf\n",lambs,lats_dip*180.0/pi);


                 }
                 fclose(f1);
                 //-----------------------------------------------------------------------
                 //---------------------------Station----------------------------
                 //Geodetic to Geocentric
                 double lat_geoc_est = lat;
                 Lat_Geod_Geoc(lat_geoc_est,lat_geoc_est,0);
                 //Decimal degrees
                 lat_temp = float(lat_geoc_est*180.0/pi);
                 lamb_temp = float(lamb*180.0/pi);
                 h_dip = h;

                 //negative for west hemisphere
                 if( lamb_temp > 180 )lamb_temp = lamb_temp-360.0;

                 //Corrected Geomagnetic Model - CGM
                 //for linux
//               geodcgeomag_(&MJD,&ut,&lat_temp,&lamb_temp,&lat_dip_temp,&lamb_dip_temp,&MLT,&ler);

                 //for windows
                 //GEODCGEOMAG(&MJD,&ut,&lat_temp,&lamb_temp,&lat_dip_temp,&lamb_dip_temp,&MLT,&ler);

                 lat_dip = (lat_dip_temp*pi/180.0); //radians

                 //Geocentric to Geodetic
                 Lat_Geod_Geoc(lat_dip,lat_dip,1);

                 //negative for west hemisphere
                 if( lamb_dip_temp > 180 )
                   lamb_dip = ((lamb_dip_temp - 360.0)*pi/180.0);
                 else
                   lamb_dip = (lamb_dip_temp*pi/180.0);

                 //Geomagnetic cartesian
                 Curv_Cart(lat_dip,lamb_dip,h_dip,Xm,Ym,Zm);

                 //-----------------------------------------------------------------------

                 //-----------------------Pierce point--------------------------
                 //Geodetic to Geocentric
                 double lat_geoc_piono = lat_ion;
                 Lat_Geod_Geoc(lat_geoc_piono,lat_geoc_piono,0);

                 //Decimal degrees
                 lat_temp = float(lat_geoc_piono*180.0/pi);
                 lamb_temp = float(lamb_ion*180.0/pi);

                 //negative for west hemisphere
                 if( lamb_temp > 180 )
                   lamb_temp = lamb_temp-360.0;

                 //Correct Geomagnetic Model - CGM

                 //for linux
//                geodcgeomag_(&MJD,&ut,&lat_temp,&lamb_temp,&lat_ion_dip_temp,&lamb_ion_dip_temp,&MLT,&ler);

                 //for windows
                 //GEODCGEOMAG(&MJD,&ut,&lat_temp,&lamb_temp,&lat_ion_dip_temp,&lamb_ion_dip_temp,&MLT,&ler);

                 lat_ion_dip = (lat_ion_dip_temp*pi/180.0); //Radians

                 //Geocentric to geodetic
                 Lat_Geod_Geoc(lat_ion_dip,lat_ion_dip,1);

                 //negative for west hemisphere
                 if( lamb_ion_dip_temp > 180 )
                   lamb_ion_dip = ((lamb_ion_dip_temp-360.0)*pi/180.0);
                 else
                   lamb_ion_dip = (lamb_ion_dip_temp*pi/180.0);


                 //Get difference among geomagnetic coordinates (Xsm - Xm) in the Geomagnetic Local system
                 Sist_Local(Xm,Ym,Zm,Xsm,Ysm,Zsm,lat_dip,lamb_dip,Em,Nm,Um);

                 //Compute elevation and azimuth in the geomagnetic local system
                 Ang_Local(Em,Nm,Um,Azm,Elevm,Zm);

                 //Geomagnetic vector
                 BTJ = BtJ(lat_ion_dip,Zm,Azm);

                 BTJ = BTJ*Beq;

                 Iono_Second_Order(TEC,BTJ);

                 Iono_Third_Order(TEC);

                 return;
         }
         else if(Geomag_Type == 2){

                 if(!IGRF11( year,month,day,hour,min,sec,
                             (lat_ion*180.0/pi), (lamb_ion*180.0/pi),(hion/1000.0),
                             Name_IGRF_Model) )
                 {
                     system("pause");
                     return;
                 }
                 double Nm = Get_X();
                 double Em = Get_Y();
                 double Um = Get_Z();
                 double V_Proj=0;
/*
cout << "#:-------------------------------------------------------" << endl;
cout << "#:MJD, fmjd, TEC: " << MJD << " " << fracOfDay << " " << TEC << endl;
cout << "#:Nm = " << setw(30) << setprecision(15) << Nm << endl;
cout << "#:Em = " << setw(30) << setprecision(15) << Em << endl;
cout << "#:Um = " << setw(30) << setprecision(15) << Um << endl;
cout << "#:Xs = " << setw(30) << setprecision(15) << Xs << endl;
cout << "#:Ys = " << setw(30) << setprecision(15) << Ys << endl;
cout << "#:Zs = " << setw(30) << setprecision(15) << Zs << endl;
*/
                 Project_Vector( Xs,Ys,Zs,Xest,Yest,Zest,lat, lamb, h, Nm, Em, Um,V_Proj);

//cout << "#:V_Proj = " << setw(30) << setprecision(15) << V_Proj << endl;


                 Iono_Second_Order(TEC,V_Proj);


//cout << "#:I2_L1 = " << setw(30) << setprecision(15) << I2_L1 << endl;
//cout << "#:I2_L2 = " << setw(30) << setprecision(15) << I2_L2 << endl;

                 Iono_Third_Order(TEC);

                 return;

         }//if

     }catch(...) {
          cout<<"Problem computing ionospheric effects..."<<endl;
      }

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Geod_Dipolo(double lat,double lamb,double &latdip,double &lambdip)
{/*-------------------------------------------------------------------
    Purpose: Transform Geodetic coordinates to dipolar geomagnetic coordinates
    -------------------------------------------------------------------
    Input: lat - latitude    (radians)  between -90 e 90
           lamb - longitude  (radians)  between -180 e 180
           lat0 - latitde - magnetic north pole (radians)
           lamb0 - longitude magnetic north pole (radians)

    Output: latdip - Geomagnetic latitude
            lambdip - Geomagnetic longitude
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    double senlat,senlat0,coslat,coslat0,deltalamb,sendeltalamb,cosdeltalamb,
    lamb_temp;

    lamb_temp = lamb;

    senlat = sin(lat);
    senlat0 = sin(lat_pn_mag);
    coslat = cos(lat);
    coslat0 = cos(lat_pn_mag);
    deltalamb = lamb_temp - lamb_pn_mag;
    sendeltalamb = sin(deltalamb);
    cosdeltalamb = cos(deltalamb);

    latdip = asin( senlat*senlat0 + coslat*coslat0*cosdeltalamb);
    lambdip = asin( (coslat*sendeltalamb)/(cos(latdip)) );

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Calcula_N(double lat,double &N)
{/*-------------------------------------------------------------------
    Purpose: To compute Geodetic Normal
    -------------------------------------------------------------------
    Input: Geodetic latitude

    Output: geodetic Normal
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation: semi-major axis (a) and first excentricity defined in
                 the constructor method
    -------------------------------------------------------------------*/

    N = a/sqrt(1.0- (pow(exc,2)*pow(sin(lat),2) ) );

  return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Curv_Cart(double lat,double longi,double h,double &X,
 	                 double &Y, double &Z)
{/*-------------------------------------------------------------------
    Purpose: Transform geodetic to cartesian coordinates
    -------------------------------------------------------------------
    Input: lat -  latitude
           longi - longitude
           h - elipsoidal altitude

    Output:X, Y e Z - Cartesian coordinates
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

 double N;

    //Geodetic normal
    Calcula_N(lat,N);

    X = (N+h)*cos(lat)*cos(longi);
    Y = (N+h)*cos(lat)*sin(longi);
    Z = ( N*(1.0-pow(exc,2))+h)*sin(lat);

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Curv_EsfeGeoc_Cart(double lat_geoc,double long_geoc,double R,
                                double &XEsf,double &YEsf,double &ZEsf)
{ /*-------------------------------------------------------------------
    Purpose: latitude and longitude to cartesian (espherical system)
    -------------------------------------------------------------------
    Input: lat_geoc - spherical geocentric latitude
           long_geoc - spherical geocentric longitude
           R - geocentric spherical axis

    Output: XEsf, YEsf e ZEsf - spherical cartesians
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    XEsf = R * cos(lat_geoc);
    YEsf = R * cos(lat_geoc) * tan(long_geoc);
    ZEsf = R * sin(lat_geoc);

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::CartEsfeGeoc_CartGeod(double XEsf,double YEsf,double ZEsf,double quisi,
                                 double &XGeod,double &YGeod,double &ZGeod)
{/*-------------------------------------------------------------------
    Purpose: geocentric spherical cartesian to geodetic cartesian
    -------------------------------------------------------------------
    Input: XEsf, YEsf e ZEsf - spherical cartesians
           quisi - difference between geocentric and geodetic latitude (radians)

    Output: XGeod, YGeod e ZGeod - Geodetic cartesians
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    XGeod = XEsf * cos(quisi) - ZEsf * sin(quisi);
    YGeod = YEsf;
    ZGeod = -XEsf * sin(quisi) + ZEsf * cos(quisi);

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Cart_Curv(double x,double y,double z,double &fi,double &lamb,double &h)
{/*-------------------------------------------------------------------
    Purpose: Cartesian to geodetic latitude, longitude and height
             Iterative way
    -------------------------------------------------------------------
    Input: x, y e z - Geodetic cartesians

    Output: fi, lamb e h - latitude, longitude e geometric (ellipsoidal)
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
double p,N,e2,fi1,valor;

    try{
           double soma = pow(x,2)+pow(y,2);
           if(soma <=0)
           {
               cout<<endl<<"Problems in Cart_Curv subroutine"<<endl<<
                           "please, contact: "<<"haroldoh2o@gmail.com"<<endl;

               system("pause");
               exit(-1);
           }

           p = sqrt(soma);

           //approximated latitude
           e2 = pow(exc,2);
           fi1 = atan ( (z/p) * ( 1.0/(1.0-e2) ) );
           fi = 0.0;
           int i=0;

           while(fabs(fi - fi1)> (0.00000000001)||(i>10))
           {   if (i>0)
                 fi1 = fi;

                Calcula_N(fi1,N);
                h = (p/cos(fi1)) - N;
                fi = atan2 ( (z/p) * 1.0,( 1.0 - ( (e2*N)/(N+h) ) ) );
                i++;
           }//while

           Calcula_N(fi1,N);

           h = (p/cos(fi1)) - N; //Geometric altitude
           lamb = 2.0* atan2((y),(x+p)); //Longitude
        }catch(...)
        {
            cout<<endl<<"Problems in Cart_Curv subroutine"<<endl<<
                        "please, contact: "<<"haroldoh2o@gmail.com"<<endl;

            system("pause");
            exit(-1);
        }

  return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Cart_Geomag(double Xest,double Yest, double Zest,
                       double &Xm,double &Ym,double &Zm)
{/*-------------------------------------------------------------------
    Purpose: Geodetic cartesian to geomagnetic cartesian coordinates
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation: Equations can be foun in ODJIK, 2002
    -------------------------------------------------------------------*/

    double senlamb = sin(lamb_pn_mag),
           coslamb =  cos(lamb_pn_mag),
           senlat = sin(lat_pn_mag),
           coslat =  cos(lat_pn_mag),
           senlatcoslamb =  senlat*coslamb,
           senlatsenlamb =  senlat*senlamb,
           coslatcoslamb = coslat*coslamb,
           coslatsenlamb = coslat*senlamb;

    Xm = senlatcoslamb*Xest+
         senlatsenlamb*Yest-
         coslat*Zest;

    Ym = -senlamb*Xest + coslamb*Yest;

    Zm = coslatcoslamb*Xest+
         coslatsenlamb*Yest+
         senlat*Zest;

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Sist_Local(double xest,double yest,double zest,
                      double xsat,double ysat,double zsat,
                      double fi, double lambda,
                      double &e,double &n,double &u)
{ /*-------------------------------------------------------------------
    Purpose: Cartesian coordinates to local system coordinates
    -------------------------------------------------------------------
    Input: xestacao, yestacao, zestacao - Station coordinates

    Output: e, n, u - east, north and up in the local system coordinates
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
  double senfi,cosfi,senlamb,coslamb,senfi_coslamb,senfi_senlamb,cosfi_coslamb,
         cosfi_senlamb,lat,lamb,h,result,result1,dx,dy,dz;

    dx = xsat - xest ;
    dy = ysat - yest ;
    dz = zsat - zest ;

    senfi = sin(fi);
    cosfi = cos(fi);
    senlamb = sin(lambda);
    coslamb = cos(lambda);
    senfi_coslamb = senfi*coslamb;
    senfi_senlamb = senfi*senlamb;
    cosfi_coslamb = cosfi*coslamb;
    cosfi_senlamb = cosfi*senlamb;

    n = (-senfi_coslamb*dx) - (senfi_senlamb*dy) + (cosfi*dz);
    e = (-senlamb*dx)+(coslamb*dy);
    u = (cosfi_coslamb*dx) + (cosfi_senlamb*dy) + (senfi*dz);

 return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Ang_Local(double E, double N, double U,
                     double &Am, double &Em, double &Zm)
{ /*-------------------------------------------------------------------
    Purpose: azimuth and elevation angle in the local system
    -------------------------------------------------------------------
    Input: E, N, U - local system coordinates

    Output: Am, Em - Azimuth and elevation angule
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation: Equations can be found in (STRANG; BORRE, 1997, p. 502)

      STRANG, G.; BORRE, k. Linear Algebra, Geodesy and GPS. Wellesley-Cambrigde
      Press, 1997, 624p
    -------------------------------------------------------------------*/
  double ro;

    //distance between receiver and satellite in the local system
    ro = sqrt( pow(N,2)+pow(E,2)+pow(U,2) );

    //Normalized vector
    N = N/ro;
    E = E/ro;
    U = U/ro;

    //azimuth
    Am = atan(E/N);

    if(E>0 && N<0)
       Am = pi - fabs(Am);                //2nd quadrant
    else if(E<0 && N<0)
            Am = Am + pi;                 //3rd quadrant
         else if(E<0 && N>0)
            Am = (2.0*pi) - fabs(Am);     //4th quadrant

    while(Am > (2.0*pi))
      Am -= (2.0*pi);

    while (Am < 0)
      Am += (2.0*pi);

   // Am = atan2(E,N);

    //Elevation angle
    Em = asin(U);

    //Zenith angle
    Zm = (pi/2.0) - Em;


return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Pierce_Point_Coord(double Xest,double Yest,double Zest,
                              double Xsat,double Ysat,double Zsat,
                              double lat,double lamb, double h)
{/*-------------------------------------------------------------------
    Purpose: Pierce point coordinates
    -------------------------------------------------------------------
    Input: lat - user latitude
           Z - zenith angle
           A - azimuth

    Output: lat_ion - Pierce point latitude
            lamb_ion - Pierce point longitude
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:  References:

      ODIJK D. Fast precise GPS positioning in the presence of ionospheric delays.
      2002. 242 f. PhD dissertation, Faculty of Civil Engineering and Geosciences,
      Delft University of Technology, Delft.

      GISAWY, M. L. Development of an Ionosphere Monitoring Technique Using GPS
      Measurements for High Latitude GPS Users. 2003. 161 p.
      Thesis. University of Calgary. Calgary.
      Available in: <http://www.geomatics.ucalgary.ca/links/GradTheses.html>
      Acess in: mar. 2007

      MATSUOKA, M. T.; CAMARGO, P. O. C�lculo do TEC usando dados de receptores GPS
      de dupla frequencia par a produ��o de mapas da ionosfera para a regi�o brasileira.
      revista Brasileira de Cartografia, n. 56/01 jul. 2004.
    -------------------------------------------------------------------*/
    double temp;
    double E,N,U;
    double zz_az, zz_elv, zz_dist, pplat, pplon, zz_zp;
    int iJMJ;

    try{

           Sist_Local(Xest,Yest,Zest,Xsat,Ysat,Zsat,lat,lamb,E,N,U);

           Ang_Local(E,N,U,Az,El,Zen);

           temp = (Re/(Re+hion)) * sin(Zen);

           //cout<<temp<<endl;
           zl = asin(temp);

           //latitud of ionospheric pierce point
           temp = sin(lat)*cos(Zen-zl) + cos(lat)*sin(Zen-zl)*cos(Az);
           lat_ion = asin(temp);

           //longitud of ionospheric pierce point
           temp = (sin(Zen - zl)*sin(Az))/cos(lat_ion);
           lamb_ion = lamb + asin(temp);


/*
cout << "=========================================================" << endl;
cout << "curr c,b,Re: " << c << " " << b << " " << Re << endl;
cout << "HMrad  Az: " << Az << endl;
cout << "HMrad  El: " << El << endl;
cout << "HMrad  Zen: " << Zen << endl;
cout << "HMrad  zl: " << zl << endl;
cout << "HMrad  lat_ion: " << lat_ion << endl;
cout << "HMrad  lamb_ion: " << lamb_ion << endl;
cout << "pi = " << setw(20) << setprecision(16) << pi << endl;
cout << "HMdeg  Az: " << Az*(180.0/pi) << endl;
cout << "HMdeg  El: " << El*(180.0/pi) << endl;
cout << "HMdeg  Zen: " << Zen*(180.0/pi) << endl;
cout << "HMdeg  zl: " << zl*(180.0/pi) << endl;
cout << "HMdeg  lat_ion: " << lat_ion*(180.0/pi) << endl;
cout << "HMdeg  lamb_ion: " << lamb_ion*(180.0/pi) << endl;

*/

 azelIPP(Xest,Yest,Zest,Xsat,Ysat,Zsat,Re,hion,
 &zz_az, &zz_elv, &zz_dist, &pplat, &pplon, &zz_zp);

/*
 cout << "SHrad  pplat: " << pplat << endl;
 cout << "SHrad  pplon: " << pplon << endl;

 cout << "SHdeg  pplat: " << pplat*(180.0/pi) << endl;
 cout << "SHdeg  pplon: " << pplon*(180.0/pi) << endl;

getchar();
*/


 lat_ion = pplat;
 lamb_ion = pplon;

 zl = zz_zp;

/*
cout << "4newHMrad  lat_ion: " << lat_ion << endl;
cout << "4newHMrad  lamb_ion: " << lamb_ion << endl;
cout << "4curr c,b,Re: " << c << " " << b << " " << Re << endl;

 cout << "Type any int to continue." << endl;
 cin >> iJMJ;
 */


           //****************************************************************************
           //It also can be computed in this way:

          /*  temp = (Re/(Re+hion)) * cos(El);
            double quisi = (pi/2.0)-El-asin( temp );
            temp = sin(lat)*cos(quisi) + cos(lat)*sin(quisi)*cos(Az);
            lat_ion = asin(temp);
            temp = (sin(quisi)*sin(Az))/cos(lat_ion);
            lamb_ion = lamb + asin(temp);
          */
           //****************************************************************************

   }catch(...)
   {
     cout<<"Problemas no calculo da inosfera"<<endl;
     system("pause");
   }
}
//------------------------------------------------------------------------------------
double Class_Iono::BtJ(double lat_ion_m,double z_m,double a_m)
{ /*-------------------------------------------------------------------
    Purpose: Inner product between Bt (geomagnetic) and J (line-of-site)
             vectors
    -------------------------------------------------------------------
    Input: lat_ion_m - Geomagneitc latitude of the pierce point
           z_m - Geomagnetic zenith angle
           a_m - Geomagnetic azimuth angle

    Output: return inner product BtJ
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
   double BtJ,temp;

    temp = pow( (Re/(Re+hion)),3 );
    BtJ = ( cos(lat_ion_m)*sin(z_m)*cos(a_m)- 2.0*sin(lat_ion_m)*cos(z_m) ) * temp;

    BtJ = ( 2.0*sin(lat_ion_m)*cos(z_m) - cos(lat_ion_m)*sin(z_m)*cos(a_m)) * temp;

   return(BtJ);
}
//------------------------------------------------------------------------------------
void Class_Iono::Iono_First_Order(double TEC)
{/*-------------------------------------------------------------------
    Purpose: First order ionospheric effect
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
 double const1 = (40.3*pow(lamb_l1,2))/pow(c,2),
         const2 = (40.3*pow(lamb_l2,2))/pow(c,2);


   I1_L1 = const1*TEC;
   I1_L2 = const2*TEC;
return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Iono_Second_Order(double TEC,double BtJ)
{/*-------------------------------------------------------------------
    Purpose: Second order ionospheric effect
    -------------------------------------------------------------------
    Input: TEC - Slant TEC in TECU
           BtJ - Inner product between Bt (geomagnetic) and J (line-of-site) vectors

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

  double  e = 1.60218E-19,  //Electron charge - Coloumb
          A = 80.6,         // Constant  m^3/s^2
          me = 9.10939E-31,  //electron Mass - kg
          const1 = ( ( e*A*pow(lamb_l1,3) ) / ( pow(c,3)*2.0*pi*me ) ),
          const2 = ( ( e*A*pow(lamb_l2,3) ) / ( pow(c,3)*2.0*pi*me ) );

    I2_L1 = const1*BtJ*TEC;
    I2_L2 = const2*BtJ*TEC;

return;
}
//------------------------------------------------------------------------------------
void Class_Iono::Iono_Third_Order(double TEC)
{/*-------------------------------------------------------------------
    Purpose: Third order ionospheric effect
    -------------------------------------------------------------------
    Input: TEC - Slant TEC - TECU

    Output:
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: July of 2007
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
 double  A = 80.6,         // Constant
         conste = ( ( 3.0*pow(A,2)*0.66 )/( 8.0*pow(c,4) ) );

    //electron density maximum in function of the TEC
    Ne_max = ( ( (20.0-6.0)*1.0E12 )/
               ( (4.55-1.38)*1.0E18) )*TEC;

    I3_L1 = conste*pow(lamb_l1,4 )*Ne_max*TEC;
    I3_L2 = conste*pow(lamb_l2,4 )*Ne_max*TEC;

return;
}
//------------------------------------------------------------------------------------


}//namespace CIONO{






