/*-------------------------------------------------------------------
    Purpose:  C++ Class to compute geomagnetic components (X, Y, Z) from IGRF 11 model
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: August of 2010
    -------------------------------------------------------------------
    Observation: This class was developed to support the source code available at:
                 http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    -------------------------------------------------------------------*/
    
#ifndef CLASSIGRF11_H
#define CLASSIGRF11_H

//SH #include <stdio.h>
//SH #include <conio.h>
//SH #include <stdlib.h>
//SH #include <string.h>
//SH #include <ctype.h>
//SH #include <math.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>

//int my_isnan(double d)
//{
//  return (d != d);              /* IEEE: only NaN is not equal to itself */
//}

#define NaN log(-1.0)
#define FT2KM (1.0/0.0003048)
#define PI 3.141592654
#define RAD2DEG (180.0/PI)


#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define IEXT 0
#define FALSE 0
#define TRUE 1                  /* constants */
#define RECL 81

#define MAXINBUFF RECL+14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF-2
/** Max to read 2 less than total size (just to be safe) **/

#define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define EXT_COEFF1 (double)0
#define EXT_COEFF2 (double)0
#define EXT_COEFF3 (double)0

#define MAXDEG 13
#define MAXCOEFF (MAXDEG*(MAXDEG+2)+1) /* index starts with 1!, (from old Fortran?) */

namespace IGRF11{
class Class_IGRF11{

   public:

           Class_IGRF11();
           ~Class_IGRF11();

           int Read_IGRF_Model(char  Name_IGRF_Model[], double &minyr, double &maxyr,
                               int   &modelI,      /* Which model (Index) */
                               int   &nmodel,      /* Number of models in file */
                               double yrmax[MAXMOD], double altmin[MAXMOD], double altmax[MAXMOD],
                               int   max1[MAXMOD], int   max2[MAXMOD], int   max3[MAXMOD],
                               long  irec_pos[MAXMOD], char  model[MAXMOD][9],
                               double epoch[MAXMOD], double yrmin[MAXMOD]);

           int IGRF11(int year, int month, int day, int hour, int min, double sec,
                      double latitude, double longitude, double h,char Name_IGRF_Model[]);

           double gh1[MAXCOEFF];
           double gh2[MAXCOEFF];
           double gha[MAXCOEFF];              /* Geomag global variables */
           double ghb[MAXCOEFF];

           double dtemp,ftemp,htemp,itemp;

           double xtemp,ytemp,ztemp;

           FILE *stream;  /* Pointer to specified model data file */

           //Set Method
           void Set_X(double value);
           void Set_Y(double value);
           void Set_Z(double value);
           void Set_d(double value);
           void Set_f(double value);
           void Set_h(double value);
           void Set_i(double value);

           //Get Method
           double Get_X();
           double Get_Y();
           double Get_Z();
           double Get_d();
           double Get_f();
           double Get_h();
           double Get_i();

   private:

           double degrees_to_decimal();

           double julday(int month,int day, int year);

           int interpsh(double date, double dte1, int   nmax1, double dte2,
                        int   nmax2, int   gh);

           int extrapsh(double date,double dte1,int   nmax1,int   nmax2,int   gh);

           int shval3( int igdgc,double flat,double flon,double elev,int   nmax,
                       int gh,int iext,double ext1,double ext2,double ext3);

           int dihf (int gh);

           int getshc(char file[PATH],int iflag,long int  strec,int nmax_of_gh,int gh);

           double x,y,z;
           double d,f,h,i;
           int need_to_read_model;
           int need_to_interp;


            /* Control variables */
           int   again;
           int   decyears;
           int   units;
           int   decdeg;
           int   range;
           int   counter;
           int   warn_H, warn_H_strong, warn_P;

           int   fileline;

           char  inbuff[MAXINBUFF];
           long  irec_pos[MAXMOD];
           char  model[MAXMOD][9];
           double epoch[MAXMOD];

           double yrmin[MAXMOD];
           int   nmax;
           int   igdgc;
           int   isyear;
           int   ismonth;
           int   isday;
           int   ieyear;
           int   iemonth;
           int   ieday;
           int   ilat_deg;
           int   ilat_min;
           int   ilat_sec;
           int   ilon_deg;
           int   ilon_min;
           int   ilon_sec;

           int  coords_from_file;
           int arg_err;


           //char  mdfile[PATH];
           char *begin;
           char *rest;
           char args[7][MAXREAD];
           int iarg;

           char coord_fname[PATH];
           char out_fname[PATH];
           FILE *coordfile,*outfile;
           int iline;
           int read_flag;

           double minyr;
           double maxyr;
           int   modelI;             /* Which model (Index) */
           int   nmodel;             /* Number of models in file */
           double yrmax[MAXMOD];
           double altmin[MAXMOD];
           double altmax[MAXMOD];
           int   max1[MAXMOD];
           int   max2[MAXMOD];
           int   max3[MAXMOD];

           double minalt;
           double maxalt;
           double alt;
           double sdate;
           double step;
           double syr;
           double edate;
           double ddot;
           double fdot;
           double hdot;
           double idot;
           double xdot;
           double ydot;
           double zdot;
           double warn_H_val, warn_H_strong_val;

};

}//namespace IGRF11

#endif
