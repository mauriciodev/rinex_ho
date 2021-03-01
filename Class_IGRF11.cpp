/*-------------------------------------------------------------------
    Purpose:  C++ Class to compute geomagnetic components (X, Y, Z) from IGRF 11 model
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Date: August of 2010
    -------------------------------------------------------------------
    Observation: This class was developed to suport the source code available at:
                 http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
    -------------------------------------------------------------------*/

#include "Class_IGRF11.h"

namespace IGRF11{
//------------------------------------------------------------------------------------
Class_IGRF11::Class_IGRF11()
{/*-------------------------------------------------------------------
     Purpose: Class constructor
     -------------------------------------------------------------------
     Input:

     Output:
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP

     Date: January of 2010
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

    stream = NULL;

    x=0;
    y=0;
    z=0;
    d=0;
    f=0;
    h=0;
    i=0;

    need_to_read_model = 1;
    need_to_interp =1;

    /* Control variables */
    again = 1;
    decyears = 3;
    units = 4;
    decdeg = 3;
    range = -1;
    counter = 0;
    igdgc=3;
    isyear=-1;
    ismonth=-1;
    isday=-1;
    ieyear=-1;
    iemonth=-1;
    ieday=-1;
    ilat_deg=200;
    ilat_min=200;
    ilat_sec=200;
    ilon_deg=200;
    ilon_min=200;
    ilon_sec=200;

    coords_from_file = 0;
    arg_err = 0;
    iline=0;
    alt=-999999;
    sdate=-1;
    step=-1;
    edate=-1;


}
//------------------------------------------------------------------------------------
Class_IGRF11::~Class_IGRF11()//destructor
{ /*-------------------------------------------------------------------
    Purpose: Class Destructor
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
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_X(double value)
{
    x = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_Y(double value)
{
    y = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_Z(double value)
{
    z = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_d(double value)
{
    d = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_f(double value)
{
    f = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_h(double value)
{
    h = value;
}
//------------------------------------------------------------------------------------
void Class_IGRF11::Set_i(double value)
{
    i = value;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_X()
{
    return x;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_Y()
{
    return y;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_Z()
{
    return z;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_d()
{
    return d;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_f()
{
    return f;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_h()
{
    return h;
}
//------------------------------------------------------------------------------------
double Class_IGRF11::Get_i()
{
    return i;
}
//------------------------------------------------------------------------------------
int Class_IGRF11::IGRF11(int year, int month, int day, int hour, int min, double sec,
                         double latitude, double longitude, double h,char Name_IGRF_Model[])
{
      //Original source code available at: http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
      //Adapted by Marques, H. A., 2010

      /*  Variable declaration  */
      warn_H = 0;
      warn_H_val = 99999.0;
      warn_H_strong = 0;
      warn_H_strong_val = 99999.0;
      warn_P = 0;

      //strcpy(mdfile,".\\IGRF11\\IGRF11.COF");

      //Obtain the desired model file and read the data  */
      if(need_to_read_model)
        if(Read_IGRF_Model(Name_IGRF_Model,minyr,maxyr,modelI,nmodel,yrmax,altmin,
                           altmax,max1,max2,max3,irec_pos, model,epoch,yrmin) )
        {
           need_to_read_model = 0;
        }
        else
           return 0;

        // if date specified in command line then warn if past end of validity
        if ((sdate>maxyr)&&(sdate<maxyr+1))
        {
            printf("\nWarning: The date %4.2f is out of range,\n", sdate);
            printf("         but still within one year of model expiration date.\n");
            printf("         An updated model file is available before 1.1.%4.0f\n",maxyr);
        }

        //----------------------------------------------------------------
        // printf("\nHow would you like to enter the date?\n");
        // printf("       1) In decimal years.\n");
        // printf("       2) In year, month, and day.\n");

        decyears = 1;
        //----------------------------------------------------------------

        // printf("\nWould you like output for a single date or for a range of dates?\n");
        //printf("       1) A single date.\n");
        //printf("       2) A range of dates.\n");

        range = 1;
        //----------------------------------------------------------------

        //if (decyears==1)
        //{
        //        printf("\nEnter the decimal date (%4.2f to %4.0f): ",minyr, maxyr);

           sdate = julday(month,day,year);
           sdate += double(hour + min/60.0+sec/3600.0)/24.0;
        //}

        if ((sdate<minyr)||(sdate>=maxyr+1))
        {
             ismonth=isday=isyear=0;
             printf("\nError: The date %4.2f is out of range.\n", sdate);
        }

        if ((sdate>maxyr)&&(sdate<maxyr+1))
        {
             printf("\nWarning: The date %4.2f is out of range,\n", sdate);
             printf("         but still within one year of model expiration date.\n");
             printf("         An updated model file is available before 1.1.%4.0f\n",maxyr);
        }

        // Pick model */
        for (modelI=0; modelI<nmodel; modelI++)
           if (sdate<yrmax[modelI]) break;

        if (modelI == nmodel) modelI--;           // if beyond end of last model use last model */

        // Get altitude min and max for selected model. */
        minalt=altmin[modelI];
        maxalt=altmax[modelI];

        //printf("\n\nEnter Coordinate Preferences:");
        //printf("\n    1) Geodetic (WGS84 latitude and altitude above mean sea level)");
        //printf("\n    2) Geocentric (spherical, altitude relative to Earth's center)\n");

        igdgc = 2;

        // If needed modify ranges to reflect coords. */
        if (igdgc==2)
        {
            minalt+=6371.2;  // Add radius to ranges. */
            maxalt+=6371.2;
        }

        // Get unit prefs */
        if (igdgc==1)
        {
            //while ((units>3)||(units<1))
            //{
               // printf("\n\nEnter Unit Preferences:");
               // printf("\n       1) Kilometers");
               // printf("\n       2) Meters");
               // printf("\n       3) Feet\n");
               // printf("\n                            ==> ");
               units = 2;
            //}
        }
        else units = 1; // geocentric always in km */

        // Do unit conversions if neccessary */
        if (units==2)
        {
            minalt*=1000.0;
            maxalt*=1000.0;
        }
        else if (units==3)
        {
           minalt*=FT2KM;
           maxalt*=FT2KM;
        }

        //if (igdgc==2) printf("\n\nEnter geocentric altitude in km (%.2f to %.2f): ", minalt, maxalt);
        //if (igdgc==1 && units==1) printf("\n\nEnter geodetic altitude above mean sea level in km (%.2f to %.2f): ", minalt, maxalt);
        //if (igdgc==1 && units==2) printf("\n\nEnter geodetic altitude above mean sea level in meters (%.2f to %.2f): ", minalt, maxalt);
        //if (igdgc==1 && units==3) printf("\n\nEnter geodetic altitude above mean sea level in feet (%.2f to %.2f): ", minalt, maxalt);


        //alt = 6378.137+450.000;
        //alt = 6378.137+(51.2408/1000);

        alt = 6378.137+h;

        /* Convert altitude to km */
        if (units==2)
        {
           alt *= 0.001;
        }
        else if (units==3)
        {
           alt /= FT2KM;
        }

        //printf("\n\nHow would you like to enter the latitude and longitude?:");
        //printf("\n       1) In decimal degrees.");
        //printf("\n       2) In degrees, minutes, and seconds.\n");

        decdeg = 1;

        if (decdeg==1)
        {
              //printf("\n\nEnter the decimal latitude (-90 to 90) (- for Southern hemisphere).\n");
              //printf("\n\nEnter the decimal longitude (-180 to 180) (- for Western hemisphere).\n");
        } /* if (decdeg==1) */


        if(need_to_interp)
        /** This will compute everything needed for 1 point in time. **/
        if (max2[modelI] == 0)
        {
           getshc(Name_IGRF_Model, 1, irec_pos[modelI], max1[modelI], 1);
           getshc(Name_IGRF_Model, 1, irec_pos[modelI+1], max1[modelI+1], 2);
           nmax = interpsh(sdate, yrmin[modelI], max1[modelI],
                           yrmin[modelI+1], max1[modelI+1], 3);
           nmax = interpsh(sdate+1, yrmin[modelI] , max1[modelI],
                           yrmin[modelI+1], max1[modelI+1],4);

           need_to_interp=0;
        }
        else
        {
           getshc(Name_IGRF_Model, 1, irec_pos[modelI], max1[modelI], 1);
           getshc(Name_IGRF_Model, 0, irec_pos[modelI], max2[modelI], 2);
           nmax = extrapsh(sdate, epoch[modelI], max1[modelI], max2[modelI], 3);
           nmax = extrapsh(sdate+1, epoch[modelI], max1[modelI], max2[modelI], 4);

           need_to_interp=0;
        }


        /* Do the first calculations */
        shval3(igdgc, latitude, longitude, alt, nmax, 3,
               IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
        dihf(3);
        shval3(igdgc, latitude, longitude, alt, nmax, 4,
               IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3);
        dihf(4);
         ddot = ((dtemp - d)*RAD2DEG);
      if (ddot > 180.0) ddot -= 360.0;
      if (ddot <= -180.0) ddot += 360.0;
      ddot *= 60.0;
      
      idot = ((itemp - i)*RAD2DEG)*60;
      d = d*(RAD2DEG);   i = i*(RAD2DEG);
      hdot = htemp - h;   xdot = xtemp - x;
      ydot = ytemp - y;   zdot = ztemp - z;
      fdot = ftemp - f;
      
      /* deal with geographic and magnetic poles */
      if (h < 100.0) /* at magnetic poles */
        {
          d = NaN;
          ddot = NaN;
          /* while rest is ok */
        }

      if (h < 1000.0)
        {
          warn_H = 0;
          warn_H_strong = 1;
          if (h<warn_H_strong_val) warn_H_strong_val = h;
        }
      else if (h < 5000.0 && !warn_H_strong) 
        {
          warn_H = 1;
          if (h<warn_H_val) warn_H_val = h;
        }
      
      if (90.0-fabs(latitude) <= 0.001) /* at geographic poles */
        {
          x = NaN;
          y = NaN;
          d = NaN;
          xdot = NaN;
          ydot = NaN;
          ddot = NaN;
          warn_P = 1;
          warn_H = 0;
          warn_H_strong = 0;
          /* while rest is ok */
        }
      
      /** Above will compute everything for 1 point in time.  **/


       //printf("\n\n\n  Model: %s \n", model[modelI]);
       //   if (decdeg==1)
       //     {
       //       printf("  Latitude: %4.2f deg\n", latitude);
       //       printf("  Longitude: %4.2f deg\n", longitude);
       //     }
       //   else
       //     {
       //       printf("  Latitude: %d deg, %d min, %d sec\n",
       //              ilat_deg,ilat_min, ilat_sec);
       //       printf("  Longitude: %d deg,  %d min, %d sec\n",
       //              ilon_deg, ilon_min, ilon_sec);
       //     }
       //   printf("  Altitude: ");
       //   if (units==1)
       //     printf("%.2f km\n", alt);
       //   else if (units==2)
       //     printf("%.2f meters\n", alt*1000.0);
       //   else
       //     printf("%.2f ft\n", (alt*FT2KM));

          //if (range==1)
          //{
            //  printf("  Date of Interest: ");
            //  if (decyears==1)
           //     printf(" %4.2f\n\n", sdate);
           //   else
           //     printf("%d-%d-%d (yyyy-mm-dd)\n\n",  isyear, ismonth, isday);

           //  printf("  Date of Interest: ");
           //   if (decyears==1)
           //     printf(" %4.2f\n\n", sdate);
           //   else
           //     printf("%d-%d-%d (yyyy-mm-dd)\n\n",  isyear, ismonth, isday);
              
              //print_header();
              //print_result(sdate,d, i, h, x, y, z, f);
              //print_long_dashed_line();
              //print_header_sv();
              //print_result_sv(sdate,ddot,idot,hdot,xdot,ydot,zdot,fdot);
              //print_dashed_line();
              
          //} /* if range == 1 */

    //system("pause");
    return 1;
}
//------------------------------------------------------------------------------------
int Class_IGRF11::Read_IGRF_Model(char  Name_IGRF_Model[],
                    double &minyr,
                    double &maxyr,
                    int   &modelI,             /* Which model (Index) */
                    int   &nmodel,             /* Number of models in file */
                    double yrmax[MAXMOD],
                    double altmin[MAXMOD],
                    double altmax[MAXMOD],
                    int   max1[MAXMOD],
                    int   max2[MAXMOD],
                    int   max3[MAXMOD],
                    long  irec_pos[MAXMOD],
                    char  model[MAXMOD][9],
                    double epoch[MAXMOD],
                    double yrmin[MAXMOD])
{
     //Original source code available at: http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
      //Adapted by Marques, H. A., 2010

char  inbuff[MAXINBUFF];
int   fileline;



  /* Initializations. */

  inbuff[MAXREAD+1]='\0';  /* Just to protect mem. */
  inbuff[MAXINBUFF-1]='\0';  /* Just to protect mem. */




    if (!(stream = fopen(Name_IGRF_Model, "rt")))
    {  printf("\nError opening file %s.", Name_IGRF_Model);
       exit;
    }

    rewind(stream);

    fileline = 0;                            /* First line will be 1 */
    modelI = -1;

          while (fgets(inbuff,MAXREAD,stream))     /* While not end of file
                                                   * read to end of line or buffer */
          {
              fileline++;                           /* On new line */


              if (strlen(inbuff) != RECL)       /* IF incorrect record size */
              {
                  printf("Corrupt record in file %s on line %d.\n", Name_IGRF_Model, fileline);
                  fclose(stream);
                  exit(5);
              }

              /* old statement Dec 1999 */
              /*       if (!strncmp(inbuff,"    ",4)){         /* If 1st 4 chars are spaces */
              /* New statement Dec 1999 changed by wmd  required by year 2000 models */
              if (!strncmp(inbuff,"   ",3))         /* If 1st 3 chars are spaces */
              {
                  modelI++;                           /* New model */

                  if (modelI > MAXMOD)                /* If too many headers */
                  {
                      printf("Too many models in file %s on line %d.", Name_IGRF_Model, fileline);
                      fclose(stream);
                      exit(6);
                  }

                  irec_pos[modelI]=ftell(stream);
                  /* Get fields from buffer into individual vars.  */
                  sscanf(inbuff, "%s%lg%d%d%d%lg%lg%lg%lg", model[modelI], &epoch[modelI],
                         &max1[modelI], &max2[modelI], &max3[modelI], &yrmin[modelI],
                         &yrmax[modelI], &altmin[modelI], &altmax[modelI]);

                  /* Compute date range for all models */
                  if (modelI == 0)                    /*If first model */
                  {
                      minyr=yrmin[0];
                      maxyr=yrmax[0];
                  }
                  else
                  {
                      if (yrmin[modelI]<minyr)
                      {
                         minyr=yrmin[modelI];
                      }
                      if (yrmax[modelI]>maxyr){
                        maxyr=yrmax[modelI];
                      }
                  } /* if modelI != 0 */

                } /* If 1st 3 chars are spaces */

            } /* While not end of model file */

            nmodel = modelI + 1;
            fclose(stream);

return 1;
}
//------------------------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/*                           Subroutine getshc                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Reads spherical harmonic coefficients from the specified             */
/*     model into an array.                                                 */
/*                                                                          */
/*     Input:                                                               */
/*           stream     - Logical unit number                               */
/*           iflag      - Flag for SV equal to ) or not equal to 0          */
/*                        for designated read statements                    */
/*           strec      - Starting record number to read from model         */
/*           nmax_of_gh - Maximum degree and order of model                 */
/*                                                                          */
/*     Output:                                                              */
/*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
/*                        harmonic coefficients                             */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Bill Flanagan                                                  */
/*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 15, 1988                                                */
/*                                                                          */
/****************************************************************************/

int Class_IGRF11::getshc(char file[PATH],int iflag,long int  strec,int nmax_of_gh,int gh)
{
  char  inbuff[MAXINBUFF];
  char irat[9];
  int ii,m,n,mm,nn;
  int ios;
  int line_num;
  double g,hh;
  double trash;
  
  stream = fopen(file, "rt");
  if (stream == NULL)
    {
      printf("\nError on opening file %s", file);
    }
  else
    {
      ii = 0;
      ios = 0;
      fseek(stream,strec,SEEK_SET);
      for ( nn = 1; nn <= nmax_of_gh; ++nn)
        {
          for (mm = 0; mm <= nn; ++mm)
            {
              if (iflag == 1)
                {
                  fgets(inbuff, MAXREAD, stream);
                  sscanf(inbuff, "%d%d%lg%lg%lg%lg%s%d",
                         &n, &m, &g, &hh, &trash, &trash, irat, &line_num);
                }
              else
                {
                  fgets(inbuff, MAXREAD, stream);
                  sscanf(inbuff, "%d%d%lg%lg%lg%lg%s%d",
                         &n, &m, &trash, &trash, &g, &hh, irat, &line_num);
                }
              if ((nn != n) || (mm != m))
                {
                  ios = -2;
                  fclose(stream);
                  return(ios);
                }
              ii = ii + 1;
              switch(gh)
                {
                case 1:  gh1[ii] = g;
                  break;
                case 2:  gh2[ii] = g;
                  break;
                default: printf("\nError in subroutine getshc");
                  break;
                }
              if (m != 0)
                {
                  ii = ii+ 1;
                  switch(gh)
                    {
                    case 1:  gh1[ii] = hh;
                      break;
                    case 2:  gh2[ii] = hh;
                      break;
                    default: printf("\nError in subroutine getshc");
                      break;
                    }
                }
            }
        }
    }
  fclose(stream);
  return(ios);
}
//------------------------------------------------------------------------------------

/****************************************************************************/
/*                                                                          */
/*                           Subroutine extrapsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Extrapolates linearly a spherical harmonic model with a              */
/*     rate-of-change model.                                                */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of base model                                  */
/*           nmax1    - maximum degree and order of base model              */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of base model                 */
/*           nmax2    - maximum degree and order of rate-of-change model    */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of rate-of-change model       */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - Schmidt quasi-normal internal spherical             */
/*                    harmonic coefficients                                 */
/*           nmax   - maximum degree and order of resulting model           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 16, 1988                                                */
/*                                                                          */
/****************************************************************************/

int Class_IGRF11::extrapsh(double date, double dte1, int nmax1, int nmax2, int gh)
{
  int   nmax;
  int   k, l;
  int   ii;
  double factor;
  
  factor = date - dte1;
  if (nmax1 == nmax2)
    {
      k =  nmax1 * (nmax1 + 2);
      nmax = nmax1;
    }
  else
    {
      if (nmax1 > nmax2)
        {
          k = nmax2 * (nmax2 + 2);
          l = nmax1 * (nmax1 + 2);
          switch(gh)
            {
            case 3:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  gha[ii] = gh1[ii];
                }
              break;
            case 4:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  ghb[ii] = gh1[ii];
                }
              break;
            default: printf("\nError in subroutine extrapsh");
              break;
            }
          nmax = nmax1;
        }
      else
        {
          k = nmax1 * (nmax1 + 2);
          l = nmax2 * (nmax2 + 2);
          switch(gh)
            {
            case 3:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  gha[ii] = factor * gh2[ii];
                }
              break;
            case 4:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  ghb[ii] = factor * gh2[ii];
                }
              break;
            default: printf("\nError in subroutine extrapsh");
              break;
            }
          nmax = nmax2;
        }
    }
  switch(gh)
    {
    case 3:  for ( ii = 1; ii <= k; ++ii)
        {
          gha[ii] = gh1[ii] + factor * gh2[ii];
        }
      break;
    case 4:  for ( ii = 1; ii <= k; ++ii)
        {
          ghb[ii] = gh1[ii] + factor * gh2[ii];
        }
      break;
    default: printf("\nError in subroutine extrapsh");
      break;
    }
  return(nmax);
}
//------------------------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/*                           Subroutine interpsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Interpolates linearly, in time, between two spherical harmonic       */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of earlier model                               */
/*           nmax1    - maximum degree and order of earlier model           */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of earlier model              */
/*           dte2     - date of later model                                 */
/*           nmax2    - maximum degree and order of later model             */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of internal model             */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - coefficients of resulting model                     */
/*           nmax     - maximum degree and order of resulting model         */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/


int Class_IGRF11::interpsh(double date, double dte1,int nmax1, double dte2,int nmax2,int gh)
{
  int   nmax;
  int   k, l;
  int   ii;
  double factor;
  
  factor = (date - dte1) / (dte2 - dte1);
  if (nmax1 == nmax2)
    {
      k =  nmax1 * (nmax1 + 2);
      nmax = nmax1;
    }
  else
    {
      if (nmax1 > nmax2)
        {
          k = nmax2 * (nmax2 + 2);
          l = nmax1 * (nmax1 + 2);
          switch(gh)
            {
            case 3:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  gha[ii] = gh1[ii] + factor * (-gh1[ii]);
                }
              break;
            case 4:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  ghb[ii] = gh1[ii] + factor * (-gh1[ii]);
                }
              break;
            default: printf("\nError in subroutine extrapsh");
              break;
            }
          nmax = nmax1;
        }
      else
        {
          k = nmax1 * (nmax1 + 2);
          l = nmax2 * (nmax2 + 2);
          switch(gh)
            {
            case 3:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  gha[ii] = factor * gh2[ii];
                }
              break;
            case 4:  for ( ii = k + 1; ii <= l; ++ii)
                {
                  ghb[ii] = factor * gh2[ii];
                }
              break;
            default: printf("\nError in subroutine extrapsh");
              break;
            }
          nmax = nmax2;
        }
    }
  switch(gh)
    {
    case 3:  for ( ii = 1; ii <= k; ++ii)
        {
          gha[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
        }
      break;
    case 4:  for ( ii = 1; ii <= k; ++ii)
        {
          ghb[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii]);
        }
      break;
    default: printf("\nError in subroutine extrapsh");
      break;
    }
  return(nmax);
}


//------------------------------------------------------------------------------------


/****************************************************************************/
/*                                                                          */
/*                           Subroutine shval3                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Calculates field components from spherical harmonic (sh)             */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           igdgc     - indicates coordinate system used; set equal        */
/*                       to 1 if geodetic, 2 if geocentric                  */
/*           latitude  - north latitude, in degrees                         */
/*           longitude - east longitude, in degrees                         */
/*           elev      - WGS84 altitude above ellipsoid (igdgc=1), or       */
/*                       radial distance from earth's center (igdgc=2)      */
/*           a2,b2     - squares of semi-major and semi-minor axes of       */
/*                       the reference spheroid used for transforming       */
/*                       between geodetic and geocentric coordinates        */
/*                       or components                                      */
/*           nmax      - maximum degree and order of coefficients           */
/*           iext      - external coefficients flag (=0 if none)            */
/*           ext1,2,3  - the three 1st-degree external coefficients         */
/*                       (not used if iext = 0)                             */
/*                                                                          */
/*     Output:                                                              */
/*           x         - northward component                                */
/*           y         - eastward component                                 */
/*           z         - vertically-downward component                      */
/*                                                                          */
/*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  */
/*     report no. 71/1, institute of geological sciences, U.K.              */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Norman W. Peddie                                               */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/

int Class_IGRF11::shval3( int igdgc,double flat,double flon,double elev,int   nmax,
                          int gh,int iext,double ext1,double ext2,double ext3)
{
  double earths_radius = 6371.2;
  double dtr = 0.01745329;
  double slat;
  double clat;
  double ratio;
  double aa, bb, cc, dd;
  double sd;
  double cd;
  double r;
  double a2;
  double b2;
  double rr;
  double fm,fn;
  double sl[14];
  double cl[14];
  double p[119];
  double q[119];
  int ii,j,k,l,m,n;
  int npq;
  int ios;
  double argument;
  double power;
  a2 = 40680631.59;            /* WGS84 */
  b2 = 40408299.98;            /* WGS84 */
  ios = 0;
  r = elev;
  argument = flat * dtr;
  slat = sin( argument );
  if ((90.0 - flat) < 0.001)
    {
      aa = 89.999;            /*  300 ft. from North pole  */
    }
  else
    {
      if ((90.0 + flat) < 0.001)
        {
          aa = -89.999;        /*  300 ft. from South pole  */
        }
      else
        {
          aa = flat;
        }
    }
  argument = aa * dtr;
  clat = cos( argument );
  argument = flon * dtr;
  sl[1] = sin( argument );
  cl[1] = cos( argument );
  switch(gh)
    {
    case 3:  x = 0;
      y = 0;
      z = 0;
      break;
    case 4:  xtemp = 0;
      ytemp = 0;
      ztemp = 0;
      break;
    default: printf("\nError in subroutine shval3");
      break;
    }
  sd = 0.0;
  cd = 1.0;
  l = 1;
  n = 0;
  m = 1;
  npq = (nmax * (nmax + 3)) / 2;
  if (igdgc == 1)
    {
      aa = a2 * clat * clat;
      bb = b2 * slat * slat;
      cc = aa + bb;
      argument = cc;
      dd = sqrt( argument );
      argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc;
      r = sqrt( argument );
      cd = (elev + dd) / r;
      sd = (a2 - b2) / dd * slat * clat / r;
      aa = slat;
      slat = slat * cd - clat * sd;
      clat = clat * cd + aa * sd;
    }
  ratio = earths_radius / r;
  argument = 3.0;
  aa = sqrt( argument );
  p[1] = 2.0 * slat;
  p[2] = 2.0 * clat;
  p[3] = 4.5 * slat * slat - 1.5;
  p[4] = 3.0 * aa * clat * slat;
  q[1] = -clat;
  q[2] = slat;
  q[3] = -3.0 * clat * slat;
  q[4] = aa * (slat * slat - clat * clat);
  for ( k = 1; k <= npq; ++k)
    {
      if (n < m)
        {
          m = 0;
          n = n + 1;
          argument = ratio;
          power =  n + 2;
          rr = pow(argument,power);
          fn = n;
        }
      fm = m;
      if (k >= 5)
        {
          if (m == n)
            {
              argument = (1.0 - 0.5/fm);
              aa = sqrt( argument );
              j = k - n - 1;
              p[k] = (1.0 + 1.0/fm) * aa * clat * p[j];
              q[k] = aa * (clat * q[j] + slat/fm * p[j]);
              sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
              cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
            }
          else
            {
              argument = fn*fn - fm*fm;
              aa = sqrt( argument );
              argument = ((fn - 1.0)*(fn-1.0)) - (fm * fm);
              bb = sqrt( argument )/aa;
              cc = (2.0 * fn - 1.0)/aa;
              ii = k - n;
              j = k - 2 * n + 1;
              p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
              q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
            }
        }
      switch(gh)
        {
        case 3:  aa = rr * gha[l];
          break;
        case 4:  aa = rr * ghb[l];
          break;
        default: printf("\nError in subroutine shval3");
          break;
        }
      if (m == 0)
        {
          switch(gh)
            {
            case 3:  x = x + aa * q[k];
              z = z - aa * p[k];
              break;
            case 4:  xtemp = xtemp + aa * q[k];
              ztemp = ztemp - aa * p[k];
              break;
            default: printf("\nError in subroutine shval3");
              break;
            }
          l = l + 1;
        }
      else
        {
          switch(gh)
            {
            case 3:  bb = rr * gha[l+1];
              cc = aa * cl[m] + bb * sl[m];
              x = x + cc * q[k];
              z = z - cc * p[k];
              if (clat > 0)
                {
                  y = y + (aa * sl[m] - bb * cl[m]) *
                    fm * p[k]/((fn + 1.0) * clat);
                }
              else
                {
                  y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
                }
              l = l + 2;
              break;
            case 4:  bb = rr * ghb[l+1];
              cc = aa * cl[m] + bb * sl[m];
              xtemp = xtemp + cc * q[k];
              ztemp = ztemp - cc * p[k];
              if (clat > 0)
                {
                  ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
                    fm * p[k]/((fn + 1.0) * clat);
                }
              else
                {
                  ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
                    q[k] * slat;
                }
              l = l + 2;
              break;
            default: printf("\nError in subroutine shval3");
              break;
            }
        }
      m = m + 1;
    }
  if (iext != 0)
    {
      aa = ext2 * cl[1] + ext3 * sl[1];
      switch(gh)
        {
        case 3:   x = x - ext1 * clat + aa * slat;
          y = y + ext2 * sl[1] - ext3 * cl[1];
          z = z + ext1 * slat + aa * clat;
          break;
        case 4:   xtemp = xtemp - ext1 * clat + aa * slat;
          ytemp = ytemp + ext2 * sl[1] - ext3 * cl[1];
          ztemp = ztemp + ext1 * slat + aa * clat;
          break;
        default:  printf("\nError in subroutine shval3");
          break;
        }
    }
  switch(gh)
    {
    case 3:   aa = x;
		x = x * cd + z * sd;
		z = z * cd - aa * sd;
		break;
    case 4:   aa = xtemp;
		xtemp = xtemp * cd + ztemp * sd;
		ztemp = ztemp * cd - aa * sd;
		break;
    default:  printf("\nError in subroutine shval3");
		break;
    }
  return(ios);
}

//------------------------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/*                           Subroutine dihf                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
/*                                                                          */
/*     Input:                                                               */
/*           x  - northward component                                       */
/*           y  - eastward component                                        */
/*           z  - vertically-downward component                             */
/*                                                                          */
/*     Output:                                                              */
/*           d  - declination                                               */
/*           i  - inclination                                               */
/*           h  - horizontal intensity                                      */
/*           f  - total intensity                                           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 22, 1988                                                */
/*                                                                          */
/****************************************************************************/

int Class_IGRF11::dihf (int gh)
{
  int ios;
  int j;
  double sn;
  double h2;
  double hpx;
  double argument, argument2;
  
  ios = gh;
  sn = 0.0001;
  
  switch(gh)
    {
    case 3:   for (j = 1; j <= 1; ++j)
        {
          h2 = x*x + y*y;
          argument = h2;
          h = sqrt(argument);       /* calculate horizontal intensity */
          argument = h2 + z*z;
          f = sqrt(argument);      /* calculate total intensity */
          if (f < sn)
            {
              d = NaN;        /* If d and i cannot be determined, */
              i = NaN;        /*       set equal to NaN         */
            }
          else
            {
              argument = z;
              argument2 = h;
              i = atan2(argument,argument2);
              if (h < sn)
                {
                  d = NaN;
                }
              else
                {
                  hpx = h + x;
                  if (hpx < sn)
                    {
                      d = PI;
                    }
                  else
                    {
                      argument = y;
                      argument2 = hpx;
                      d = 2.0 * atan2(argument,argument2);
                    }
                }
            }
        }
		break;
    case 4:   for (j = 1; j <= 1; ++j)
        {
          h2 = xtemp*xtemp + ytemp*ytemp;
          argument = h2;
          htemp = sqrt(argument);
          argument = h2 + ztemp*ztemp;
          ftemp = sqrt(argument);
          if (ftemp < sn)
            {
              dtemp = NaN;    /* If d and i cannot be determined, */
              itemp = NaN;    /*       set equal to 999.0         */
            }
          else
            {
              argument = ztemp;
              argument2 = htemp;
              itemp = atan2(argument,argument2);
              if (htemp < sn)
                {
                  dtemp = NaN;
                }
              else
                {
                  hpx = htemp + xtemp;
                  if (hpx < sn)
                    {
                      dtemp = PI;
                    }
                  else
                    {
                      argument = ytemp;
                      argument2 = hpx;
                      dtemp = 2.0 * atan2(argument,argument2);
                    }
                }
            }
        }
		break;
    default:  printf("\nError in subroutine dihf");
		break;
    }
  return(ios);
}
//------------------------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/*                           Subroutine julday                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the decimal day of year from month, day, year.              */
/*     Supplied by Daniel Bergstrom                                         */
/*                                                                          */
/* References:                                                              */
/*                                                                          */
/* 1. Nachum Dershowitz and Edward M. Reingold, Calendrical Calculations,   */
/*    Cambridge University Press, 3rd edition, ISBN 978-0-521-88540-9.      */
/*                                                                          */
/* 2. Claus Tøndering, Frequently Asked Questions about Calendars,          */
/*    Version 2.9, http://www.tondering.dk/claus/calendar.html              */
/*                                                                          */
/****************************************************************************/

double Class_IGRF11::julday(int month,int day, int year)
{
  int days[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

  int leap_year = (((year % 4) == 0) &&
                   (((year % 100) != 0) || ((year % 400) == 0)));

  double day_in_year = (days[month - 1] + day + (month > 2 ? leap_year : 0));

  return ((double)year + (day_in_year / (365.0 + leap_year)));
}


}//namespace IGRF11
