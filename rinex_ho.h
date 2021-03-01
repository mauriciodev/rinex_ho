
/*-------------------------------------------------------------------
    Purpose:  Program to read a RINEX file and apply the second and
              third order ionospheric effects corrections
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             haroldoh2o@gmail.com
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Advisor: João Francisco Galera Monico
             galera@fct.unesp.br
             Departamento de Cartografia
             FCT/UNESP - Presidente Prudente - SP, Brazil

    Date: July of 2010
    -------------------------------------------------------------------
    Observation: Use C++ classes to read and save RINEX file:
                 http://www.ngs.noaa.gov/gps-toolbox/rinex.htm
                 The authors would like to thanks CAPES, FAPESP by financial support
    -------------------------------------------------------------------*/
//c::1304.13, SAH, fix the paths, use just one directory for all *.cpp/*.h files


#include <iostream>
#include <ios>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cctype>
#include <time.h>
//#include <unistd.h>  //must be uncommented for linux compilation

//#include <conio.h> //must be commented for linix compilation

#include <stdio.h>
#include <stdlib.h>
//#include "Rinex_Class/rinex.h"    //class to read and save RINEX (available at: http://www.ngs.noaa.gov/gps-toolbox/rinex.htm)
//#include "Class_TEC/Class_TEC.h"
//#include "Class_Iono/Class_Iono.h"
//#include "Cycle_Slip/CycleSlip.h"
//#include "FilterCode/Filter_Code.h"

#include "rinex.h"    //class to read and save RINEX (available at: http://www.ngs.noaa.gov/gps-toolbox/rinex.htm)
#include "Class_TEC.h"
#include "Class_Iono.h"
#include "CycleSlip.h"
#include "Filter_Code.h"
#include "Comp_Bias.h"

//#include ".\Fortran_Class\fortran.h" //Class to mix C/C++ and Fortran (available at: http://arnholm.org/software/index.htm)

using namespace NGSrinex; //namespace related to the RINEX class
using namespace std;

using namespace CSLIP;
using namespace FILTEROBS;
using namespace CTEC;
using namespace CIONO;




double TECM[NSAT]={0.0};

double i1m_l1[NSAT]={0.0},i1m_l2[NSAT]={0.0},
       i2m_l1[NSAT]={0.0},i2m_l2[NSAT]={0.0},
       i3m_l1[NSAT]={0.0},i3m_l2[NSAT]={0.0};


bool Cycle_Slip_Flag=false;
/*void CorrectIono(int numobstype,ObsEpoch &currentObsEpoch,SatObsAtEpoch tempSatObsAtEpoch[],
                 double I1_L1,double I1_L2,double I2_L1,
                 double I2_L2,double I3_L1,double I3_L2,double ph1[],
                 double ph2[],double ca[],double p1[],double p2[],
                 int pos_currentObsEpoch,int pos_tempSatObsAtEpoch);
 */
void CorrectIono(int numobstype,ObsEpoch &currentObsEpoch,SatObsAtEpoch tempSatObsAtEpoch[],
                 double I1_L1,double I1_L2,double I1_L5,  //not being used
                 double I2_L1,double I2_L2,double I2_L5,
                 double I3_L1,double I3_L2,double I3_L5,
                 double ph1[],double ph2[],double ph5[],   //phase
                 double ca[],double p1[],double p2[],double c2[],double c5[], //code
                 int pos_currentObsEpoch,int pos_tempSatObsAtEpoch);


void ReadInput(string InpFileName,string &filenameobs,string &filenamenav,string &newobsfile,
               string &outfile,double &X0,double &Y0,double &Z0,double &Mask_Ele,
               int &read_tec,int &save_file,char GIM_File[],
               double &br,char arqu_p1c1[],char arqu_p1p2_sat[]);

int ReadNavFile(PRNBlock  currentPRNBlock[],RinexFile *myFileNav,string filenamenav,
                int &cont_efe_nav, fstream &log);

void StoreObsVector(int numobstype,ObsEpoch &currentObsEpoch,double ph1[],
                    double ph2[],double ca[],double p1[],double p2[],
                    double c2[],double c5[],double ph5[],int n_epoch);

void ReadParam(double &sigca,double &sigp2,double &sigph1,
               double &sigph2,double &h_ion,double &R_e,double &B_eq,double &Ne_max);

void sat_pos_vel(double  tr, double  toe,double  dt, double  a, double  ri0,
	         double  dri, double  dn, double  cm0, double  e, double  w,
	         double  cus, double  cuc, double  crs, double  crc, double  cis,
       	         double  cic, double  w0, double  wd, double  &x1, double  &y1, double  &z1,
                double  &dx1, double  &dy1, double  &dz1);

void azelIPP(double staX, double staY, double staZ,
     double svX, double svY, double svZ,
     double Re, double htIono,
     double *az, double *elv, double *dist,
     double *pplat, double *pplon, double *zp);

int xyz_to_geod ( double a, double finv, double x, double y, double z,
                  double *dlat, double *dlon, double *h);


void Sat_Angle (double lat,double longi,double xs,double ys,double zs,
                double Xest,double Yest,double Zest,double &Ele, double &Az );

void zera_vetor(double V[], int n);
void zera_vetor(int V[], int n);

/*typedef struct Mean_Phase_Code{
  double Last_Diff_Li_Pi;
  double Mean_Diff_Li_Pi;
  int Count_Mean;
}MEAN_PHASE_CODE;

*/
//---------------------------------------------------------------------------
void Sat_Angle (double lat,double longi,double xs,double ys,double zs,
                double Xest,double Yest,double Zest,double &Ele, double &Az )

{/*-------------------------------------------------------------------
     Purpose: to compute elevation angle and azimuth of satellite
              in the Geodetic Local System
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

 double dx,dy,dz,ro,coslat,coslon,sinlat,sinlon,de,dn,du;


  dx = xs - Xest;
  dy = ys - Yest;
  dz = zs - Zest;

  //satellite-receiver distance
  ro = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));

  //Local system centered at receiver position

  coslat = cos(lat);
  coslon = cos(longi);
  sinlat = sin(lat);
  sinlon = sin(longi);

  dn = (-1.0*sinlat*coslon*dx) + (-sinlat*sinlon*dy) + (coslat*dz);
  de = (-sinlon*dx)+(coslon*dy);
  du = (coslat*coslon*dx) + (coslat*sinlon*dy) + (sinlat*dz);

  //Local normalized vector
  dn /= ro;
  de /= ro;
  du /= ro;

  // Em = atan2(U,temp);    //Elevation angle
  Ele = asin(du);
 // Zen = (pi/2.0) - Ele;   //Zenithal angle

  //Azimuth
  Az = atan(de/dn);

  //quadrant analisys
  if(de >=0 && dn <0) Az = PI - fabs(Az);                      //2nd quadrant
   else if(de <= 0 && dn <0) Az = Az + PI;                     //3rd quadrant
         else if(de <= 0 && dn > 0) Az = (2.0*PI) - fabs(Az);  //4th quadrant

  if(Az > (2.0*PI))       //if azimuth greater than 360 degrees
    while(Az > (2.0*PI))
      Az -= (2.0*PI);
  else if (Az < 0)       //if azimuth less than 0
    while (Az < 0)
      Az += (2.0*PI);

return;
}
//---------------------------------------------------------------------------
void zera_vetor(double V[], int n)
{/*-------------------------------------------------------------------
    Purpose: Start array
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1

    Date: July of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

 for(int i=0;i<n;i++) V[i] = 0.0;
}
//---------------------------------------------------------------------------
void zera_vetor(int V[], int n)
{/*-------------------------------------------------------------------
    Purpose: Start array
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1

    Date: July of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

 for(int i=0;i<n;i++) V[i] = 0;
}
//---------------------------------------------------------------------------
void ReadInput(string InpFileName,string &filenameobs,string &filenamenav,string &newobsfile,
               string &outfile,double &X0,double &Y0,double &Z0,double &Mask_Ele,
               int &read_tec,int &save_file,char GIM_File[],
               double &br,char arqu_p1c1[],char arqu_p1p2_sat[],
               int &Geomag_Type,char Name_IGRF_Model[])
{/*-------------------------------------------------------------------
    Purpose: Read rinex_ha.inp input file
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
char temp[500]={"\0"},temp1[500]={"\0"};
 string strtemp;

    ifstream inpfile(InpFileName.c_str());

    if(inpfile.fail() )
    {
        cout<<"It was not possible to open the file: "<<InpFileName.c_str()<<endl;
        getchar();
        exit(-1);
    }//if

    //Observation rinex file name
    inpfile.getline(temp,500);
    sscanf(temp,"%s",temp1);
    filenameobs = temp1;

    //broadcast ephemeris file name
    inpfile.getline(temp,500);
    sscanf(temp,"%s",temp1);
    filenamenav = temp1;

    //new observation file name
    inpfile.getline(temp,500);
    sscanf(temp,"%s",temp1);
    newobsfile = temp1;

    //output files
    inpfile.getline(temp,500);
    sscanf(temp,"%s",temp1);
    outfile = temp1;

     inpfile.getline(temp,500);  //Receiver coordinates (if 0.0 => Try to read coordinates from RINEX)
    sscanf(temp,"%lf %lf %lf",&X0,&Y0,&Z0);

    //Elevation Mask (means that observables under this mask won't be corrected)
    inpfile.getline(temp,500);
    sscanf(temp,"%lf",&Mask_Ele);

    //save files (yes =1; no =0)
    inpfile.getline(temp,500);
    sscanf(temp,"%d",&save_file);


    //0 = Tec from raw pseudorange; 1 = TEC from smothed pseudorange by phase; 2 = TEC from GIM
    inpfile.getline(temp,500);
    sscanf(temp,"%d",&read_tec);

    //ionex GIM file name. If not using GIM, please insert any name, just to read
    inpfile.getline(temp,500);
    sscanf(temp,"%s",GIM_File);

    //receiver DCB (P1-P2)
    inpfile.getline(temp,500);
    sscanf(temp,"%lf",&br);

    //DCB (P1-C1) file name - default from CODE: P1C1.DCB
    inpfile.getline(temp,500);
    sscanf(temp,"%s",arqu_p1c1);

    //DCB (P1-P2) file name - default from CODE: DCBSAT.DCB
    inpfile.getline(temp,500);
    sscanf(temp,"%s",arqu_p1p2_sat);

    inpfile.getline(temp,500);  //0 = Dipolar model; 1 = CGM from PIM; 2 = IGRF model
    sscanf(temp,"%d",&Geomag_Type);

    inpfile.getline(temp,500);
    sscanf(temp,"%s",Name_IGRF_Model); //Name of the IGRF coefficients

return;
}
//---------------------------------------------------------------------------
void ReadParam(double &sigca,double &sigp2,double &sigph1,
               double &sigph2,double &h_ion,double &R_e,double &B_eq,double &Ne_max)
{/*-------------------------------------------------------------------
    Purpose: Read parameters from rinex_ha_param.dat
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
char temp[500]={"\0"};


    ifstream fileparam ("rinex_ha_param.dat");

    if(fileparam.fail() )
    {
        cout<<"It was not possible to open the rinex_ha_param.dat file"<<endl;
        getchar();
        exit(-1);
    }//if

    //CA precision
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&sigca);

    //P2 precision
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&sigp2);

    //phase L1 precision
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&sigph1);

    //phase L2 precision
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&sigph2);

    //ionospheric layer height (m)
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&h_ion);

    //Equatorial Earth axis (m)
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&R_e);

    //geomagneitc induction magnitude at equator (Tesla)
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&B_eq);

    //Electron density maximum
    fileparam.getline(temp,500);
    sscanf(temp,"%lf",&Ne_max);

return;
}
//---------------------------------------------------------------------------
void CorrectIono(int numobstype,ObsEpoch &currentObsEpoch,SatObsAtEpoch tempSatObsAtEpoch[],
                 double I1_L1,double I1_L2,double I1_L5,  //not being used
                 double I2_L1,double I2_L2,double I2_L5,
                 double I3_L1,double I3_L2,double I3_L5,
                 double ph1[],double ph2[],double ph5[],   //phase
                 double ca[],double p1[],double p2[],double c2[],double c5[], //code
                 int pos_currentObsEpoch,int pos_tempSatObsAtEpoch)
{/*-------------------------------------------------------------------
    Purpose: Correct GPS observable from 2nd and 3rd order ionospheric effects
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
    Observation: Corrections for L5 must be implemented
    -------------------------------------------------------------------*/

    double corrected_obs(0);
    double f1 = 1575420000.0;   //L1 frequency - Hz
    double f2 = 1227600000.0;   //L2 frequency - Hz
    double f5 = 1176450000.0;    //L5 frequency - Hz

    double c = 299792458.0;     //velocity of light
    double lamb_l1 = c/f1;      //wavelenght - L1
    double lamb_l2 = c/f2;      //wavelenght - L2
    double lamb_l5 = c/f5;      //wavelenght - L2

    for( int j = 0; j < numobstype; j++ )
    {
            //name associated to each type of observation
            int var = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;


            switch(var)
            {
                /* In the rinex.h file you can find enum variable (It was changed from original version by Marques, H. A.)
                enum OBSTYPE { NOOBS = 0, L1 = 1, L2 = 2, C1 = 3, P1 = 4, P2 = 5,
                               C2 = 6, D1 = 7, D2 = 8, T1 = 9, T2= 10, S1 = 11, S2 = 12,
                               L5 =13, C5 = 14, D5 = 15, S5 = 16  };
              */

                case 1: //Phase - L1
                        if(ph1[pos_tempSatObsAtEpoch]!=0.0)
                           corrected_obs = ph1[pos_tempSatObsAtEpoch] + ((I2_L1/2.0) + (I3_L1/3.0))/lamb_l1;
                        else corrected_obs = ph1[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 2: //Phase - L2
                        if(ph2[pos_tempSatObsAtEpoch]!=0.0)
                          corrected_obs = ph2[pos_tempSatObsAtEpoch] + ((I2_L2/2.0) + (I3_L2/3.0))/lamb_l2;
                        else corrected_obs = ph2[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 3: //CA
                        if(ca[pos_tempSatObsAtEpoch]!=0.0)
                           corrected_obs = ca[pos_tempSatObsAtEpoch] - I2_L1 - I3_L1;
                        else corrected_obs = ca[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 4: //P1
                        if(p1[pos_tempSatObsAtEpoch]!=0.0)
                          corrected_obs = p1[pos_tempSatObsAtEpoch] - I2_L1 - I3_L1;
                        else corrected_obs = p1[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 5: //P2
                        if(p2[pos_tempSatObsAtEpoch]!=0.0)
                          corrected_obs = p2[pos_tempSatObsAtEpoch] - I2_L2 - I3_L2;
                        else corrected_obs = p2[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 6:  //C2
                        if(c2[pos_tempSatObsAtEpoch]!=0.0)
                          corrected_obs = c2[pos_tempSatObsAtEpoch] - I2_L2 - I3_L2;
                        else corrected_obs = c2[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 7: //D1
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 8: //D2
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 9: //T1
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 10://T2
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
               case 11: //S1
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 12: //S2
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                case 13:  //Phase - L5
                        if(ph5[pos_tempSatObsAtEpoch]!=0.0)
                        {  printf("To be imlemented for L5\n");
                         //corrected_obs = ph5[pos_tempSatObsAtEpoch] + ((I2_L5/2.0) + (I3_L5/3.0))/lamb_l5;
                        }
                        else corrected_obs = ph5[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;

                        break;
                case 14:  //C5
                        if(c5[pos_tempSatObsAtEpoch]!=0.0)
                        {  printf("To be imlemented for L5\n");
                          // corrected_obs = c5[pos_tempSatObsAtEpoch] - I2_L5 - I3_L5;
                        }
                        else corrected_obs = c5[pos_tempSatObsAtEpoch];

                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = corrected_obs;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                  case 15: //D5
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;
                   case 16: //S5
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].observation = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].observation;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsPresent = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsPresent;;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].obsType = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].obsType;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satCode = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satCode;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].satNum = currentObsEpoch.getSatListElement(pos_currentObsEpoch).satNum;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].LLI = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].LLI;
                        tempSatObsAtEpoch[pos_tempSatObsAtEpoch].obsList[j].sigStrength = currentObsEpoch.getSatListElement(pos_currentObsEpoch).obsList[j].sigStrength;
                        break;

                default: break;
            }//switch
       // }//if
    }//for

    return;
}
//---------------------------------------------------------------------------
void StoreObsVector(int numobstype,ObsEpoch &currentObsEpoch,double ph1[],
                    double ph2[],double ca[],double p1[],double p2[],
                    double c2[],double c5[],double ph5[],int n_epoch)
{/*-------------------------------------------------------------------
    Purpose: store observables from currentObsEpoch in arrays
    -------------------------------------------------------------------
    Input: numobstype - number of obs type
           currentObsEpoch - object from class ObsEpoch
           n_epoch - number of epoch to be extracted data

    Output: ph1, ph2, ca, p2 - arrays with observables
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1

    Date: July of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

    for( int j = 0; j < numobstype; j++ )
    {
        if ( j != 0 &&  j % 5 == 0 )  // Max. 5 obs. per line
        {
            cout << endl;
        }
       if ( currentObsEpoch.getSatListElement(n_epoch).obsList[ j ].obsPresent ) //Check if available observable
        {
            int var = currentObsEpoch.getSatListElement(n_epoch).obsList[j].obsType;

            switch(var)
            { /* In the rinex.h file you can find enum variable
                enum OBSTYPE { NOOBS = 0, L1 = 1, L2 = 2, C1 = 3, P1 = 4, P2 = 5,
                               C2 = 6, D1 = 7, D2 = 8, T1 = 9, T2= 10, S1 = 11, S2 = 12,
                               L5 =13, C5 = 14, D5 = 15, S5 = 16  };
              */
                case 1: ph1[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;   break;
                case 2: ph2[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;   break;
                case 3: ca[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
                case 4: p1[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
                case 5: p2[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
                case 6: c2[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;

               // case 7: dopp1[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
               // case 8: dopp2[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
               // case 11: snr1[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
               // case 12: snr2[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;

                case 13: ph5[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;
                case 14: c5[n_epoch] = currentObsEpoch.getSatListElement(n_epoch).obsList[j].observation;    break;

                default: break;




            }//switch
         }//if
    }//for
}
//---------------------------------------------------------------------------
int ReadNavFile(PRNBlock  currentPRNBlock[],RinexFile *myFileNav,string filenamenav,
                int &cont_efe_nav, fstream &log)
{/*-------------------------------------------------------------------
    Purpose: Read broadcast orbit file
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
string theNavFileType;

    try
    {   //opening file
        myFileNav->setPathFilenameMode( filenamenav, ios::in );

        //type of file
        theNavFileType = myFileNav->getRinexFileType();
        delete myFileNav;
        myFileNav = 0;
    }//try
    catch( RinexFileException &openExcep )
    {
        cout << "Error opening file: " << filenamenav << endl
             << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
             << "Exiting program... " << endl << endl;

        log << "Error opening file: " << filenamenav << endl
            << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
            << "Exiting program..." << endl << endl;

        exit (-1);
    }//catch


    //If the file is an NAV file
    if( theNavFileType[0] == 'N' )
    {
        RinexNavFile mynav;

        try
        {
            mynav.setPathFilenameMode(filenamenav,ios_base::in);
        }
        catch( RinexFileException &openExcep )
        {
            cout << "Error opening file: " << filenamenav << endl;
            cout << " RinexFileException is: " << endl <<
            openExcep.getMessage() << endl;
        }//catch

        try
        {
            mynav.readHeader();

        }
        catch( RequiredRecordMissingException &headerExcep )
        {
            cout << " RequiredRecordMissingException is: " << endl <<
            headerExcep.getMessage() << endl;
        }

        //read the PRN Blocks
        try
        {   cont_efe_nav=0; //number of ephemeris block

            while( mynav.readPRNBlock( currentPRNBlock[cont_efe_nav] ) != 0 )
            {
                cont_efe_nav++;
            }//while
        }//try
        catch( RinexReadingException &readingExcep )
        {
            cout << " RinexReadingException is: " << endl <<
            readingExcep.getMessage() << endl;
        }//catch

        log << endl << endl;
        log << "Number of OBS with WARNINGS = " << mynav.getNumberWarnings() << endl;
        log << endl << endl;
        log << "*** Messages ERRORS:" << endl;
        log << mynav.getErrorMessages() << endl;
        log << "*** Messages WARNINGS:" << endl;
        log << mynav.getWarningMessages() << endl;


    }//if

return (1);

}
//---------------------------------------------------------------------------
void sat_pos_vel( double  tr, double  toe,double  dt, double  a, double  ri0,
	              double  dri, double  dn, double  cm0, double  e, double  w,
	              double  cus, double  cuc, double  crs, double  crc, double  cis,
       	          double  cic, double  w0, double  wd, double  &x1, double  &y1, double  &z1,
                  double  &dx1, double  &dy1, double  &dz1)
{/*-------------------------------------------------------------------
    Purpose: Compute the satellite position and velocity at tr time
    -------------------------------------------------------------------
    Input: toe --> TIME ORIGIN EPHEMERIS:
           a   --> SQUARE ROOT SEMI MAJOR AXIS AT TOE:
           ri0 --> INCLINATION AT TOE:
           dri --> RATE  OF INCLINATION:
           dn  --> CORRECTION  OF MEAN MOTION:
           cm0 --> MEAN ANOMALY AT TOE:
           e   --> ECCENTRICITY:
           w   --> ARGUMENT OF PERIGEE:
           cus --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE ARGUMENT OF LATITUDE:
           cuc --> IBID COSINE:
           crs --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE ORBITS RADIUS:
           crc --> IBID COSINE:
           cis --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE INCLINATION ANGLE:
           cic --> IBID COSINE:
           w0  --> RIGHT ASCENSION AT REFERENCE TIME:
           wd  --> RATE OF RIGHT ASCENCION:

    Output: x1, y1, z1 - satellite coordinates
            dx1, dy1, dz1 - satellite velocity
    -------------------------------------------------------------------
    Authors: Joao Francisco Galera Monico

             Adapted and converted to C by Marques, H. A. 2007

             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1

    Date: July of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

   int i=0,j=0;

   double argl, arglc, arglc1, arglcc, argls, argls1, arglsc, b,
	  cargl, cicw, cinc, cisw, cmk, cn, cno, cosi, cosw, crad, cta,
	  dargl, deca, den, dinc, drad, drad1,dta, dw, dxop, dyop,
	  ec[10], eca, pi, sicw, sini, sinw, sisw, sta, sta1, tang, tta,
	  ww, xop, yop;

    const double u = 3.986004418E+14;     // Gravitational Constant - m/s2
    const double we = 7.2921151467E-5;    //Earth rotation rate

    //time diff from toe
    pi = 4.0*atan(1.0);
    dt = 0.0;
    dt = tr - toe;

    //Accounting for begning or end of week crossovers
    if(dt > 302400.0)
      dt -= 604800.0;
    else if (dt < -302400)
      dt += 604800.0;

    //semi-minor axis
    b = a;
    b = pow(b,2);

    //corrected mean anomaly
    cno = sqrt(u/pow(b,3));
    cn = cno + dn;
    cmk = cm0 + cn * (dt);

    //excentry anomaly and latitude argument
    for (j=0;j<10;j++) ec[j]=0.0;

    i = 1;
    ec[0] = cmk;

    do{ ec[i] = cmk + e * sin(ec[i-1]);
        i++;
      }while( ((fabs(ec[i]-ec[i-1])) <= 0.000000000001)||(i<10));

    eca = ec[i-1];
    den = 1.0 - e * cos(eca);
    cta = (cos(eca) - e)/den;
    sta1 = sqrt( fabs(1.0- pow(e,2)) );
    sta = sta1 * sin(eca)/den;

    if (cta!=0) tta = sta/cta;

    //quadrant
    if( (sta >= 0.0) && (cta >= 0.0) )
         tang = atan(tta);
   	else if( (sta >= 0.0) && (cta <= 0.0) )
            tang = pi + atan(tta);
	     else if( (sta < 0.0) && (cta < 0.0) )
		      tang = pi + atan(tta);
	     else //if( (sta < 0.0) && (cta >= 0.0) )
		      tang = 2 * pi + atan(tta);

    argl = tang + w;

    //excentry anomaly derivative
    deca = cn/den;

    //true anomaly derivative
    dta = (sta1 * deca)/den;

    //latitude argument correction, ray and inclination
    argls = sin(argl*2.e0);
    arglc = cos(argl*2.e0);
    cargl = cus*argls + cuc*arglc;
    crad = crc*arglc + crs*argls;
    cinc = cic*arglc + cis*argls;
    cargl += argl;
    crad += b * den;
    cinc += ri0 + dri * (dt);

    //latitude argument derivative
    argls1 = argls*dta;
    arglc1 = arglc*dta;
    dargl = dta - 2.0 * cuc * argls1 + 2.0 * cus * arglc1;

    //ray derivative
    drad1 = b*e*sin(eca)*deca;
    drad = drad1 - 2.0*crc*argls1 + 2.0*crs*arglc1;

    //inclination derivative
    dinc = dri - 2.0*cic*argls1 + 2.0*cis*arglc1;

    //orbital plan coordinates
    arglsc = sin(cargl);
    arglcc = cos(cargl);
    xop = crad*arglcc;
    yop = crad*arglsc;

    //orbital plan velocity
    dxop = drad*arglcc - crad*dargl*arglsc;
    dyop = drad*arglsc + crad*dargl*arglcc;

    //latitude correction of ascendent node
    dw = wd - we;
    ww = w0 + dw * dt - we* toe;

    //terrestrial coordinates - WGS84
    sinw = sin(ww);
    cosw = cos(ww);
    cosi = cos(cinc);
    sini = sin(cinc);
    cicw = cosi*cosw;
    cisw = cosi*sinw;
    sisw = sini*sinw;
    sicw = sini*cosw;

    x1 = ((xop * cosw) - (yop * cisw));
    y1 = ((xop * sinw) + (yop * cicw));
    z1 = (yop * sini);

    //velocities - WGS84
    dx1 = cosw*dxop - cisw*dyop - sinw*dw*xop + (-cicw*dw + sisw*dinc)* yop;
    dy1 = sinw*dxop + cicw*dyop + cosw*dw*xop + (-cisw*dw - sicw*dinc)* yop;
    dz1 = sini*dyop + cosi*dinc*yop;
    dt =  toe - tr;

return;
} //satpove
//---------------------------------------------------------------------------

void azelIPP(double staX, double staY, double staZ,
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

   du = -1.0*srlat*crlon*dx - srlat*srlon*dy + crlat*dz;
   dv = -1.0*srlon*dx + crlon*dy;
   dw = crlat*crlon*dx + crlat*srlon*dy + srlat*dz;

   *dist = sqrt( du*du + dv*dv + dw*dw );

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

// cout << "^^^^^ from azelIPP ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
// cout << "PP lat  (deg)" << setw(30) << setprecision(15) << geod_pp_lat << endl;
// cout << "PP long (deg)" << setw(30) << setprecision(15) << geod_pp_lon << endl;
// cout << "PP hgt  (m)  " << setw(30) << setprecision(5)  << geod_pp_hgt << endl;
// cout << "azimuth (rad)" << setw(30) << setprecision(15)  << *az << endl;
// cout << "azimuth (deg)" << setw(30) << setprecision(15)  << *az * (180.0/ppii) << endl;
// cout << "elev    (deg)" << setw(30) << setprecision(15)  << *elv * (180.0/ppii) << endl;
// cout << "SV xyz  (m)" << setw(16) << setprecision(4) << svX << " "
//                       << setw(16) << setprecision(4) << svY << " "
//                       << setw(16) << setprecision(4) << svZ << endl;
// cout << "zp,htIon,Re" << setw(16) << setprecision(9) << *zp * (180.0/ppii) << " "
//                       << setw(16) << setprecision(4) << htIono << " "
//                       << setw(16) << setprecision(4) << Re  << endl;
// cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

   // compute geocentric latitude and longitude for IPP
   dequa = sqrt( xipp*xipp + yipp*yipp);
   *pplat = atan2( zipp, dequa);
   *pplon = atan2( yipp, xipp);

   *pplat = *pplat * (180.0/ppii);  // convert rad to deg
   *pplon = *pplon * (180.0/ppii);
}

/*---------------------------------------------------------------------*/

 int xyz_to_geod ( double a, double finv, double x, double y, double z,
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


