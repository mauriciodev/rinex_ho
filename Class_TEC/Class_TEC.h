/*-------------------------------------------------------------------
    Purpose: Class to Read GIM and interpoalte TEC,
             to compute TEC from raw pseudorange or from
             pseudorange smoothed by carrier phase
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP - Brazil

    Date: January of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

#ifndef CLASSTEC_H
#define CLASSTEC_H
#include <iostream>
#include <ios>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cctype>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "../Include/constant.h"

using namespace std;


namespace CTEC{
class Class_TEC
{
	public:


		Class_TEC(); //constructor

        //constructor
        Class_TEC( double sigca,double sigp2,
                   double sigph1,double sigph2,
                   double dcbr,double sdcbr);

		~Class_TEC(); //destructor

        //Read GIM file
        int ReadGIM(char *NameGIM); 

        //Read P1-C1 and P1-P2 files
        void Read_DCB_Sat(char arqu_p1c1[],char arqu_p1p2[]); 

        //VTEC from GIM
        bool CalcTecEpoch(double mjd,double fracOfDay,
                           double latuser,double longuser);

        //Linear interpolation
        void Linear_Interp(double VTEC1,double VTEC2,
                           double coord1,double coord2,
                           double coord_user,double &VTEC );

        void Pierce_Point_Coord( double Xest,double Yest,double Zest,
                                 double Xsat,double Ysat,double Zsat,
                                 double lat,double lamb, double h);

         void Sist_Local(double xest,double yest,double zest,
                          double xsat,double ysat,double zsat,
                          double fi, double lambda,
                          double &e,double &n,double &u);

         void Ang_Local(double E, double N, double U,
                         double Xsat,double Ysat,double Zsat,
                         double &Am, double &Em, double &Zm);

         //STEC from pseudorange
         double Calc_STEC_PR(double CA,double P2,int prn);


         void Set_var_ca(double value);
         void Set_var_p2(double value);          //variances for CA and P2
         void Set_var_ph1(double value);
         void Set_var_ph2(double value);         //variances for phase L1 and phase L2
         void Set_br(double value);              //receiver hardware delay
         void Set_sbr(double value);

         double Get_P1P2(int prn);
         double Get_P1C1(int prn);
         double Get_S_TEC(int prn);
         double Get_br();                        //receiver hardware delay
         double Get_sbr();

         double Get_LatIon();
         double Get_LambIon();
         double Get_ZL();
         double GetVTEC();
         int Get_Gim_Flag();

         double Get_Epoch_First_MAP();
         double Get_Epoch_Last_Map();
         double Get_Interval();
         double Get_Ini_Lat();
         double Get_End_Lat();
         double Get_Lat_Increment();
         double Get_Ini_Lon();
         double Get_End_Lon();
         double Get_Lon_Increment();
         double Get_Elevation_Cuttof();
         double Get_Base_Radius();
         double Get_Ini_Height();
         double Get_End_Height();
         double Get_Height_Incremet();

    private:


         void ReadMap(int numval,int iexp, double tecval[],int &irc,ifstream &f1);
         double DJUL(int j1,int m1, double t);

         void ReadHeaderGim(ifstream &f1, int &IEXP, int &IRC);

         double **TECMAP;           //Matrix to store Ionospheric Maps;
         double Az,El,Zen,zl;
         double lat_ion,lamb_ion;
         double VTEC;

         //To use in the reading GIM file
         int Num_Value;             //number of Maps read in the GIM file
         int Gim_Flag;              //GIM File Flag

         double Epoch_First_MAP;    //in MJD
         double Epoch_Last_Map;     //in MJD
         double Interval;           //Seconds
         double Ini_Lat;            //From Latitude
         double End_Lat;            //To Latitude
         double Lat_Increment;      //in Degrees
         double Ini_Lon;            //From Longitude
         double End_Lon;            //To Longitude
         double Lon_Increment;      //in Degrees
         double Elevation_Cuttof;   //in Degrees
         double Base_Radius;        //in Km
         double Ini_Height;         //From height
         double End_Height;         //to Height
         double Height_Incremet;    //in km
         double hion;               //Ionospheric layer height in m
         double Re;                 //Equatorial Earth axis in m

         //Used for STEC from pseudorange
         double var_ca,var_p2,      //variances for CA and P2
                var_ph1,var_ph2;    //variances for phase L1 and phase L2

         //interfrequency bias to be read in th input files
         double p1p2[NSAT][2],      //Satellite P1-P2 hardware delay and standard deviation (ns)
                p1c1[NSAT][2],      //Satellite P1-C1 hardware delay
                br,                 //receiver hardware delay
                sbr;                //standard deviation of receiver hardware delay

         double s_tec[NSAT];        //Standard Deviation of TEC


};

}//namespace

#endif // CLASSGIM_H
