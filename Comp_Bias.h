/*-------------------------------------------------------------------
    Purpose:  To compute phase bias for all satellite arcs data in the rinex file
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques

    Date: Sept. of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
#ifndef COMPBIAS
#define COMPBIAS

#include <iomanip>
#include <cmath>
#include <cstring>
#include <cctype>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
//#include<dir.h>
#include <sys/stat.h>

#include "CycleSlip.h"
#include "Class_TEC.h"

#include "rinex.h"


using namespace std;

using namespace NGSrinex; //namespace related to the RINEX class
using namespace CTEC;
using namespace CSLIP;

//Struct to be used in each epoch
typedef struct Mean_Phase_Code{
  double Last_Diff_Li_Pi;
  double Mean_Diff_Li_Pi;   //mean of phase bias coputed recursevely
  double Start_Time;
  double End_Time;
  int Count_Mean;
  int arc_number;
}MEAN_PHASE_CODE;

//struct to store data for each arc
typedef struct Phase_Bias{
  double Mean_Phase_Bias;  //final mean phase bias for the arc
  double Start_Time;
  double End_Time;
  int Count_Mean;         //number of data in each arc
  int arc_number;         //number of the arc
}MEAN_PHASE_BIAS;

const int max_arc = 100;




//functions
void Comp_Bias( string filenameobs, Class_TEC Calc_Tec,MEAN_PHASE_BIAS Phase_Bias[][max_arc],
                ofstream &File_Cycleslip);

void Close_PRN_Files();

bool directory_exists(string pathname);

int Create_Dir(string Directory);

double Get_Phase_Bias(int PRN, MEAN_PHASE_BIAS Phase_Bias[][max_arc], double MJD);

void Check_Sat_Out( int NS,
                    int isat[],
                    bool Sat_Epoch_Ant[],
                    bool Sat_History[],
                    bool SatOut[],
                    bool &flag_sat_out);

void Check_New_Sat( int NS, int isat[],
                    bool Sat_Epoch[],
                    bool Sat_History[],
                    bool New_Sat[],
                    bool &flag_sat_in);

int OpenPRNFile();
void Print_Data_PRN(int PRN, MEAN_PHASE_CODE *Sat_Average_Phase_Code);


#endif
