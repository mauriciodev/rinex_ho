/*-------------------------------------------------------------------
    Purpose: To compute Pseudorange smothed by phase
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             FAPESP PROCESS: 05/03522-1
    Date: July of 2010
    Update by Marques, H. A. in May of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/

//c::1304.13, SAH, change path, use one directory for all *.cpp/*.h files


#ifndef FCode
#define FCode

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
//#include "../Include/constant.h"
#include "constant.h"

using namespace std;

namespace FILTEROBS{

class FilterCode{

    public:

             FilterCode(); //constructor

             void Start_Smooth_Code(double TGPS,double P1,double P2, double PH1,
                                      double PH2,double var_p1,double var_p2,int iono);

             void Filter_Obs(double TGPS,double P1,double P2,double PH1,double PH2,double var_p1,
                             double var_p2,double var_ph1,double var_ph2,int iono );

             //Get methods
             double Get_ph1_bef();
             double Get_ph2_bef();
             double Get_p1_filt_bef();
             double Get_p2_filt_bef();
             double Get_p1_pred();
             double Get_p2_pred();
             double Get_p1_filt();
             double Get_p2_filt();
             double Get_var_p1_filt();
             double Get_var_p2_filt();

             int Get_quality();

             int Get_cont_epoch();

             //set Methods
             void Set_ph1_bef(double value);
             void Set_ph2_bef(double value);
             void Set_p1_filt_bef(double value);
             void Set_p2_filt_bef(double value);
             void Set_p1_pred(double value);
             void Set_p2_pred(double value);
             void Set_p1_filt(double value);
             void Set_p2_filt(double value);
             void Set_var_p1_filt(double value);
             void Set_var_p2_filt(double value);

             void Set_quality(int value);

             void Set_cont_epoch(int value);
             void Set_var_if_code(double value);
             void Set_var_if_phase(double value);

    private:
             void Prediction(double p1_k,double p2_k,double ph1_k,double ph2_k,int iono);

             void Smoothing(double p1,double p2,double ph1,double ph2,double var_p1,double var_p2,
                            double var_ph1,double var_ph2,int iono);

             double ph1_bef,
                    ph2_bef,
                    p1_filt_bef,
                    p2_filt_bef,
                    p1_pred,
                    p2_pred,
                    p1_filt,
                    p2_filt,
                    var_p1_filt,
                    var_p2_filt;

             double  PR_Filt_Bef,
                     PR_IF,PR_IF_Pred,PR_IF_Filt_Bef;

             double var_if_code,var_if_phase;

             double dtime;
             int quality;

             int cont_epoch;


};


}//namespace

#endif
