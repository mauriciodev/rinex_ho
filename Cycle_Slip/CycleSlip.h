/*-------------------------------------------------------------------
    Purpose: Class to detect cycle slip
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

#ifndef Cycle
#define Cycle

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "../Include/constant.h"

using namespace std;

using namespace std;

namespace CSLIP{

class Class_Cycle_Slip{

 public:

    int ICslip;
    Class_Cycle_Slip(); //constructor

	void Wide_Cycle_Slip(double ph1,double ph2,double ca, double p2, bool initial_epoch);

    void WN_Cycle_Slip(double ph1,double ph2,double p1,double p2,bool initial_epoch);

    void Start_CycleSlip();

    //Get Method
    double Get_wlbias();        //wide lane bias
    double Get_sig_wlbias();    //RMS wide lane bias
    double Get_diff_wl();       //diff betwen phase wide-lane and code wide-lane
    bool Get_flag();            //flag to indicate a potential cycle slip
    int Get_count();            //current point in the data arc

    //Set Method
    void Set_wlbias(double value);
    void Set_sig_wlbias(double value);
    void Set_diff_wl(double value);
    void Set_flag(bool value);
    void Set_count(double value);

  private:
    double wlbias;        //wide lane bias
    double sig_wlbias;    //RMS wide lane bias
    double diff_wl;       //diff betwen phase wide-lane and code wide-lane
    bool flag;            //flag to indicate a potential cycle slip
    int count;            //current point in the data arc

    double old_ph1, old_ph2, old_p1, old_p2;

};

}//namespace CSlip

#endif
