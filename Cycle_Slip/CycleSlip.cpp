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

#include "CycleSlip.h"

namespace CSLIP{

//---------------------------------------------------------------------------
Class_Cycle_Slip::Class_Cycle_Slip()
{
 Start_CycleSlip();


 ICslip = 0;

 old_ph1 = 0.0;
 old_ph2 = 0.0;
 old_p1  = 0.0;
 old_p2  = 0.0;

}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Start_CycleSlip()
{
   wlbias = 0.0;
   diff_wl = 0.0;
   sig_wlbias = 0.0;
   flag = false;
   count = 0;
   ICslip = 0;

}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Wide_Cycle_Slip(double ph1,double ph2,double ca,double p2,bool initial_epoch)
{/*-------------------------------------------------------------------
     Purpose: To compute cycle slip based on wide-lane combination
     -------------------------------------------------------------------
     Input:  CycleSlip - struct of WLCycleSlip type
             ph1 - L1 phasse in cycles
             ph2 - L2 phase in cycles
             ca  - CA in meters
             p2  - P2 in meters
             initial_epoch - to indicate initial epoch of the data arc

     Output:  bias and sigma wide-lane inside CycleSlip struct
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              FAPESP PROCESS: 05/03522-1
     Date: July of 2010
     -------------------------------------------------------------------
     Observation: We are following Turbo Edit algorithm:
                  Blewitt, G. (1990), An Automatic Editing Algorithm for GPS data, Geophys. Res. Lett., 17(3),
                  199–202, doi:10.1029/GL017i003p00199
     -------------------------------------------------------------------*/

  double Lw=0.0,Pw=0.0;
  double bwi=0.0,wlbias_before=0.0,sig_before=0.0;

    Lw = wl1p*ph1 + wl2p*ph2;   // Wide lane phase (meters)
    Pw = wl1r*ca + wl2r*p2;     // 'NL' range (meters)

    bwi = (1.0/wlwl)*(Lw - Pw); //time-average phase widelane - code widelane

    if (initial_epoch || count == 0.0) //first epoch or starting new computation
    {
      wlbias = bwi;
      sig_wlbias = pow(0.5*wlwl,2); //set intial sigma as 0.5 wide-lane cycles
      count=1;

    }//if
    else
    {
      count++;  //current number of points in the data arc
    }//else

    if(flag == true) //if there is cycle slip in the epoch k-1
    {  sig_wlbias = pow(0.5*wlwl,2); //set intial sigma as 0.5 wide-lane cycles
       count = 1;
       flag = false;

    }//if

    if(!initial_epoch && count >= 1 )
    {

        wlbias_before = wlbias; //bias from epoch k-1

        diff_wl = bwi - wlbias_before;

        wlbias = wlbias_before + (diff_wl/count); //filter for wide-lane

        sig_before =  sig_wlbias; //variance from epoch k-1

        //computation of recursive mean
        if(diff_wl !=0 )
          sig_wlbias = sig_before+(1.0/(count))*( pow(diff_wl,2) - sig_before);


        //Subsequent epoch estimates bw(i+1) are required to lie within 4sig(i)
        //if diff_wl greather than 4.0*sig from epoch k-1
        if( fabs(diff_wl) > (4.0*sqrt(sig_before)) )
        {
           flag = true;  //to indicate a cycle slip
          //cout<<"prn: "<<prn<<" sig_wlbias: "<<sig_wlbias<<endl;

        }//if
    }//if

    //WN_Cycle_Slip(ph1,ph2,ca,p2,initial_epoch);

  return;

}
//---------------------------------------------------------------------------
void  Class_Cycle_Slip::WN_Cycle_Slip(double ph1,double ph2,double p1,double p2,bool initial_epoch)
{

  double Curnl=0.0, Oldnl=0.0,Difnl=0.0, Difwl=0.0;

  double DMAXWL = 2.0,  //2 m
         DMAXNL = 0.2;  //0.2
  /*
    if(flag == true) //if there is cycle slip in the epoch k-1
    {
      count = 0.0;
      flag = false;

    }//if

    if (initial_epoch || count == 0.0) //first epoch or starting new computation
    {
      count=1;

    }//if
    else
    {
      count++;  //current number of points in the data arc
    }//else
   */
    if(!initial_epoch && count >= 1
       && old_ph1 && old_ph2 && old_p1 && old_p2)
    {
       Curnl = ( ph2*wlf2 - ph1*wlf1);
       Oldnl = ( old_ph2*wlf2 - old_ph1*wlf1);

       if(Oldnl!=0)
         Difnl = Curnl - Oldnl;

       Difwl = ( ph1 - ph2 )*wlwl - ( old_ph1 - old_ph2 )*wlwl
             - ( p1 + p2 )/2.0 + ( old_p1 + old_p2 )/2.0;

       // if TOLERANCE FAILURE
       if( fabs(Difnl) > DMAXNL )
        flag=true;

       if(fabs(Difwl) > DMAXWL)
         flag=true;

     }//if

     old_ph1 = ph1;
     old_ph2 = ph2;
     old_p1  = p1;
     old_p2  = p2;

  return;
}
//---------------------------------------------------------------------------
double Class_Cycle_Slip::Get_wlbias()        //wide lane bias
{
     return (wlbias);
}
//---------------------------------------------------------------------------

double Class_Cycle_Slip::Get_sig_wlbias()    //RMS wide lane bias
{
     return (sig_wlbias);
}
//---------------------------------------------------------------------------
double Class_Cycle_Slip::Get_diff_wl()       //diff betwen phase wide-lane and code wide-lane
{
     return (diff_wl);
}
//---------------------------------------------------------------------------
bool Class_Cycle_Slip::Get_flag()            //flag to indicate a potential cycle slip
{
     return (flag);
}
//---------------------------------------------------------------------------
int Class_Cycle_Slip::Get_count()            //current point in the data arc
{
     return (count);
}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Set_wlbias(double value)        //wide lane bias
{
 wlbias = value;
}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Set_sig_wlbias(double value)    //RMS wide lane bias
{
 sig_wlbias = value;
}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Set_diff_wl(double value)       //diff betwen phase wide-lane and code wide-lane
{
  diff_wl = value;
}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Set_flag(bool value)            //flag to indicate a potential cycle slip
{
   flag = value;
}
//---------------------------------------------------------------------------
void Class_Cycle_Slip::Set_count(double value)            //current point in the data arc
{
  count = value;
}
//---------------------------------------------------------------------------

}//namespace CSlip
