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
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
#include "Filter_Code.h"

namespace FILTEROBS{

FilterCode::FilterCode()//Constructor
{ /*-------------------------------------------------------------------
      Purpose: Initializa variables
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
    ph1_bef = 0.0;
    ph2_bef = 0.0;
    p1_filt_bef = 0.0;
    p2_filt_bef = 0.0;

    p1_pred = 0.0,
    p2_pred = 0.0,
    p1_filt = 0.0;
    p2_filt = 0.0;

    var_p1_filt = 0.0;
    var_p2_filt = 0.0;

    quality = 1;

    cont_epoch=0;


}
//------------------------------------------------------------------------------------
void FilterCode::Filter_Obs(double TGPS,double P1,double P2,double PH1,double PH2,double var_p1,
                            double var_p2,double var_ph1,double var_ph2,int iono )
{/*-------------------------------------------------------------------
     Purpose: To smoth code by phase
     -------------------------------------------------------------------
     Input:  TGPS - current GPS time (seconds of week)
             P1 - P1 in meters
             P2 - P2 in meters
             PH1 - Phase L1 in cycles
             PH2 - Phase L2 in cycles
             var_p1 - Variance of P1
             var_p2 - Variance of P2
             var_ph1 - Variance of PH1
             var_ph2 - Variance of PH2
             iono - if iono=2 => iono-free, otherwise => L1 and L2

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
    dtime = TGPS - dtime;
    if(dtime==0)
    {
       Start_Smooth_Code(TGPS,P1,P2, PH1, PH2,var_p1,var_p2,iono);
    }
    else
    {

       double difrng = P2 - P1;
       if(difrng == 0.0)
       {
	   quality = 4;
           PR_Filt_Bef = P1;
           return;
       }
       else{

              Prediction(P1,P2,PH1,PH2,iono);
              Smoothing(P1,P2,PH1,PH2,var_p1,var_p2,var_ph1,var_ph2,iono );
       }
    }

return;
}
//------------------------------------------------------------------------------------
void FilterCode::Start_Smooth_Code(double TGPS,double P1,double P2, double PH1, double PH2,
                                     double var_p1,double var_p2,int iono)
{/*-------------------------------------------------------------------
     Purpose: Start the algorithm of code smothed by phase
     -------------------------------------------------------------------
     Input: TGPS - current GPS time (seconds of week)
             P1 - P1 in meters
             P2 - P2 in meters
             PH1 - Phase L1 in cycles
             PH2 - Phase L2 in cycles
             var_p1 - Variance of P1
             var_p2 - Variance of P2
             var_ph1 - Variance of PH1
             var_ph2 - Variance of PH2
             iono - if iono=2 => iono-free, otherwise => L1 and L2

     Output:
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              FAPESP PROCESS: 05/03522-1
     Date: July of 2010
     Updated in may of 2013
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/
    dtime = TGPS;

    if (iono ==1)
    {   //Start variables to smooth P1 and P2 by carrier
        quality = 10;

        double rngdif = P2 - P1;

        if( rngdif == 0.0 ) quality = 0;

        //smoothed code at initial epoch is equal the own pseudorange
        p1_filt = P1;
        p2_filt = P2;

        //Store observables to use in the next epoch
        ph1_bef = PH1;          //phase L1
        ph2_bef = PH2;          //phase L2
        p1_filt_bef = p1_filt;  //smothed P1
        p2_filt_bef = p2_filt;  //smoothed P2

        var_p1_filt = var_p1;
        var_p2_filt = var_p2;

        p1_pred = 0.0;
        p2_pred =0.0;
    }
    else if(iono==2)
    {
        //Start variables to smooth ionosphere free pseudorange by carrier
        double var_if_code = (pow(m1,2)*var_p1) + (pow(m2_code,2)*var_p2);

        quality = 10;

        double rngdif = P2 - P1;

        if( rngdif == 0.0 ) quality = 0;

        PR_IF = m1*P1 - m2_code*P2;  //ionosphere free pseudorange
        p1_filt = PR_IF;             //smoothed pseudorange

        //Store observables (last epoch) to use in the next epoch
        ph1_bef = PH1;
        ph2_bef = PH2;

        PR_IF_Filt_Bef = PR_IF;  //smoothed pseudorange from last epoch

        var_p1_filt = var_if_code;

        PR_IF_Pred = 0.0;
    }

    cont_epoch=1;   //number of points in the phase data arc

return;
}
//------------------------------------------------------------------------------------
void FilterCode::Prediction(double p1_k,double p2_k,double ph1_k,double ph2_k,int iono)
{/*-------------------------------------------------------------------
     Purpose:To compute prediction of PR smoothed
     -------------------------------------------------------------------
     Input: p1_k - P1 of current epoch
            p2_k - P2 of current epoch
            ph1_k - Ph1 of current epoch
            ph2_k - Ph2 of current epoch
            iono - if iono=2 => iono-free, otherwise => L1 and L2

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


 double D1_kbef_k(0.0),D2_kbef_k(0.0),D_kbef_k(0.0),M1_kbef_k(0.0),M2_kbef_k(0.0),
        curr_IF_Ph(0.0),old_IF_Ph(0.0);


     if (iono == 1)
     {
        //integrated doppler measurement - relative distance
        D1_kbef_k = wlf1*(ph1_k - ph1_bef);    //L1
        D2_kbef_k = wlf2*(ph2_k - ph2_bef);    //L2

        //ion-free relative distance
        D_kbef_k  = m1*D1_kbef_k - m2_code*D2_kbef_k;

       /* //for test.. MHA
        double ionfree_k = m1*ph1_k - m2_phase*ph2_k;
        double ionfree_kbef = m1*ph1_bef -   m2_phase*ph2_bef;
        D_kbef_k =  wlf1*(ionfree_k -  ionfree_kbef);
       */

        //Pseudorange difference projection
        M1_kbef_k = 2.0*D_kbef_k - D1_kbef_k;
        M2_kbef_k = 2.0*D_kbef_k - D2_kbef_k;

        //Prediction for P1 and P2
        p1_pred = p1_filt_bef + M1_kbef_k;
        p2_pred = p2_filt_bef + M2_kbef_k;
      }
      else if(iono==2)
      {  //Despite smoothing for ion-free pseudorange is implemented
        //we are using smothing for P1 and P2 for TEC computation
        curr_IF_Ph = m1*(ph1_k*wlf1) - m2_code*(ph2_k*wlf2);    //Current ion-free phase (meters)
        old_IF_Ph = m1*(ph1_bef*wlf1) - m2_code*(ph2_bef*wlf2); //ion-free phase from last epoch (meters)

        PR_IF_Pred = PR_IF_Filt_Bef + (curr_IF_Ph - old_IF_Ph);  //Predicted Pseudorange ion-free

      }

return;
}
//------------------------------------------------------------------------------------
void FilterCode::Smoothing(double p1,double p2,double ph1,double ph2,double var_p1,double var_p2,
                           double var_ph1,double var_ph2,int iono)
{/*-------------------------------------------------------------------
     Purpose: To compute the smothing of pseudorange
     -------------------------------------------------------------------
     Input: P1 - P1 in meters
            P2 - P2 in meters
            PH1 - Phase L1 in cycles
            PH2 - Phase L2 in cycles
            var_p1 - Variance of P1
            var_p2 - Variance of P2
            var_ph1 - Variance of PH1
            var_ph2 - Variance of PH2
            iono - if iono=2 => iono-free, otherwise => L1 and L2

     Output:
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              FAPESP PROCESS: 05/03522-1
     Date: July of 2010
     Updated in may of 2013
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

 double weight(0.0);

   cont_epoch++; //number of points in the phase data arc

   if (iono == 1)
   {   //P1 smoothed by phase
       weight = (var_p1 + cont_epoch*var_ph1)/(cont_epoch*(var_p1+var_ph1));
       p1_filt = p1_pred + weight*(p1 - p1_pred);

       //P2 smoothed by phase
       weight = (var_p2 + cont_epoch*var_ph2)/(cont_epoch*(var_p2+var_ph2));
       p2_filt = p2_pred + weight*(p2 - p2_pred);

       //Variance of smoothed codes
       var_p1_filt = (var_p1*(var_p1 + cont_epoch*var_ph1))/(cont_epoch*(var_p1+var_ph1));
       var_p2_filt = (var_p2*(var_p2 + cont_epoch*var_ph2))/(cont_epoch*(var_p2+var_ph2));

       //update observables to used in the next epoch
       ph1_bef = ph1;
       ph2_bef = ph2;

       p1_filt_bef = p1_filt;
       p2_filt_bef = p2_filt;

    }
    else if(iono==2)
    {
       var_if_code = (pow(m1,2)*var_p1) + (pow(m2_code,2)*var_p2);
       var_if_phase = (pow(m1,2)*var_ph1) + (pow(m2_code,2)*var_ph2);

       PR_IF = m1*p1 - m2_code*p2;  //current ion-free pseudorange

       weight = (var_if_code + cont_epoch*var_if_phase)/(cont_epoch*(var_if_code+var_if_phase));

       p1_filt = PR_IF_Pred + weight*(PR_IF - PR_IF_Pred); //smoothed pseudorange

       var_p1_filt = ((cont_epoch-1)*var_if_phase + var_if_code)/(cont_epoch);

       //update observables to be used in the next epoch
       ph1_bef = ph1;
       ph2_bef = ph2;

       PR_IF_Filt_Bef = p1_filt;

    }
return;
}
//---------------------------------------------------------------------------
double FilterCode::Get_ph1_bef()
{
    return(ph1_bef);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_ph2_bef()
{
    return(ph2_bef);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p1_filt_bef()
{
    return(p1_filt_bef);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p2_filt_bef()
{
    return(p2_filt_bef);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p1_pred()
{
    return(p1_pred);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p2_pred()
{
    return(p2_pred);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p1_filt()
{
    return(p1_filt);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_p2_filt()
{
    return(p2_filt);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_var_p1_filt()
{
    return(var_p1_filt);
}
//------------------------------------------------------------------------------------
double FilterCode::Get_var_p2_filt()
{
    return(var_p2_filt);
}
//------------------------------------------------------------------------------------
int FilterCode::Get_quality()
{
    return(quality);
}
//------------------------------------------------------------------------------------
int FilterCode::Get_cont_epoch()
{
    return(cont_epoch);
}
//------------------------------------------------------------------------------------
void FilterCode::Set_ph1_bef(double value)
{
    ph1_bef = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_ph2_bef(double value)
{
    ph2_bef = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p1_filt_bef(double value)
{
    p1_filt_bef = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p2_filt_bef(double value)
{
    p2_filt_bef = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p1_pred(double value)
{
    p1_pred = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p2_pred(double value)
{
    p2_pred = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p1_filt(double value)
{
    p1_filt = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_p2_filt(double value)
{
    p2_filt = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_var_p1_filt(double value)
{
    p1_filt = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_var_p2_filt(double value)
{
    p2_filt = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_quality(int value)
{
    quality = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_cont_epoch(int value)
{
    cont_epoch = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_var_if_code(double value)
{
    var_if_code = value;
}
//------------------------------------------------------------------------------------
void FilterCode::Set_var_if_phase(double value)
{
  var_if_phase = value;
}
//------------------------------------------------------------------------------------

}//namespace

