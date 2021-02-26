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
           PD_Filt_Bef = P1;
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
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/
    dtime = TGPS;

    double var_if_code = (pow(m1,2)*var_p1) + (pow(m2_code,2)*var_p2);

    if(iono==2)
    {
      quality = 10;
 	  double rngdif = P2 - P1;
	  if( rngdif == 0.0 )
		quality = 0;

      Phase_Range = 0.0;
      Sum_Iono_Diff = 0.0;
      Sum_Ambiguity = 0.0;
      Filter_Correction = 0.0;

      PD_IF = m1*P1 - m2_code*P2; //ion-free code

      PD_IF_Pred = PD_IF;      //Predicted smoothed Pseudorange
      PD_IF_Filt_Bef = PD_IF;  //smoothed pseudorange from last epoch
      p1_filt = PD_IF;         //smoothed pseudorange

      P1_Filt_Ini = PD_IF;     //Pseudorange at initial epoch

      //Store observables (last epoch) to use in the next epoch
      ph1_bef = PH1;
      ph2_bef = PH2;
      p1_bef = P1;
      p2_bef = P2;

      var_p1_filt = var_if_code;

    }
    else
    {
        p1_filt_bef = P1;
        p2_filt_bef = P2;

        p1_filt = P1;
        p2_filt = P2;

        ph1_bef = PH1;
        ph2_bef = PH2;

        var_p1_filt = var_p1;
        var_p2_filt = var_p2;

        p1_pred = 0.0;
        p2_pred =0.0;
    }

    cont_epoch=1;

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


 double D1_kbef_k,D2_kbef_k,D_kbef_k,M1_kbef_k,M2_kbef_k;

     if(iono==2)
     {  double curr_IF_Ph =  m1*(ph1_k*wlf1) - m2_code*(ph2_k*wlf2);    //Current ion-free phase (meters)

        double old_IF_Ph =  m1*(ph1_bef*wlf1) - m2_code*(ph2_bef*wlf2); //phase ion-free phase from last epoch

        diff_phase = curr_IF_Ph - old_IF_Ph;

        double curr_IF_PD = m1*p1_k - m2_code*p2_k;   //current ion-free for code
        double old_IF_PD = m1*p1_bef - m2_code*p2_bef; //code ion-free from last epoch

        double diff_PD = curr_IF_PD - old_IF_PD;

        Delta_Code_Phase =  diff_PD - diff_phase;  //Difference between code and phase (ambiguity)

        //------------------------------------------------------------------------

        PD_IF_Pred = PD_IF_Filt_Bef + (curr_IF_Ph - old_IF_Ph);  //Predicted Pseudorange

     }
     else
     {
        //integrated doppler measurement - relative distance
        D1_kbef_k = wlf1*(ph1_k - ph1_bef);    //L1
        D2_kbef_k = wlf2*(ph2_k - ph2_bef);    //L2

        //ion-free relative distance
        D_kbef_k  = m1*D1_kbef_k - m2_code*D2_kbef_k;

        //Pseudorange difference projection
        M1_kbef_k = 2.0*D_kbef_k - D1_kbef_k;
        M2_kbef_k = 2.0*D_kbef_k - D2_kbef_k;

        //Prediction for P1 and P2
        p1_pred = p1_filt_bef + M1_kbef_k;
        p2_pred = p2_filt_bef + M2_kbef_k;
      }

return;
}
//------------------------------------------------------------------------------------
void FilterCode::Smoothing(double p1,double p2,double ph1,double ph2,double var_ca,double var_p2,
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
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

 double weight;
 double var_if_code = (pow(m1,2)*var_ca) + (pow(m2_code,2)*var_p2);
 double var_if_phase = (pow(m1,2)*var_ph1) + (pow(m2_code,2)*var_ph2);

 const double Max_Ambiguity = 3.0;

   cont_epoch++;

   if(iono==2)
   {

       Sum_Iono_Diff += diff_phase; //Sum of ion-free carrier/code from initial epoch

       Sum_Ambiguity += Delta_Code_Phase; //Sum of ambiguities
       Last_correction = (Sum_Ambiguity - Filter_Correction)/cont_epoch; //Last correction to ambiguity

       Filter_Correction += Last_correction; //Filtered correction

       if( fabs(Filter_Correction) > Max_Ambiguity )//Quality control
       {
		   quality = 0;
		   return;
       }

       //ion-free phase taking into account ambiguity (diff code-phase) since from last epoch
       Phase_Range = P1_Filt_Ini + Sum_Iono_Diff;

      /* p1_filt = P1_Filt_Ini + Filter_Correction + Sum_Iono_Diff;

       var_p1_filt = fabs(Last_correction);

       ph1_bef = ph1;
       ph2_bef = ph2;
       p1_bef = p1;
       p2_bef = p2;
       */

       //------------------------------------------------------------------------
       PD_IF = m1*p1 - m2_code*p2;  //ion-free code

       weight = (var_if_code + cont_epoch*var_if_phase)/(cont_epoch*(var_if_code+var_if_phase));

       p1_filt = PD_IF_Pred + weight*(PD_IF - PD_IF_Pred); //smoothed pseudorange

       PD_IF_Filt_Bef = p1_filt;

       var_p1_filt = ((cont_epoch-1)*var_if_phase + var_if_code)/(cont_epoch);

       //update observables to next epoch
       ph1_bef = ph1;
       ph2_bef = ph2;
       p1_bef = p1;
       p2_bef = p2;

   }
   else
   {   //P1 smoothed by phase
       weight = (var_ca + cont_epoch*var_ph1)/(cont_epoch*(var_ca+var_ph1));
       p1_filt = p1_pred + weight*(p1 - p1_pred);

       //P2 smoothed by phase
       weight = (var_p2 + cont_epoch*var_ph2)/(cont_epoch*(var_p2+var_ph2));
       p2_filt = p2_pred + weight*(p2 - p2_pred);

       //Variance of smoothed codes
       var_p1_filt = (var_ca*(var_ca + cont_epoch*var_ph1))/(cont_epoch*(var_ca+var_ph1));
       var_p2_filt = (var_p2*(var_p2 + cont_epoch*var_ph2))/(cont_epoch*(var_p2+var_ph2));
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

}//namespace

