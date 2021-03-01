/*-------------------------------------------------------------------
    Purpose: class to Read GIM interpolate TEC or compute TEC from Pseudorange
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
#include "Class_TEC.h" // class's header file


namespace CTEC{

Class_TEC::Class_TEC()
{   // class constructor
    //Default values
    Re = 6371000.0;
    hion = 450000.0;

    Num_Value = 0;  //Number of data read in the GIM file

    for(int i =0;i<NSAT;i++)
    {  p1p2[i][0] = 0.0;
       p1c1[i][0] = 0.0;

       p1p2[i][1] = 0.0;
       p1c1[i][1] = 0.0;

       s_tec[i] = 0.0;
    }//for

    br=0.0; //receiver hardware delay - DCB(P1-P2) (m)
    sbr=0.0;  //standard deviation of the reciever DCB(P1-P2)

}
//------------------------------------------------------------------------------------
Class_TEC::Class_TEC( double sigca = 0.6,double sigp2=0.8,
                      double sigph1=0.006,double sigph2=0.008,
                      double dcbr=0,double sdcbr=0)
{/*-------------------------------------------------------------------
     Purpose: Class constructor
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

    //Default values
    Re = 6371000.0;
    hion = 450000.0;

    Num_Value = 0;  //Number of data read in the GIM file

    var_ca = pow(sigca,2);           //variance of CA - meter
    var_p2 = pow(sigp2,2);           //variance of P2 - meter
    var_ph1 = pow(sigph1,2);         //variance of Phase L1 - meter
    var_ph2 = pow(sigph2,2);         //variance of Phase L2 - meter

    br = dcbr;                       //receiver hardware delay
    sbr = sdcbr;                     //standard deviation of receiver hardware delay

    for(int i =0;i<NSAT;i++)
    {  p1p2[i][0] = 0.0;
       p1c1[i][0] = 0.0;

       p1p2[i][1] = 0.0;
       p1c1[i][1] = 0.0;

       s_tec[i] = 0.0;
    }//for
}
//------------------------------------------------------------------------------------
Class_TEC::~Class_TEC()
{  // class destructor
	// insert your code here
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_ZL()
{ //zenith angle at the pierce point
  //unit: rad
  return(zl);
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_LatIon()
{
     return(lat_ion);      //Latitude at pierce point
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_LambIon()
{
     return(lamb_ion);     //Longitude at pierce point
}
//------------------------------------------------------------------------------------
int Class_TEC::Get_Gim_Flag()
{
     return(Gim_Flag);
}
//------------------------------------------------------------------------------------
double Class_TEC::GetVTEC()
{
 return(VTEC);  //Vertical TEC
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Epoch_First_MAP()
{
    return (Epoch_First_MAP); //in MJD
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Epoch_Last_Map()
{
    return (Epoch_Last_Map); //in MJD
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Interval()
{
    return (Interval);  //Seconds
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Ini_Lat()
{
    return (Ini_Lat); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_End_Lat()
{
    return (End_Lat); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Lat_Increment()
{
    return (Lat_Increment); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Ini_Lon()
{
    return (Ini_Lon); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_End_Lon()
{
    return (End_Lon); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Lon_Increment()
{
    return (Lon_Increment); //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Elevation_Cuttof()
{
    return (Elevation_Cuttof);  //in Degrees
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Base_Radius()
{
    return (Base_Radius);   //in Km
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Ini_Height()
{
    return (Ini_Height);   //in Km
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_End_Height()
{
    return (End_Height);   //in Km
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_Height_Incremet()
{
    return (Height_Incremet);  //in km
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_P1P2(int prn)
{/*-------------------------------------------------------------------
     Purpose: Return satellite P1-P2 hardware delay (ns)
     -------------------------------------------------------------------
     Input: prn - Satellite number

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

    return(p1p2[prn][0]);
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_P1C1(int prn)
{/*-------------------------------------------------------------------
     Purpose: Return satellite P1-C1 hardware delay (ns)
     -------------------------------------------------------------------
     Input: prn - Satellite number

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

    return(p1c1[prn][0]);
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_S_TEC(int prn)
{/*-------------------------------------------------------------------
     Purpose: Return standard deviation of STEC
     -------------------------------------------------------------------
     Input: prn - Satellite PRN

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

    return (s_tec[prn]);
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_br()                        //receiver hardware delay
{
    return(br);
}
//------------------------------------------------------------------------------------
double Class_TEC::Get_sbr()
{
    return(sbr);
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_var_ca(double value)
{
    var_ca = value;
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_var_p2(double value)          //variances for CA and P2
{
    var_p2 = value;
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_var_ph1(double value)
{
    var_ph1 = value;
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_var_ph2(double value)         //variances for phase L1 and phase L2
{
    var_ph2 = value;
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_br(double value)              //receiver hardware delay
{
    br = value;
}
//------------------------------------------------------------------------------------
void Class_TEC::Set_sbr(double value)
{
    sbr = value;
}
//------------------------------------------------------------------------------------
double Class_TEC::Mean_Phase_Code(double ph1, double ph2, double p1, double p2,
                       double &Mean_Diff_Li_Pi,double &Last_Diff_Li_Pi,int &Count_Mean,
                       int &arc_number,int initial_epoch,bool Cycle_Slip_Flag,int prn )
{  double Li(0.0),Pi(0.0),Curr_Diff_Li_Pi(0.0),Curr_Last_Diff(0.0);
   double f1_quad = pow(freq1,2),
       f2_quad = pow(freq2,2),
       alfa = 1.0/( 40.309*(1.0/f2_quad-1.0/f1_quad) ),
       //values in meters
       Bsm = c*(p1p2[prn][0]*1.0e-9),   //satellite p1-p2 DCB
       Brm = c*(br*1.0e-9);             //receiver DCB

   /*   di2d1=f1_quad/(f1_quad-f2_quad),
      di2d2=f2_quad/(f1_quad-f2_quad),

      d1 = di2d1*(Brm+Bsm),
      d2 = di2d2*(Brm+Bsm),
      di = d2-d1;
   */

    Li = wlf1*ph1 - wlf2*ph2;         // Geomery free phase (meters)
    Pi = (p2 - p1);                   // Geometry free pseudorange (meters)

    Curr_Diff_Li_Pi = Li - Pi - (Brm+Bsm);   //time-average phase - code - dcb

    if(Cycle_Slip_Flag)
    { //it means that mean filter will be re-initialized in the next epoch
      Count_Mean=0;
      Mean_Diff_Li_Pi = Curr_Diff_Li_Pi = Last_Diff_Li_Pi = 0.0;
    }
    else if (initial_epoch==0 || Count_Mean==0) //first epoch or starting new computation
    {
      Mean_Diff_Li_Pi = Curr_Diff_Li_Pi;
      Last_Diff_Li_Pi = Mean_Diff_Li_Pi; //Diff Phase and code to be used in the next epoch
      Count_Mean=1;
      arc_number++;

    }//if
    else
    {
        Count_Mean++;  //current number of points in the data arc

        Curr_Last_Diff = Curr_Diff_Li_Pi - Last_Diff_Li_Pi;  //difference form current and last epoch

        Mean_Diff_Li_Pi = Last_Diff_Li_Pi + (Curr_Last_Diff/Count_Mean); //Mean filter <Li - Pi>

        Last_Diff_Li_Pi = Mean_Diff_Li_Pi; //Diff Phase and code to be used in the next epoch
   }//else

return (Mean_Diff_Li_Pi);
}
//------------------------------------------------------------------------------------
long double Class_TEC::Calc_STEC_PH_PR(double ph1, double ph2,double p1,double p2,double Mean_Diff_Li_Pi,int prn )
{
  double Li(0.0),Pi(0.0),Curr_Diff_Li_Pi(0.0),STEC(0.0);

double f1_quad = pow(freq1,2),
       f2_quad = pow(freq2,2),
       alfa = 1.0/( 40.309*(1.0/f2_quad-1.0/f1_quad) );

/*       //values in meters
       Bsm = c*(p1p2[prn][0]*1.0e-9),   //satellite p1-p2 DCB
       Brm = c*(br*1.0e-9),             //receiver DCB

      di2d1=f2_quad/(f1_quad-f2_quad),
      di2d2=f1_quad/(f1_quad-f2_quad),

      d1 = di2d1*(Brm+Bsm),
      d2 = di2d2*(Brm+Bsm),
      di = d2-d1;
  */
     Li = wlf1*ph1 - wlf2*ph2;    // Geomery free phase (meters)

     //slant TEC in the first epoch from =pseudorange
    // STEC = alfa*(  Li - Mean_Diff_Li_Pi + Brm + Bsm  );
     STEC = alfa*(  Li - Mean_Diff_Li_Pi  );
    // double teste = Calc_STEC_PR(p1,p2,prn);

 return (STEC);      //return STEC value
}
//------------------------------------------------------------------------------------
long double Class_TEC::Calc_STEC_PR(double P1,double P2,int prn)
{/*-------------------------------------------------------------------
    Purpose: Compute slant TEC STEC from raw pseudorange CA and P1
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
    Observation: DCB P1-C1 applied to CA to compatibilize CA with P1
    -------------------------------------------------------------------*/

double c = 299792458.0;

double f1_quad = pow(freq1,2),
       f2_quad = pow(freq2,2),
       alfa = (f1_quad*f2_quad)/(40.309*(f2_quad-f1_quad)),
       //values in meters
       Bsm = c*(p1p2[prn][0]*1.0e-9),   //satellite p1-p2 DCB
       Brm = c*(br*1.0e-9);             //receiver DCB

double var_bsm = pow(p1p2[prn][1],2), //P1-P2 rms
       var_brm = pow(sbr,2);          //DCBr rms

long double Stec=0.0,var_tec=0.0;

    //slant TEC from pseudorange
    Stec = alfa*(  (P1 - P2) - Brm - Bsm  );

    //Computing Standard Deviation of TEC
    var_tec = pow(alfa,2)*(var_ca + var_p2 + var_brm + var_bsm  ); //(ele/m^2)
    s_tec[prn] = sqrt(var_tec);


return (Stec);    //return STEC value

}
//------------------------------------------------------------------------------------
void Class_TEC::Read_DCB_Sat(char arqu_p1c1[],char arqu_p1p2[])
{/*-------------------------------------------------------------------
    Purpose: Read DCBs (P1-P2) and (P1-C1) from CODE files
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
 char temp[500]={"\0"};
 string linha;
 int prn;
 int Start_Col_DCB,Start_Col_RMS;

    ifstream arqu1(arqu_p1c1);
    ifstream arqu2(arqu_p1p2);

    if(arqu1.fail())
    {
        cout<<"It was not possible to open DCB File: "<< arqu_p1c1<< " in Read_DCB_Sat subroutine."<<endl;
        exit(-1);
    }//if

    if(arqu2.fail() )
    {
        cout<<"It was not possible to open DCB File: "<< arqu_p1p2<< " in Read_DCB_Sat subroutine."<<endl;
        exit(-1);
    }//if

    while(!arqu1.eof())
    {
      //get line from P1-C1 file
      arqu1.getline(temp,500);
      linha = temp;
      int j = linha.substr(0,3).compare("PRN");
      if(j==0)
      { //Find the collun where are the DCBs!
        Start_Col_DCB = linha.find("VALUE (NS)");
        Start_Col_RMS = linha.find("RMS (NS)");
      }//if

      int i = linha.substr(0,1).compare("G");
      if(i==0)
      {
        prn = atoi(linha.substr(1,2).c_str());

        if(prn>=1&&prn<=32){ //only GPS
           p1c1[prn][0]=atof(linha.substr(Start_Col_DCB,10).c_str());   //Read P1-C1
           p1c1[prn][1]=atof(linha.substr(Start_Col_RMS,10).c_str());   //Read P1-C1 RMS
        }//if
      }
    }//while

    while(!arqu2.eof())
    {
      //get line from P1-P2 file
      arqu2.getline(temp,500);
      linha = temp;

      int j = linha.substr(0,3).compare("PRN");
      if(j==0)
      { //Find the collun where are the DCBs!
        Start_Col_DCB = linha.find("VALUE (NS)");
        Start_Col_RMS = linha.find("RMS (NS)");
      }//if

      int i = linha.substr(0,1).compare("G");
      if(i==0)
      {
        prn = atoi(linha.substr(1,2).c_str());
        if(prn>=1&&prn<=32){ //only GPS
           p1p2[prn][0] = atof(linha.substr(Start_Col_DCB,10).c_str());  //Read P1-P2
           p1p2[prn][1] = atof(linha.substr(Start_Col_RMS,10).c_str()); //Read P1-P2 RMS
        }//if
      }


    }

}
//------------------------------------------------------------------------------------
double Class_TEC::DJUL(int Year,int Month, double t)
{/*-------------------------------------------------------------------
     Purpose:  To Compute Julian Date
     -------------------------------------------------------------------
     Input:  Year
             Month
             t = DAY+HOUR/24.0+MIN/1440.0+SEC/86400.0

     Output: Julian date
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP

     Date: January of 2010
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

	int i, j, k, m;
	double djul;

	j = Year;
	m = Month;
	if( m <= 2 )
	{
	   j -= 1;
	   m += 12;
    }

	i = j/100;
	k = 2 - i + i/4;
	djul = (365.25*j - fmod(365.25*j,1.0)) - 679006.0;
	djul += int(30.6001*(m+1)) + t + k;

	return( djul );
}
//------------------------------------------------------------------------------------
void  Class_TEC::ReadMap(int numval,int iexp, double tecval[],int &irc,ifstream &f1)
{/*-------------------------------------------------------------------
     Purpose: Read Ionex Data
     -------------------------------------------------------------------
     Input:

     Output:
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP

     Date: January of 2010
     -------------------------------------------------------------------
     Observation: The C++ subroutine ReadMap was developed based on the
                  original fortran subroutine RDIXDT from author (S.SCHAER, 1997)
     -------------------------------------------------------------------*/
char LINE[82]={"\0"};
string line1;

int auxval[16]={0}, iaux, iaux_, ilin, ilin_, ival, ival1, ival2,
    lfnerr, naux, nlin;

    try{
	       // COMPUTE NUMBER OF DATA LINES TO BE READ
	       nlin = (numval - 1)/16 + 1;

	       // READ SINGLE DATA BLOCK
	       for( ilin=1; ilin <= nlin; ilin++ )
           {
	       	   ival1 = 16*ilin - 15;
	       	   ival2 = 16*ilin;

		       if( ival2 > numval )
                  ival2 = numval;

               for(int i=0;i<82;i++) LINE[i]='\0';

               line1="";

               f1.getline(LINE,82);
               line1=LINE;

		       naux = ival2 - ival1 + 1;

               int i=0;

               for (int i=0;i<16;i++) auxval[i]=0.0;


		       for( iaux=1; iaux <= naux; iaux++ )
               {
                   sscanf(line1.substr(i,i+5) .c_str(),"%d",&auxval[iaux-1]);
                   i+=5;
		       }//for

		       // CONVERT TEC/RMS VALUES FROM ACTIVE UNIT INTO TECU
		       for( iaux=1; iaux <= naux; iaux++ )
               {
			        iaux_ = iaux - 1;
			        ival = ival1 + iaux - 1;

			        if( auxval[iaux_] != 9999 )
             		  tecval[ival-1] = (pow(10.0,iexp)*auxval[iaux_]);
			        else
				      tecval[ival-1] = 999.9;
			   }//for

           }//for
     }catch(...){
        irc = 0;
        return;
     }

	 irc = 1;
	 return;

}
//------------------------------------------------------------------------------------
void Class_TEC::ReadHeaderGim(ifstream &f1, int &Iexp, int &Header_Flag)

{/*-------------------------------------------------------------------
     Purpose: Read Header of GIM file
     -------------------------------------------------------------------
     Input: f1 - ifstream file

     Output: Iexp - Exponent
             Header_Flag = 1 -> ok
             Header_Flag = 0 -> Header Problem
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP

     Date: January of 2010
     -------------------------------------------------------------------
     Observation: The C++ subroutine ReadHeaderGim was developed based on the
                  original fortran subroutine RDIXHD from author (S.SCHAER, 1997)
     -------------------------------------------------------------------*/

char Line[83]="\0", AuxDataBlock0[60]="\0";
double GridFactor[4],G_Version0 = 1,G_Version = 0;

string line1,AuxDataBlock="\0", SatSystem="\0", Program="\0", Agency="\0", Date="\0",
       MappFunction="\0", ObsUsed="\0";

const int NumReg = 22;
int RegFlag[NumReg]={0},iint,NMAP,idim;
int itemp;

    try{
            Iexp=-1;

            Header_Flag=0;

            f1.seekg(0);

            //Read and check the first ionex record
            f1.getline(Line,82);
            line1=Line;
            if(line1.substr(60,20).compare("IONEX VERSION / TYPE") != 0)
            {
                 printf("Invalid First Ionex record");
                 Header_Flag=0;
                 return;
            }//if

            if(line1.substr(20,20).compare("IONOSPHERE MAPS     ") != 0)
            {
                printf("Invalid ionex Label");
                Header_Flag=0;
                return;
            }

            //SATELLITE SYSTEM
            SatSystem = line1.substr(40,20);

            //IONEX VERSION NUMBER
            sscanf(line1.c_str(),"%lf",&G_Version);
            if(G_Version != G_Version0)
            {
                printf("Invalid Ionex Version. Greater than 1.0");
                Header_Flag=0;
                return;
            }

            RegFlag[0] = 1;

            do
            {
                  for(int i=0;i<83;i++) Line[i]='\0';

                  f1.getline(Line,83);

                  //cout<<Line<<endl;

                  line1=Line;

                  //cout<<line1.substr(60,21)<<endl;

                  // LOOK FOR CURRENT KEYWORD AND STORE CORRESPONDING INFORMATION
                  if(line1.substr(60,20).compare("PGM / RUN BY / DATE ") == 0) // PROGRAM, AGENCY, DATE
                  {
                      Program = line1.substr(0,20);
                      Agency = line1.substr(20,20);
                      Date = line1.substr(40,20);
                      RegFlag[1]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("COMMENT             ") == 0)  //SKIP COMMENT LINES
                  {
                      RegFlag[2]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("DESCRIPTION         ") == 0)  //SKIP DESCRIPTION LINES
                  {
                      RegFlag[3]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("EPOCH OF FIRST MAP  ") == 0)  //EPOCH OF FIRST MAP
                  {
                      int IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC;
                      sscanf(line1.c_str(),"%6d %6d %6d %6d %6d %6d",&IYEAR,&IMONTH,&IDAY,&IHOUR,&IMIN,&ISEC);
                      RegFlag[4]=1;

                      Epoch_First_MAP=DJUL(IYEAR,IMONTH,(IDAY+IHOUR/24.0+IMIN/1440.0+ISEC/86400.0) );

                  }//else if
                  else if(line1.substr(60,20).compare("EPOCH OF LAST MAP   ") == 0) //EPOCH OF LAST MAP
                  {
                      int IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC;
                      sscanf(line1.c_str(),"%6d %6d %6d %6d %6d %6d",&IYEAR,&IMONTH,&IDAY,&IHOUR,&IMIN,&ISEC);
                      RegFlag[5]=1;

                      Epoch_Last_Map=DJUL(IYEAR,IMONTH,(IDAY+IHOUR/24.0+IMIN/1440.0+ISEC/86400.0) );
                  }//else if
                  else if(line1.substr(60,20).compare("INTERVAL            ") == 0)
                  {
                      sscanf(line1.c_str(),"%6d",&iint);
                      RegFlag[6]=1;
                      Interval=double(iint);
                  }//else if
                  else if(line1.substr(60,20).compare("# OF MAPS IN FILE   ") == 0)
                  {
                      sscanf(line1.c_str(),"%6d",&NMAP);
                      RegFlag[7]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("MAPPING FUNCTION    ") == 0)
                  {
                      MappFunction = line1.substr(2,4);
                     //sscanf(line1.c_str(),"%s",MappFunction.c_str());
                      RegFlag[8]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("ELEVATION CUTOFF    ") == 0)
                  {
                      sscanf(line1.c_str(),"%lf",&Elevation_Cuttof);
                      RegFlag[9]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("OBSERVABLES USED    ") == 0)
                  {
                      ObsUsed = line1.substr(0,60);
                      RegFlag[10]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("# OF STATIONS       ") == 0)
                  {
                      sscanf(line1.c_str(),"%6d",&itemp);
                      RegFlag[11]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("# OF SATELLITES     ") == 0)
                  {
                      sscanf(line1.c_str(),"%6d",&itemp);
                      RegFlag[12]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("BASE RADIUS         ") == 0)
                  {
                      sscanf(line1.c_str(),"%lf",&Base_Radius);
                      RegFlag[13]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("MAP DIMENSION       ") == 0)
                  {
                      sscanf(line1.c_str(),"%6d",&idim);
                      RegFlag[14]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("HGT1 / HGT2 / DHGT  ") == 0)
                  {
                      sscanf(line1.c_str(),"%lf %lf %lf",&Ini_Height,&End_Height,&Height_Incremet);
                      RegFlag[15]=1;
                 }//else if
                  else if(line1.substr(60,20).compare("LAT1 / LAT2 / DLAT  ") == 0)
                  {
                      sscanf(line1.c_str(),"%lf %lf %lf",&Ini_Lat,&End_Lat,&Lat_Increment);
                      RegFlag[16]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("LON1 / LON2 / DLON  ") == 0)
                  {
                      sscanf(line1.c_str(),"%lf %lf %lf",&Ini_Lon,&End_Lon,&Lon_Increment);
                      RegFlag[17]=1;
                  }//else if
                  else if(line1.substr(60,20).compare("EXPONENT            ") == 0) // EXPONENT
                  {
                      sscanf(line1.c_str(),"%d",&Iexp);
                      RegFlag[18]=1;
                  } //else if
                  else if(line1.substr(60,20).compare("START OF AUX DATA   ") == 0)  // AUXILIARY DATA BLOCK FOUND
                  {
                      line1.copy(AuxDataBlock0,60,0);
                      RegFlag[19]=1;

                      int NADT=0;
                      do{
                            //read next AUX DATA line
                            for(int i=0;i<81;i++) Line[i]='\0';
                            f1.getline(Line,82);
                            line1=Line;

                            if (strcmp(AuxDataBlock.c_str(),AuxDataBlock0)==0)
                            {
                                 NADT=NADT+1;
                                 if (NADT > 16)
                                 {
                                     printf("So many Aux Data Lines");
                                     Header_Flag=0;
                                     return;
                                 }//if
                            }//if
                      }while (line1.substr(60,20).compare("END OF AUX DATA     ") != 0);

                      RegFlag[20] = 1;
                  }//else


            }while( line1.substr(60,20).compare("END OF HEADER       ") != 0 &&
                    line1.substr(60,20).compare("END OF FILE         ") != 0);


            if(line1.substr(60,20).compare("END OF FILE         ") ==0)
            {
	            printf("End of File reached");
		        Header_Flag = 0;
                return;
            }

            RegFlag[21] = 1;

            // Check ionex Header
	        for( int i=0; i < NumReg; i++ )
            {
		        if( RegFlag[i] == 0 )
                {
		            printf(" Record Missing in the GIM file header");
			        Header_Flag = 0;
                    return;
                }//if
            }//for


	        if( iint == 0 )
            {
	  	        printf(" Interval 0 Not Supported" );
	            Header_Flag = 0;
      	        return;
      		}//if
	        else
            {
	        	if( NMAP > 1 )
                {
	        	    int iint0 = int((Epoch_Last_Map-Epoch_First_MAP)/(NMAP-1)*86400.0);

	                if( iint != iint0 )
      	    	    {
	  		            printf("Incosistent Time Arguments");
	        	        Header_Flag = 0;
            	        return;
	                }//if
                }//if
            }//else

            // SATELLITE SYSTEM
            if( SatSystem.compare(0,3,"BEN") != 0 &&
                SatSystem.compare(0,3,"ENV") != 0 &&
                SatSystem.compare(0,3,"ERS") != 0 &&
                SatSystem.compare(0,3,"GEO") != 0 &&
                SatSystem.compare(0,3,"GLO") != 0 &&
                SatSystem.compare(0,3,"GNS") != 0 &&
                SatSystem.compare(0,3,"GPS") != 0 &&
                SatSystem.compare(0,3,"IRI") != 0 &&
                SatSystem.compare(0,3,"MIX") != 0 &&
	            SatSystem.compare(0,3,"NNS") != 0 &&
                SatSystem.compare(0,3,"TOP") != 0 )
            {
   	            printf("Satellite System Unknown");
                Header_Flag = 0;
                return;
	        }//if


            // MAPPING FUNCTION
	        if( MappFunction.compare("NONE") != 0 &&
                MappFunction.compare("COSZ") != 0 &&
                MappFunction.compare("QFAC") != 0 )
            {
         	      printf("Mapping Function Unknown");
                Header_Flag = 0;
                return;
       	    }//if


      	    if( idim == 2 ) // MAP DIMENSION
            {
      	   	    if( (Ini_Height != End_Height) || (Height_Incremet != 0.0) )
                {
         	         printf("Ionosphere Height Unknown");
                   Header_Flag = 0;
                   return;
                }//if
            }//if
      	    else
            {
      		    if( (Ini_Height == End_Height) || (End_Height == 0.0) )
                {
         	         printf(" Error related to Map Dimension");
                     Header_Flag = 0;
                     return;
      		    }//if
            }//else

            // Check Latitude and Longitude
	        GridFactor[0] = Ini_Lat/Lat_Increment;
	        GridFactor[1] = End_Lat/Lat_Increment;
	        GridFactor[2] = Ini_Lon/Lon_Increment;
	        GridFactor[3] = End_Lon/Lon_Increment;

	        for( int i=0; i < 4; i++ )
            {
	             double xdif = fabs(GridFactor[i]-int(GridFactor[i]));
		         if( xdif > 1.0E-4 )
                 {	Header_Flag = 0;
                    return;
                 }//if

            }//for
    }catch(...){
        Header_Flag = 0;
        return;
    }

     Header_Flag = 1;

     return;

}
//---------------------------------------------------------------------------
int Class_TEC::ReadGIM(char *NameGIM)
{/*-------------------------------------------------------------------
     Purpose: Read Ionex File
     -------------------------------------------------------------------
     Input:

     Output:
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP

     Date: January of 2010
     -------------------------------------------------------------------
     Observation: The C++ subroutine ReadGIM was developed based on the
                  original fortran subroutine RDIXFL from author (S.SCHAER, 1997)
     -------------------------------------------------------------------*/
// 02/05/2013:MHA - Small changes in the source code


const int MaxMap = 67379;
int NHGT,IMAP,IVAL,iaux=0;
int NUMMAP[3]={0};
double XEPO=0.0;

char Line[83]="\0";
string line1;


     ifstream f1;
     int Iexp,Header_Flag=0;
     Gim_Flag=0;


     f1.open(NameGIM); //OPEN IONEX INPUT FILE

     if(f1.fail())
     {
        printf("Error Opening GIM FILE");
        return Gim_Flag;
     }

     //Read GIM Header
     ReadHeaderGim( f1,Iexp,Gim_Flag);

     Re = Base_Radius*1.0E3;  //Equatorial Earth axis in m
     hion = ((End_Height+Ini_Height)/2.0)*1.0E3; //Ionospheric layer height in m

     if(Gim_Flag==0)
       return Gim_Flag;

     double XINT=Interval/86400.0;
     int NMAP=int ((Epoch_Last_Map-Epoch_First_MAP)/XINT)+1;

     int NLAT=int((End_Lat-Ini_Lat)/Lat_Increment)+1;

     int NLON=int((End_Lon-Ini_Lon)/Lon_Increment)+1;

     if(Height_Incremet==0.0)
       NHGT=1;
     else
       NHGT=int((End_Height-Ini_Height)/Height_Incremet)+1;


     Num_Value = NMAP*NHGT*NLAT;  //Number of data read

     if (Num_Value > MaxMap)
     {
        printf ("TOO MANY TEC/RMS VALUES");
        Gim_Flag=0;
         return Gim_Flag;
     }//if

     try { //Allocating memory
            TECMAP = new double*[Num_Value];       // STEP 1: SET UP THE ROWS.

            for (int j = 0; j < Num_Value; j++)
              TECMAP[j] = new double[NLON+2];      // STEP 2: SET UP THE COLUMNS
         }
         catch (std::bad_alloc)
         {  // ENTER THIS BLOCK ONLY IF bad_alloc IS THROWN.
            // YOU COULD REQUEST OTHER ACTIONS BEFORE TERMINATING
            printf("Could not allocate memory for TECMAP...");
            return(0);
         }

     // INITIALIZE TEC/RMS MAPS
     // -----------------------

     for (int i = 0; i < Num_Value; i++)
       for (int j = 0; j < NLON+2; j++)
          TECMAP[i][j] = -999.9;            // ARBITRARY INITIALIZATION

     int ITYP=-1;

     do     //LOOP OVER ALL DATA LINES
     {
          for(int i=0;i<82;i++) Line[i]='\0';
          f1.getline(Line,82);
          line1=Line;

          // LOOK FOR CURRENT KEYWORD AND STORE CORRESPONDING INFORMATION
          // ------------------------------------------------------------
          if(line1.substr(60,20).compare("START OF TEC MAP    ") == 0) //START OF TEC MAP
          {
              sscanf(line1.c_str(),"%6d",&IMAP);
              ITYP=0;
              NUMMAP[ITYP]++;
          }//if
          else if(line1.substr(60,20).compare("START OF RMS MAP    ") == 0)
          {
              Header_Flag = 1;
              break;
          }//else if
          else if(line1.substr(60,20).compare("START OF HEIGHT MAP ") == 0) //START OF HEIGHT MAP
          {
              if (ITYP==-1)
                sscanf(line1.c_str(),"%6d",&IMAP);

              ITYP=2;
              NUMMAP[ITYP]++;

              do{     //SKIP HEIGHT MAP
                   for(int i=0;i<82;i++) Line[i]=0;
                   f1.getline(Line,82);
                   line1=Line;
                }while(line1.substr(60,20).compare("END OF HEIGHT MAP   ") != 0);

          }//else
          else if (line1.substr(60,20).compare("EPOCH OF CURRENT MAP") == 0)
          {
              int IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC;

              //EPOCH OF CURRENT MAP
              sscanf(line1.c_str(),"%6d %6d %6d %6d %6d %6d",&IYEAR,&IMONTH,&IDAY,&IHOUR,&IMIN,&ISEC);

              XEPO = DJUL(IYEAR,IMONTH,IDAY+IHOUR/24.0+IMIN/1440.0+ISEC/86400.0);

              double XEPO0 = Epoch_First_MAP+(IMAP-1)*Interval/86400.0;
              int IDIF = int(fabs((XEPO-XEPO0))*86400.0);

              if (IDIF !=0.0)
              {
                  string error = "EPOCH OF CURRENT MAP WRONG";
                  error.append("MAP NUMBER: ");
                  error += IMAP;

                  cout<< error.c_str()<<endl;

                  Gim_Flag=0;
                  return Gim_Flag;
              }//if
          }//else
          else if (line1.substr(60,20).compare("EXPONENT            ") == 0)
          {
              //EXPONENT
               sscanf(line1.c_str(),"%6d",&Iexp);
          }
          else if (line1.substr(60,20).compare("COMMENT             ") == 0)
          {
              //COMMENT
          }
          else if (line1.substr(60,20).compare("                    ") == 0)
          {
              //SKIP BLANK LINES
          }
          else if (line1.substr(60,20).compare("LAT/LON1/LON2/DLON/H") == 0)
          {
              if (ITYP == -1)
              {
                  string error = " MAP TYPE UNDEFINED";
                  error.append(Line);

                  cout<< error.c_str()<<endl;

                  Gim_Flag=0;
                  return Gim_Flag;
              }//if

              double XLAT,XLON1,XLON2,XDLON,XHGT;

              sscanf(line1.substr(2,6).c_str(),"%lf",&XLAT);
              sscanf(line1.substr(8,6).c_str(),"%lf",&XLON1);

              ReadMap(NLON,Iexp,&TECMAP[iaux][0],Header_Flag,f1);

              TECMAP[iaux][NLON] = XLAT;    //Lat
              TECMAP[iaux][NLON+1] = XEPO;  //Date
              iaux++;

              if(Header_Flag==0)
              {
                 printf("DATA BLOCK NOT CLOSED");
                // return Header_Flag;
              }//if

          }//else
          else
          {
             //IONEX DATA RECORD";
          }


     }while(line1.substr(60,20).compare("END OF FILE         ") != 0);

     f1.close();

     return 1;
}
//---------------------------------------------------------------------------
bool Class_TEC::CalcTecEpoch(double mjd,double fracOfDay,
                           double latuser,double longuser)
{/*-------------------------------------------------------------------
     Purpose: Find out in the GIM the nearest map of the position (lat, long)
              and MJD. Next step, interpolate the VTEC using bilinear interpolation
     -------------------------------------------------------------------
     Input:  mjd - modified julian date
                   fracOfDay - fraction of the day in mjd
                   latuser, longuser - user latitude and longitude

     Output: VTEC - Interpolated Vertical TEC
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              FAPESP PROCESS: 05/03522-1
     Date: July of 2010
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/

double  longini = Ini_Lon,           //initial long
        DLON0 = Lon_Increment,       //Increment in Longitude (in degrees)
        diflat0=0, diflat1=0,
        latgrid=0,
        DLON=0,
        hora_grid=0,dif_time0=0,dif_time1=0,
        longrid[75]={0.0},diflon0=0,diflon1=0,
        lat1=0,lat2=0,lon1=0,lon2=0,VTEC1=0,VTEC2=0,VTEC3=0,VTEC4=0,
        VTEC_Epo1=0,VTEC_Epo2=0,
        Ga=0,Gb=0;

int NLON = int((End_Lon-Ini_Lon)/Lon_Increment)+1, //number of longitudes in each epoch
    NLAT=int( (End_Lat-Ini_Lat) /Lat_Increment)+1; //number of latitudes in each epoch

int k=0;
bool time_gim = false;

    mjd += fracOfDay; //MJD + a fraction of the day (Tr in MJD)

    //Find out nearest time and grid of mjd
   while(k < Num_Value)
    {
        hora_grid = TECMAP[k][NLON+1];  //hour of the grid
        dif_time0 = mjd - hora_grid;    //hour difference in mjd

        if(dif_time0 < 0)
        { time_gim = true;
          break;
        }//if

        if(k==0)
          dif_time1 = dif_time0;        //only int he first time
        else if(dif_time1>dif_time0)
          dif_time1 = dif_time0;

        k = k + NLAT; //to indicate the position in the matrix TECMAP relative to the wanted mjd

    }

    if(time_gim == false)  //MJD in the GIM was not found
      return(false);

    k = k - NLAT; //get back 1 map

    if(k<0) k = 0;

    int start_mapa = k;  //start map found

    //Find out the nearest latitudess
    while(k < Num_Value)
    {
       latgrid = TECMAP[k][NLON];    //Grid Latitude
       diflat0 = latgrid - latuser;  //difference in latitude

       if(diflat0 < 0) break;

       if(k==start_mapa)              //at the first time
         diflat1 = diflat0;
       else if(diflat1>diflat0)
         diflat1 = diflat0;

       k = k + 1;                  //k refers to the position of latitude found
    }//while

    //Find out the nearest lonngitude
    int ijk=0;

    while(ijk<NLON)
    {
       longrid[ijk] =  longini + DLON;

       diflon0 = longuser - longrid[ijk];

       if(diflon0 < 0) break;

       if(ijk==0)
         diflon1 = diflon0;
       else if(diflon1>diflon0)
         diflon1 = diflon0;

       ijk += 1;

       DLON = DLON + DLON0;
    }//while

    //----------Fist Grid-----------------------------------
    lat1 = TECMAP[k-1][NLON];      //initial lat od the grid
    lat2 = TECMAP[k][NLON];        //final lat of the grid
    lon1 = longrid[ijk-1];         //initial long of the grid
    lon2 = longrid[ijk];           //final long of the grid

    VTEC1 = TECMAP[k-1][ijk-1];    //VTEC1
    VTEC2 = TECMAP[k-1][ijk];      //VTEC2
    VTEC3 = TECMAP[k][ijk-1];      //VTEC3
    VTEC4 = TECMAP[k][ijk];        //VTEC4

    Ga = Gb = VTEC_Epo1 = 0.0;

    Linear_Interp(VTEC1,VTEC2,lon1,lon2,longuser, Ga );
    Linear_Interp(VTEC3,VTEC4,lon1,lon2,longuser, Gb );
    Linear_Interp(Ga,Gb,lat1,lat2,latuser, VTEC_Epo1 );

    //----------Second Grid-----------------------------------
    lat1 = TECMAP[k-1+NLAT][NLON];   //initial lat od the grid
    lat2 = TECMAP[k+NLAT][NLON];     //final lat of the grid
    lon1 = longrid[ijk-1];           //initial long of the grid
    lon2 = longrid[ijk];             //final long of the grid

    VTEC1 = TECMAP[k-1+NLAT][ijk-1]; //VTEC1
    VTEC2 = TECMAP[k-1+NLAT][ijk];   //VTEC2
    VTEC3 = TECMAP[k+NLAT][ijk-1];   //VTEC3
    VTEC4 = TECMAP[k+NLAT][ijk];     //VTEC4

    //----------Interpolation-----------------------------------
    Ga = Gb = VTEC_Epo2 = 0.0;

    Linear_Interp(VTEC1,VTEC2,lon1,lon2,longuser, Ga );
    Linear_Interp(VTEC3,VTEC4,lon1,lon2,longuser, Gb );
    Linear_Interp(Ga,Gb,lat1,lat2,latuser, VTEC_Epo2 );

    double hora1 = TECMAP[k][NLON+1],
           hora2 = TECMAP[k+NLAT][NLON+1];

    //final interpolation
    Linear_Interp(VTEC_Epo1,VTEC_Epo2,hora1,hora2,mjd, VTEC );

 return true;

}
//---------------------------------------------------------------------------
void Class_TEC::Linear_Interp(double VTEC1,double VTEC2,
                              double coord1,double coord2,
                              double coord_user,double &VTEC )
{/*-------------------------------------------------------------------
    Purpose:  Linear interpolation
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             haroldoh2o@gmail.com

    Date: January of 2010
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
double Ga,Gb;

   if(VTEC1 >= VTEC2)
    VTEC = ( (VTEC1 - VTEC2)*(coord1 - coord_user)- (coord1 - coord2)*VTEC1 )
         / -(coord1 - coord2);
   else
    VTEC = ( (VTEC2 - VTEC1)*(coord2 - coord_user)- (coord2 - coord1)*VTEC2 )
         / -(coord2 - coord1);

}
//---------------------------------------------------------------------------
void Class_TEC::Sist_Local(double xest,double yest,double zest,
                          double xsat,double ysat,double zsat,
                          double fi, double lambda,
                          double &e,double &n,double &u)
{/*-------------------------------------------------------------------
     Purpose: Geodetic Cartesian Coordinates to Geodetic Local system
     -------------------------------------------------------------------
     Input:  xestacao, yestacao, zestacao - station coordinates

     Output:  e, n, u - coord. in the local system
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP
             Processo FAPESP: 05/03522-1

     Date: January of 2010
     -------------------------------------------------------------------
     Observation:
     -------------------------------------------------------------------*/
  double senfi,cosfi,senlamb,coslamb,senfi_coslamb,senfi_senlamb,cosfi_coslamb,
         cosfi_senlamb,lat,lamb,h,result,result1,dx,dy,dz;

   dx = xsat - xest ;
   dy = ysat - yest ;
   dz = zsat - zest ;

   senfi = sin(fi);
   cosfi = cos(fi);
   senlamb = sin(lambda);
   coslamb = cos(lambda);
   senfi_coslamb = senfi*coslamb;
   senfi_senlamb = senfi*senlamb;
   cosfi_coslamb = cosfi*coslamb;
   cosfi_senlamb = cosfi*senlamb;

   n = (-senfi_coslamb*dx) - (senfi_senlamb*dy) + (cosfi*dz);
   e = (-senlamb*dx)+(coslamb*dy);
   u = (cosfi_coslamb*dx) + (cosfi_senlamb*dy) + (senfi*dz);

}
//------------------------------------------------------------------------------------
void Class_TEC::Ang_Local(double E, double N, double U,
                         double Xsat,double Ysat,double Zsat,
                         double &Am, double &Em, double &Zm)
{/*-------------------------------------------------------------------
     Purpose: to compute Azimuth and elevation angle in the local system
     -------------------------------------------------------------------
     Input: E, N, U - Coord. in the local system

     Output: Am, Em - Azimuth and elevation angle
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
              Processo FAPESP: 05/03522-1

     Date: January of 2010
     -------------------------------------------------------------------
     Observation:  The equations can be found in (STRANG; BORRE, 1997, p. 502)

      STRANG, G.; BORRE, k. Linear Algebra, Geodesy and GPS. Wellesley-Cambrigde
      Press, 1997, 624p
     -------------------------------------------------------------------*/

double ro;


    //satellite-receiver distance in the local system
    ro = sqrt( pow(N,2)+pow(E,2)+pow(U,2) );

    //normalized vectors in the local system
    N = N/ro;
    E = E/ro;
    U = U/ro;

    //azimuth in the local system
    Am = atan(E/N);

    if(E>0 && N<0) Am = M_PI - fabs(Am);                  //2nd quadrant
    else if(E<0 && N<0) Am = Am + M_PI;                   //3rd quadrant
         else if(E<0 && N>0) Am = (2.0*M_PI) - fabs(Am);  //4th quadrant

    while(Am > (2.0*M_PI))
      Am -= (2.0*M_PI);
    while (Am < 0)
      Am += (2.0*M_PI);

    Em = asin(U);            //Elevation angle
    Zm = (M_PI/2.0) - Em;    //Zenith angle

return;
}
//------------------------------------------------------------------------------------
void Class_TEC::Pierce_Point_Coord(double Xest,double Yest,double Zest,
                                   double Xsat,double Ysat,double Zsat,
                                   double lat,double lamb, double h)
{/*-------------------------------------------------------------------
     Purpose: compute geodetic coordinates at the pierce point
     -------------------------------------------------------------------
     Input:

     Output:  lat_ion = latitude of pierce point
             lamb_ion = longitude of pierce point
     -------------------------------------------------------------------
     Authors: Haroldo Antonio Marques
              Programa de Pos-Graduacao em Ciencias Cartograficas
              FCT/UNESP - Presidente Prudente - SP
               Processo FAPESP: 05/03522-1

     Date: January of 2010
     -------------------------------------------------------------------
     Observation:

      References:

      GISAWY, M. L. Development of an Ionosphere Monitoring Technique Using GPS
      Measurements for High Latitude GPS Users. 2003. 161 p.
      Thesis. University of Calgary. Calgary.
      Disponível em: <http://www.geomatics.ucalgary.ca/links/GradTheses.html>
      Acesso em: mar. 2007

      MATSUOKA, M. T.; CAMARGO, P. O. Cálculo do TEC usando dados de receptores GPS
      de dupla frequencia par a produção de mapas da ionosfera para a região brasileira.
      revista Brasileira de Cartografia, n. 56/01 jul. 2004.
     -------------------------------------------------------------------*/
    double temp;
    double E,N,U;

    try{
          Sist_Local(Xest,Yest,Zest,Xsat,Ysat,Zsat,lat,lamb,E,N,U);

          Ang_Local(E,N,U,Xsat,Ysat,Zsat,Az,El,Zen);

          temp = (Re/(Re+hion)) * sin(Zen);

          zl = asin(temp);

          //pierce point latitude
          temp = sin(lat)*cos(Zen-zl) + cos(lat)*sin(Zen-zl)*cos(Az);
          lat_ion = asin(temp);

          //pierce point longitude
          temp = (sin(Zen - zl)*sin(Az))/cos(lat_ion);
          lamb_ion = lamb + asin(temp);

   }catch(...)
   {
      printf("Problems in the Pierce_Point_Coord (pierce point) subroutine");
      return;
   }
}
//------------------------------------------------------------------------------------

}//namespace

