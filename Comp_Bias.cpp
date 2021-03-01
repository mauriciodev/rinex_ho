/*-------------------------------------------------------------------
    Purpose:  To compute phase bias for all satellite arcs data in the rinex file
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques

    Date: Sept. of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
#include <iostream>
#include <ios>
#include <fstream>


#include "Comp_Bias.h"


//file to save data fo each satellite
FILE *PRN1,*PRN2,*PRN3,*PRN4,*PRN5,*PRN6,*PRN7,*PRN8,*PRN9,*PRN10,*PRN11,*PRN12,
        *PRN13,*PRN14,*PRN15,*PRN16,*PRN17,*PRN18,*PRN19,*PRN20,*PRN21,*PRN22,*PRN23,
        *PRN24,*PRN25,*PRN26,*PRN27,*PRN28,*PRN29,*PRN30,*PRN31,*PRN32;

void Comp_Bias(string filenameobs, Class_TEC Calc_Tec,MEAN_PHASE_BIAS Phase_Bias[][max_arc],
               ofstream &File_Cycleslip)
{ /*-------------------------------------------------------------------
    Purpose:  Read rinex file and compute phase bias for all satellite arc data
              The results fr each arc is stored in the matrix typedef struct Phase_Bias[PRNi][arc]
    -------------------------------------------------------------------
    Author: Haroldo Antonio Marques

    Date: Sept. of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/


 double ph1[MAXPRNID]={0.0},    //Phase - L1
        ph2[MAXPRNID]={0.0},    //Phase - L2
        ca[MAXPRNID]={0.0},     //code - c1
        p1[MAXPRNID]={0.0},     //code - p1
        p2[MAXPRNID]={0.0},     //code - p2
        c2[MAXPRNID]={0.0},     //code - L2C
        c5[MAXPRNID]={0.0},     //code C5 - L5
        ph5[MAXPRNID]={0.0};

 double STEC=0.0;

 double MJD, Interval;
 string theObsFileType,SatCode;  //type of observation

 int Isat[NSAT],NS,PRN,cont_epoch,Ns_temp,count_sat_out;
 int numobstype;
 bool Cycle_Slip_Flag=false,wlslip=false, nlslip=false,epoca_ini_cs,
       Sat_History[NSAT]={false},Curr_Sat_Epoch[NSAT]={false},
       Sat_Out[NSAT]={false},New_Sat[NSAT]={false};

 MEAN_PHASE_CODE Sat_Average_Phase_Code[NSAT]={0.0};


    Class_Cycle_Slip CycleSlip[NSAT];       //cycle slip class
    RinexFile *myFile = new RinexFile();   //object from class RinexFile

    OpenPRNFile();   //Open one file for each PRN

    //start structure Sat_Average_Phase_Code
    memset(Sat_Average_Phase_Code,0,sizeof(Sat_Average_Phase_Code));


    try {  //try to open rinex file
            myFile->setPathFilenameMode( filenameobs, ios::in );

            //Type of file...
            theObsFileType = myFile->getRinexFileType();
            delete myFile;
            myFile = 0;
        }//try
       catch( RinexFileException &openExcep )
        {
            cout << "Error opening file: " << filenameobs << endl
                 << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
                << "Exiting the program... " << endl << endl;

            getchar();
            exit (-1);
        }//catch

    if( theObsFileType[0] == 'O' ) //check if it is observation file
    {
        RinexObsFile myobs;        //RinexObsFile class
        RinexObsFile mynewobs;     //RinexObsFile class
        ObsEpoch  currentObsEpoch; //ObsEpoch class


        try
        {    //try to open observation file
             myobs.setPathFilenameMode(filenameobs,ios_base::in);
        }//try
        catch( RinexFileException &openExcep )
        {
            cout << "Error opening file: " << filenameobs << endl
                 << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
                 << "Exiting the program... " << endl << endl;

            exit (-1);
        }//catch

        //Read Header
         myobs.readHeader();

         Interval = myobs.getObsInterval(); //interval of data

        while( myobs.readEpoch( currentObsEpoch ) != 0 )   //read epochs
        {
            if(cont_epoch==0)   //To indicate initial epoch
              epoca_ini_cs = true;
            else
              epoca_ini_cs = false;

            //GPS Reception time in MJD
            MJD = currentObsEpoch.getEpochTime().GetMJD().mjd + currentObsEpoch.getEpochTime().GetMJD().fracOfDay;

            int hour = currentObsEpoch.getEpochTime().GetYMDHMS().hour;
            int min =  currentObsEpoch.getEpochTime().GetYMDHMS().min;
            double sec =  currentObsEpoch.getEpochTime().GetYMDHMS().sec;

            //set 0 to  GPS obs. vector
            memset (ca,0,sizeof(ca));
            memset (p1,0,sizeof(p1));
            memset (c2,0,sizeof(c2));
            memset (c5,0,sizeof(c5));
            memset (ph1,0,sizeof(ph1));
            memset (ph2,0,sizeof(ph2));
            memset (Isat,0,sizeof(Isat));


            //number of satelites from file
            NS =  (unsigned short) currentObsEpoch.getNumSat();
            Ns_temp = NS;
            count_sat_out=0;

            for (int i = 0; i < NS; i++ )
            {
                 SatCode = "";
                 SatCode = currentObsEpoch.getSatListElement(i).satCode;  //Satellite Code

                 int ptr = strcmp(SatCode.c_str(),"G"); //G for GPS satellites
                 int ptr1 = strcmp(SatCode.c_str()," ");

                 if( (ptr == 0 || ptr1 ==0) && currentObsEpoch.getSatListElement(i).satNum != 9999)//Only GPS satellites
                 {
                     int pos = i-count_sat_out;
                     Isat[pos] = currentObsEpoch.getSatListElement(i).satNum;


                     numobstype = myobs.getNumObsTypes();
                     //Store observations in arrays
                     for( int j = 0; j < numobstype; j++ )
                     {
                         if ( j != 0 &&  j % 5 == 0 )  // Max. 5 obs. per line
                         {
                             cout << endl;
                         }
                         if ( currentObsEpoch.getSatListElement(pos).obsList[ j ].obsPresent ) //Check if available observable
                         {
                             int var = currentObsEpoch.getSatListElement(pos).obsList[j].obsType;

                             switch(var)
                             { /* In the rinex.h file you can find enum variable
                                  enum OBSTYPE { NOOBS = 0, L1 = 1, L2 = 2, C1 = 3, P1 = 4, P2 = 5,
                                                 C2 = 6, D1 = 7, D2 = 8, T1 = 9, T2= 10, S1 = 11, S2 = 12,
                                                 L5 =13, C5 = 14, D5 = 15, S5 = 16  };
                               */
                                 case 1: ph1[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;   break;
                                 case 2: ph2[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;   break;
                                 case 3: ca[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;
                                 case 4: p1[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;
                                 case 5: p2[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;
                                 case 6: c2[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;
                                 //...
                                 case 13: ph5[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;
                                 case 14: c5[pos] = currentObsEpoch.getSatListElement(pos).obsList[j].observation;    break;

                                 default: break;
                             }//switch
                          }//if ( currentObsEpoch.getSatListElement(pos).obsList[ j ].obsPresent )
                     }//for( int j = 0; j < numobstype; j++ )
                 }//if( (ptr == 0 || ptr1 ==0) && currentObsEpoch.getSatListElement(i).satNum != 9999)
                 else
                 {
                     Ns_temp--;
                     count_sat_out++;
                     //...
                     //if you want to store GLONASS for example insert your code here

                  }//else     */
            }//for

            NS = Ns_temp; //number of GPS satelite

            //-----------Check if the satellite stoped to be tracked-----------------------------
            if(!epoca_ini_cs)
            {
               bool flag_sat_out = false;

               //Check out if satellite stop to be tracked
               Check_Sat_Out( NS,Isat,Curr_Sat_Epoch,Sat_History,Sat_Out,flag_sat_out);

               if(flag_sat_out)
                for(int i=0;i<NSAT;i++)
                  if(Sat_Out[i])   //means that satellite stop to be tracked (phase break)
                  {  PRN = i;

                    if(Sat_Average_Phase_Code[PRN].Count_Mean > 0)
                      if(Sat_Average_Phase_Code[PRN].Start_Time)
                      {   Sat_Average_Phase_Code[PRN].End_Time = MJD;

                          //Print data in the specific satellite file
                          Print_Data_PRN(PRN, &Sat_Average_Phase_Code[PRN]);

                          cout<<PRN<<" "<<Sat_Average_Phase_Code[PRN].Start_Time<<" "
                                   <<Sat_Average_Phase_Code[PRN].End_Time<<" "
                                   <<Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi<<endl;

                          //Store data for the arc
                          int arc = Sat_Average_Phase_Code[PRN].arc_number-1;
                          Phase_Bias[PRN][arc].Mean_Phase_Bias = Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi;
                          Phase_Bias[PRN][arc].Start_Time = Sat_Average_Phase_Code[PRN].Start_Time;
                          Phase_Bias[PRN][arc].End_Time = Sat_Average_Phase_Code[PRN].End_Time;
                          Phase_Bias[PRN][arc].Count_Mean = Sat_Average_Phase_Code[PRN].Count_Mean;
                          Phase_Bias[PRN][arc].arc_number = Sat_Average_Phase_Code[PRN].arc_number;

                          Sat_Average_Phase_Code[PRN].Start_Time = 0.0;
                          Sat_Average_Phase_Code[PRN].End_Time = 0.0;
                          Sat_Average_Phase_Code[PRN].Count_Mean = 0;  //force to start filter to compute phase bias
                      }//if

                      Sat_Out[i]=false;
                  }
            }//if


            //Monitoring in and out of satellites
            memset (Curr_Sat_Epoch,false,sizeof(Curr_Sat_Epoch));

            for(int i=0;i<NS;i++)
            {
                 PRN = Isat[i];

                 if(PRN !=-1)
                   Curr_Sat_Epoch[PRN] = true; //store satellites of the epoch which will take part of the processing
            }

            if(epoca_ini_cs==true )
            {
                 for(int i=0;i<NS;i++)
                 {  PRN = Isat[i];

                    if(PRN !=-1)
                      Sat_History[PRN] = true; //store visible satellite in the first epoch
                 }//for
            }//if
            else{

                  bool flag_sat_in = false;

                  //Check if the satellite is satarting to be tracked
                  Check_New_Sat( NS, Isat, Curr_Sat_Epoch, Sat_History, New_Sat,flag_sat_in);

                  if(flag_sat_in)
                     for(int i=0;i<NSAT;i++)
                       if(New_Sat[i])   //means that satellite start to be tracked (phase break)
                       {   PRN = i;

                           CycleSlip[PRN].Start_CycleSlip();

                           //start structure to compute phase bias
                           Sat_Average_Phase_Code[PRN].Start_Time = 0.0;
                           Sat_Average_Phase_Code[PRN].End_Time = 0.0;
                           Sat_Average_Phase_Code[PRN].Count_Mean = 0;
                           Sat_Average_Phase_Code[PRN].Last_Diff_Li_Pi = 0.0;
                           Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi = 0.0;

                           New_Sat[PRN] = false;  //start vector for new epoch
                       }//if
            }//else


           //-------DCB corrections for P1-----------------------------------
            //Applying DCB P1-C1 to compatibilize C1 with P1
            for(int i=0;i<NS;i++)
            {
               int PRN = Isat[i];
               if(PRN && p1[i]==0.0 && ca[i]!=0.0)
               {
                   double p1c1 = Calc_Tec.Get_P1C1(PRN);
                   p1[i] = ca[i]+c*p1c1*1.0e-9;
               }
            }//for

            for (int i = 0; i < NS; i++ )
            {   //int PRN = currentObsEpoch.getSatListElement(i).satNum;
                PRN = Isat[i];
                if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                {
                   // cout<<"Computing cycle slip..."<<endl;

                    wlslip = nlslip = false;
                    CycleSlip[PRN].Wide_Cycle_Slip( ph1[i],  //PHL1[PRN],
                                                    ph2[i],  //PHL2[PRN],
                                                    p1[i],   //p1[PRN],
                                                    p2[i],   //P2[PRN],
                                                    epoca_ini_cs);

                      CycleSlip[PRN].WN_Cycle_Slip( ph1[i],  //PHL1[PRN],
                                                    ph2[i],  //PHL2[PRN],
                                                    p1[i],   //p1[PRN],
                                                    p2[i],   //P2[PRN],
                                                    epoca_ini_cs);

                     File_Cycleslip<<setw(2)<<setfill('0')<<hour<<":"
                                    <<setw(2)<<setfill('0')<<min<<":"
                                    <<setw(4)<<setfill('0')<<sec;

                     File_Cycleslip<<setfill(' ');

                     for(int j=1;j<=32;j++)
                     {
                         int PRN = j;

                         wlslip = CycleSlip[PRN].Get_flag_wl();   //flag for wide lane
                         nlslip = CycleSlip[PRN].Get_flag_nl();   //flag for narrow lane

                         if(wlslip==true)
                             File_Cycleslip<<setw(13)<<setprecision(6)<<CycleSlip[PRN].Get_diff_wl();
                         else if(nlslip==true)
                             File_Cycleslip<<setw(13)<<setprecision(6)<<CycleSlip[PRN].Get_diff_nl();
                         else  File_Cycleslip<<setw(13)<<setprecision(6)<<9999999;
                     }

                     File_Cycleslip<<endl;
                } //if
            }//for

            for (int i = 0; i < NS; i++ )
            {   //int PRN = currentObsEpoch.getSatListElement(i).satNum;
                PRN = Isat[i];
                if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0){

                //PRN = Isat[i];
                if(Sat_Average_Phase_Code[PRN].Start_Time==0)
                   Sat_Average_Phase_Code[PRN].Start_Time = MJD;

                if( CycleSlip[PRN].Get_flag_wl()==true || CycleSlip[PRN].Get_flag_nl()==true ) //flag for wide lane or flag for narrow lane
                {  Cycle_Slip_Flag  = true;
                   CycleSlip[PRN].Start_CycleSlip();

                  if (Sat_Average_Phase_Code[PRN].Count_Mean > 0)  //if there are computed value for phase bias
                  {
                     Sat_Average_Phase_Code[PRN].End_Time = MJD;

                     //Print data in the specific satellite file
                     Print_Data_PRN(PRN, &Sat_Average_Phase_Code[PRN]);

                     cout<<PRN<<" "<<Sat_Average_Phase_Code[PRN].Start_Time<<" "<<Sat_Average_Phase_Code[PRN].End_Time<<" "
                         <<Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi<<endl;

                     int arc = Sat_Average_Phase_Code[PRN].arc_number-1;
                     Phase_Bias[PRN][arc].Mean_Phase_Bias = Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi;
                     Phase_Bias[PRN][arc].Start_Time = Sat_Average_Phase_Code[PRN].Start_Time;
                     Phase_Bias[PRN][arc].End_Time = Sat_Average_Phase_Code[PRN].End_Time;
                     Phase_Bias[PRN][arc].Count_Mean = Sat_Average_Phase_Code[PRN].Count_Mean;
                     Phase_Bias[PRN][arc].arc_number = Sat_Average_Phase_Code[PRN].arc_number;

                     Sat_Average_Phase_Code[PRN].Start_Time = 0.0;
                     Sat_Average_Phase_Code[PRN].End_Time = 0.0;
                  }//if
                }
                else
                  Cycle_Slip_Flag=false;

                Calc_Tec.Mean_Phase_Code(ph1[i],ph2[i],p1[i],p2[i],
                                         Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi,
                                         Sat_Average_Phase_Code[PRN].Last_Diff_Li_Pi,
                                         Sat_Average_Phase_Code[PRN].Count_Mean,
                                         Sat_Average_Phase_Code[PRN].arc_number,
                                         cont_epoch,Cycle_Slip_Flag,PRN);

            } //for



        }//while( myobs.readEpoch( currentObsEpoch ) != 0 )   //read epochs


           cont_epoch++; //counter of epochs

           currentObsEpoch.initializeData();
      }

    }//if( theObsFileType[0] == 'O' )


    for(int i=1;i<=NSAT;i++)
    {   PRN = i;

       if(Sat_Average_Phase_Code[PRN].Count_Mean)
         if(Sat_Average_Phase_Code[PRN].Start_Time)
         {  Sat_Average_Phase_Code[PRN].End_Time = MJD;

            Sat_Average_Phase_Code[PRN].End_Time += Interval/86400.0;

            int arc = Sat_Average_Phase_Code[PRN].arc_number-1;
            Phase_Bias[PRN][arc].Mean_Phase_Bias = Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi;
            Phase_Bias[PRN][arc].Start_Time = Sat_Average_Phase_Code[PRN].Start_Time;
            Phase_Bias[PRN][arc].End_Time = Sat_Average_Phase_Code[PRN].End_Time;
            Phase_Bias[PRN][arc].Count_Mean = Sat_Average_Phase_Code[PRN].Count_Mean;
            Phase_Bias[PRN][arc].arc_number = Sat_Average_Phase_Code[PRN].arc_number;

             //Print data in the specific satellite file
            Print_Data_PRN(PRN, &Sat_Average_Phase_Code[PRN]);

            cout<<PRN<<" "<<Sat_Average_Phase_Code[PRN].Start_Time<<" "
                <<Sat_Average_Phase_Code[PRN].End_Time<<" "
                <<Sat_Average_Phase_Code[PRN].Mean_Diff_Li_Pi<<endl;

         }
    }//for

}//end of function
//---------------------------------------------------------------------------
double Get_Phase_Bias(int PRN, MEAN_PHASE_BIAS Phase_Bias[][max_arc], double MJD)
{
   int arc = Phase_Bias[PRN][0].arc_number;
   int cont =0;
   while(arc != 0)
   {
      if( MJD >= Phase_Bias[PRN][cont].Start_Time && MJD < Phase_Bias[PRN][cont].End_Time )
      {
         if(Phase_Bias[PRN][cont].Count_Mean >= 1)
          return(Phase_Bias[PRN][cont].Mean_Phase_Bias);
      }
      cont++;
      arc = Phase_Bias[PRN][cont].arc_number;

   }

 return (0.0);
}
//---------------------------------------------------------------------------
void Check_Sat_Out( int NS,
                    int isat[],
                    bool Sat_Epoch_Ant[],
                    bool Sat_History[],
                    bool SatOut[],
                    bool &flag_sat_out)
{/*-------------------------------------------------------------------
    Purpose: Check if satellite stop to be observed
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
  int sat_out,i,j;
  bool ok;


    for(i=0;i<NSAT;i++)
    {
        if(Sat_Epoch_Ant[i] == 1)
        {  //Chek if satellite was tracked
           j=0;
           ok=false;
           while(j<NS)
           {
              if(isat[j] == i)
              { ok = true;
                break;
              }//if
              j++;
           }//while

           if(!ok)
           {  //Reorg vectors

              sat_out = i;
              Sat_Epoch_Ant[sat_out] = 0;
              Sat_History[sat_out] = 0;

              SatOut[i]=true;   //store PRN satellite getting out

              flag_sat_out=true;  //flag to indicates satellite out

           }//if (!ok)
        }//if
    }//for

return;

}
//---------------------------------------------------------------------------
void Check_New_Sat( int NS,int isat[],bool Sat_Epoch[],bool Sat_History[],bool New_Sat[],bool &flag_sat_in)
{/*-------------------------------------------------------------------
    Purpose: Check if a new satellite started to be observed
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques

    Date: Sept of 2013
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/


    memset(New_Sat,false,sizeof(New_Sat));

    for(int i=0;i<NS;i++)
    {  if(isat[i]!=-1)
       {  int PRN = isat[i];

          if( Sat_History[PRN] == 0 && Sat_Epoch[PRN] == 1 ) //indicates new satellite
          {
              Sat_History[PRN] = Sat_Epoch[PRN];

              New_Sat[PRN]=true;    //store PRN satellite getting in

              flag_sat_in=true;  //flag to indicates satellite in
          }//if( Sat_History[pos] == 0 && Sat_Epoch_Ant[pos] == 1 )
       }//if(isat[i]!=-1)
    }//for


return;
}
//---------------------------------------------------------------------------
bool directory_exists(string pathname)
{

   struct stat sb;

   if (stat(pathname.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
   {

      printf("Directory %s already exist\n",pathname.c_str());
      return true; //it is a directory...
   }

  return false;
}
//---------------------------------------------------------------------------
int Create_Dir(string Directory)
{
    int check = mkdir(Directory.c_str());

    if (!check)
    {
         printf("Directory created: %s \n",Directory.c_str());
         return 1;
    }
    else
    {
        printf("Unable to create directory: %s\n",Directory.c_str());
        return 0;
    }

}
//---------------------------------------------------------------------------
int OpenPRNFile()
{

   if(! directory_exists(".\\PhaseBias"))
        Create_Dir(".\\PhaseBias");


  if ((PRN1 = fopen(".\\PhaseBias\\PRN1.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN1.dat");
        exit (-1);
   }//if
   if ((PRN2 = fopen(".\\PhaseBias\\PRN2.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN2.dat");
        exit (-1);
   }//if
   if ((PRN3 = fopen(".\\PhaseBias\\PRN3.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN3.dat");
        exit (-1);
   }//if
   if ((PRN4 = fopen(".\\PhaseBias\\PRN4.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN4.dat");
        exit (-1);
   }//if

   if ((PRN5 = fopen(".\\PhaseBias\\PRN5.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN5.dat");
        exit (-1);
   }//if
   if ((PRN6 = fopen(".\\PhaseBias\\PRN6.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN6.dat");
        exit (-1);
   }//if
   if ((PRN7 = fopen(".\\PhaseBias\\PRN7.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN7.dat");
        exit (-1);
   }//if
   if ((PRN8 = fopen(".\\PhaseBias\\PRN8.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN8.dat");
        exit (-1);
   }//if
   if ((PRN9 = fopen(".\\PhaseBias\\PRN9.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN9.dat");
        exit (-1);
   }//if
   if ((PRN10 = fopen(".\\PhaseBias\\PRN10.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN10.dat");
        exit (-1);
   }//if
   if ((PRN11 = fopen(".\\PhaseBias\\PRN11.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN11.dat");
        exit (-1);
   }//if
   if ((PRN12 = fopen(".\\PhaseBias\\PRN12.dat", "w+"))== NULL) {
           fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN12.dat");
           exit (-1);
   }//if
   if ((PRN13 = fopen(".\\PhaseBias\\PRN13.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN13.dat");
        exit (-1);
   }//if
   if ((PRN14 = fopen(".\\PhaseBias\\PRN14.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN14.dat");
        exit (-1);
   }//if
   if ((PRN15 = fopen(".\\PhaseBias\\PRN15.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN15.dat");
        exit (-1);
   }//if
   if ((PRN16 = fopen(".\\PhaseBias\\PRN16.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN16.dat");
        exit (-1);
   }//if
   if ((PRN17 = fopen(".\\PhaseBias\\PRN17.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN17.dat");
        exit (-1);
   }//if
   if ((PRN18 = fopen(".\\PhaseBias\\PRN18.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN18.dat");
        exit (-1);
   }//if
   if ((PRN19 = fopen(".\\PhaseBias\\PRN19.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN19.dat");
        exit (-1);
   }//if
   if ((PRN20 = fopen(".\\PhaseBias\\PRN20.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN20.dat");
        exit (-1);
   }//if
   if ((PRN21 = fopen(".\\PhaseBias\\PRN21.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN21.dat");
        exit (-1);
   }//if
   if ((PRN22 = fopen(".\\PhaseBias\\PRN22.dat", "w+"))== NULL) {
           fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN22.dat");
           exit (-1);
   }//if
   if ((PRN23 = fopen(".\\PhaseBias\\PRN23.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN23.dat");
        exit (-1);
   }//if
   if ((PRN24 = fopen(".\\PhaseBias\\PRN24.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN24.dat");
        exit (-1);
   }//if

   if ((PRN25 = fopen(".\\PhaseBias\\PRN25.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN25.dat");
        exit (-1);
   }//if
   if ((PRN26 = fopen(".\\PhaseBias\\PRN26.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN26.dat");
        exit (-1);
   }//if
   if ((PRN27 = fopen(".\\PhaseBias\\PRN27.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN27.dat");
        exit (-1);
   }//if
   if ((PRN28 = fopen(".\\PhaseBias\\PRN28.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN28.dat");
        exit (-1);
   }//if
   if ((PRN29 = fopen(".\\PhaseBias\\PRN29.dat", "w+"))== NULL) {
        fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN29.dat");
        exit (-1);
   }//if
   if ((PRN30 = fopen(".\\PhaseBias\\PRN30.dat", "w+"))== NULL) {
          fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN30.dat");
          exit (-1);
   }//if
   if ((PRN31 = fopen(".\\PhaseBias\\PRN31.dat", "w+"))== NULL) {
          fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN31.dat");
          exit (-1);
   }//if
   if ((PRN32 = fopen(".\\PhaseBias\\PRN32.dat", "w+"))== NULL) {
          fprintf(stderr, "Cannot open output file %s.\n",".\\PhaseBias\\PRN32.dat");
          exit (-1);
   }//if
}
//---------------------------------------------------------------------------
void Close_PRN_Files()
{
    if(PRN1 != NULL) fclose(PRN1);
    if(PRN2 !=  NULL) fclose(PRN2);
    if(PRN3 != NULL) fclose(PRN3);
    if(PRN4 != NULL) fclose(PRN4);
    if(PRN5 != NULL) fclose(PRN5);
    if(PRN6 != NULL) fclose(PRN6);
    if(PRN7 != NULL) fclose(PRN7);
    if(PRN8 != NULL) fclose(PRN8);
    if(PRN9 != NULL) fclose(PRN9);
    if(PRN10 != NULL) fclose(PRN10);
    if(PRN11 != NULL) fclose(PRN11);
    if(PRN12 != NULL) fclose(PRN12);
    if(PRN13 != NULL) fclose(PRN13);
    if(PRN14 != NULL) fclose(PRN14);
    if(PRN15 != NULL) fclose(PRN15);
    if(PRN16 != NULL) fclose(PRN16);
    if(PRN17 != NULL) fclose(PRN17);
    if(PRN18 != NULL) fclose(PRN18);
    if(PRN19 != NULL) fclose(PRN19);
    if(PRN20 != NULL) fclose(PRN20);
    if(PRN21 != NULL) fclose(PRN21);
    if(PRN22 != NULL) fclose(PRN22);
    if(PRN23 != NULL) fclose(PRN23);
    if(PRN24 != NULL) fclose(PRN24);
    if(PRN25 != NULL) fclose(PRN25);
    if(PRN26 != NULL) fclose(PRN26);
    if(PRN27 != NULL) fclose(PRN27);
    if(PRN28 != NULL) fclose(PRN28);
    if(PRN29 != NULL) fclose(PRN29);
    if(PRN30 != NULL) fclose(PRN30);
    if(PRN31 != NULL) fclose(PRN31);
    if(PRN32 != NULL) fclose(PRN32);

return;
}
//---------------------------------------------------------------------------
void Print_Data_PRN(int PRN, MEAN_PHASE_CODE *Sat_Average_Phase_Code)
{

   switch(PRN)
   {
          case 	1: fprintf(PRN1,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN1,"\n");
                   break;
          case 	2: fprintf(PRN2,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN2,"\n");
                   break;

          case 	3: fprintf(PRN3,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN3,"\n");
                   break;

          case 	4: fprintf(PRN4,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN4,"\n");
                   break;

          case 	5: fprintf(PRN5,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN5,"\n");
                   break;

          case 	6: fprintf(PRN6,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN6,"\n");
                   break;

          case 	7: fprintf(PRN7,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN7,"\n");
                   break;

          case 	8: fprintf(PRN8,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN8,"\n");
                   break;

          case 	9: fprintf(PRN9,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN9,"\n");
                   break;

          case 	10: fprintf(PRN10,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN10,"\n");
                   break;

          case 	11: fprintf(PRN11,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN11,"\n");
                   break;

          case 	12: fprintf(PRN12,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN12,"\n");
                   break;

          case 	13: fprintf(PRN13,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN13,"\n");
                   break;

          case 	14: fprintf(PRN14,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN14,"\n");
                   break;

          case 	15: fprintf(PRN15,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN15,"\n");
                   break;

          case 	16: fprintf(PRN16,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN16,"\n");
                   break;

          case 	17: fprintf(PRN17,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN17,"\n");
                   break;

          case 	18: fprintf(PRN18,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN18,"\n");
                   break;

          case 	19: fprintf(PRN19,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN19,"\n");
                   break;

          case 	20: fprintf(PRN20,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN20,"\n");
                   break;

          case 	21: fprintf(PRN21,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN21,"\n");
                   break;

          case 	22: fprintf(PRN22,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN22,"\n");
                   break;

          case 	23: fprintf(PRN23,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN23,"\n");
                   break;

          case 	24: fprintf(PRN24,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN24,"\n");
                   break;

          case 	25: fprintf(PRN25,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN25,"\n");
                   break;

          case 	26: fprintf(PRN26,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN26,"\n");
                   break;

          case 	27: fprintf(PRN27,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN27,"\n");
                   break;

          case 	28: fprintf(PRN28,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN28,"\n");
                   break;

          case 	29: fprintf(PRN29,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN29,"\n");
                   break;

          case 	30: fprintf(PRN30,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN30,"\n");
                   break;

          case 	31: fprintf(PRN31,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN31,"\n");
                   break;

          case 	32: fprintf(PRN32,"%6d %15.4lf %15.4lf %15.4lf %5d ",Sat_Average_Phase_Code->arc_number,Sat_Average_Phase_Code->Start_Time,
                                          Sat_Average_Phase_Code->End_Time,
                                          Sat_Average_Phase_Code->Mean_Diff_Li_Pi,Sat_Average_Phase_Code->Count_Mean);
                   fprintf(PRN32,"\n");
                   break;
   }

}
//---------------------------------------------------------------------------
