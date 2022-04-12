
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

    Reviewes and updates by Steve Hilla
    -------------------------------------------------------------------
    Observation: Use C++ classes to read and save RINEX file:
                 http://www.ngs.noaa.gov/gps-toolbox/rinex.htm
                 The authors would like to thanks CAPES, FAPESP by financial support
    -------------------------------------------------------------------*/

//c::1304.13, SAH, change path, use one directory for all *.cpp/*.h files
//c::1304.13, SAH, comment out:  printf("\033[20;4f") and system()
//c::1304.17, SAH, write cycle slip detection info to the log file.

//MHA, I've changed some system("pause") by getchar();
//insert a file to print cycle slip

//comment out: FilterObs[prn].Set_cont_epoch(0); near line 1023



#include "rinex_ho.h"

int main(int argc, char* argv[]) {

//    #include "Include/Variable.inc"
    #include "Variable.h"
    string SatCode;
    int Isat[NSAT]={0},NS,Ns_temp(0);
    int count_sat_out(0),Pos_GPS_CurrObs[NSAT]={0};
    bool wlslip, nlslip; //added by SH
    double dx(0.0),dy(0.0),dz(0.0),dsv_sq(0.0),sv_range(0.0),omega_e(0.0),propagationDelay(0.0);
    double dtom(0.0),usx(0.0),usy(0.0),svx(0.0),svy(0.0),svz(0.0);
    double sh_az(0.0), sh_elv(0.0), sh_dist(0.0), sh_pplat(0.0), sh_pplon(0.0), sh_zp(0.0);

    //lat. and long. of IPP in degrees
    double lat_IPP_deg,lon_IPP_deg;

    MEAN_PHASE_BIAS Phase_Bias[NSAT][max_arc]={0.0};

//    memset(Phase_Bias,0,sizeof(Phase_Bias));   //set 0 to  Phase_Bias struct

//    system("clear");  //for linux
    //clrscr();           //for windows

    //check-out input parameters
    if(argc == 1)
      InpFileName = "rinex_ha.inp";
    else
       InpFileName = argv[1];


    // InpFileName = "rinex_ha_brft_full.inp";  //TEST
    // InpFileName = "rinex_ha_ASC1_full.inp";

    //Imput parameters from rinex_ha.inp
    //reading input parameters
    ReadInput( InpFileName,filenameobs,filenamenav,newobsfile,
               outfile,X0,Y0,Z0,Mask_Ele,read_tec,save_file,GIM_File,dcbr,
               arqu_p1c1_sat,arqu_p1p2_sat,Geomag_Type,Name_IGRF_Model);

    //Read data from rinex_ha_param.dat file
    ReadParam(sigca,sigp2,sigph1,sigph2,h_ion,R_e,B_eq,Ne_max);

    //object of class T_TEC
    Class_TEC Calc_Tec(sigca,sigp2,sigph1,sigph2,dcbr,0.0);

    double var_ca = pow(sigca,2),
           var_p2 = pow(sigp2,2),
           var_ph1 = pow(sigph1,2),
           var_ph2 = pow(sigph2,2);

    //Intial latitude and longitude for north pole
    double latpn = 79.8,
           lonpn = 288.1;

    //object from Class_Iono class
    Class_Iono calc_iono( latpn,lonpn,
                          h_ion, R_e, B_eq, Ne_max);

    //-----------Opening Files--------------------------------------------------
    fstream log; //log file
    log.open("rinex.log",ios::out);
    if( log.fail() )
    {
        cout << "It was not possible to open the file: rinex.log!" << endl;
        getchar();
        exit (-1);
    }//if

    string temp = outfile+"_tec.txt";
    ofstream tecfile(temp.c_str(),ios::out);
    if( tecfile.fail() )
    {
        cout << "It was not possible to open the file: tec.txt!" << endl;
        getchar();
        exit (-1);
    }//if
    else
    {
      tecfile.setf(ios::fixed);
      tecfile.setf(ios::showpoint);
      tecfile.setf(ios::right);
    }

    temp = outfile+"_I1L1.txt";
    temp = outfile+"_I2L1.txt";
    ofstream fstionoI2l1(temp.c_str());
    if( fstionoI2l1.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if

    temp = outfile+"_I2L2.txt";
    ofstream fstionoI2l2(temp.c_str());
    if( fstionoI2l2.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if

    temp = outfile+"_I3L1.txt";
    ofstream fstionoI3l1(temp.c_str());
    if( fstionoI3l1.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if

    temp = outfile+"_I3L2.txt";
    ofstream fstionoI3l2(temp.c_str());
    if( fstionoI3l2.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if


    //Satellite elevation angle
    temp = outfile+"_EleAngleSat.txt";   //MHA (may of 2013):File to print Satellite Elevation Angle
    ofstream File_sateleangle(temp.c_str());

    if( File_sateleangle.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if

    //MHA (may of 2013):File to print Cycle Slip
    temp = outfile+"_CycleSlip.txt";
    ofstream File_Cycleslip(temp.c_str());
    if( File_Cycleslip.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        getchar();
        exit (-1);
    }//if

    //---------------Configure output files-------------------------------------
    fstionoI2l1.setf(ios::fixed);       fstionoI2l2.setf(ios::fixed);
    fstionoI2l1.setf(ios::showpoint);   fstionoI2l2.setf(ios::showpoint);
    fstionoI2l1.setf(ios::right);       fstionoI2l2.setf(ios::right);

    fstionoI3l1.setf(ios::fixed);       fstionoI3l2.setf(ios::fixed);
    fstionoI3l1.setf(ios::showpoint);   fstionoI3l2.setf(ios::showpoint);
    fstionoI3l1.setf(ios::right);       fstionoI3l2.setf(ios::right);

    File_sateleangle.setf(ios::fixed);
    File_sateleangle.setf(ios::showpoint);
    File_sateleangle.setf(ios::right);

    File_Cycleslip.setf(ios::fixed);
    File_Cycleslip.setf(ios::showpoint);
    File_Cycleslip.setf(ios::right);

    if(save_file)  //MHA (may of 2013): Insert the print header file here
    {
        if(cont_epoch == 0)
        {   //In the first epoch print the output file headers
            string nome = filenameobs.substr(0,8);

            tecfile<<setw(8)<<"h:m:s";
            fstionoI2l1<<setw(8)<<"h:m:s";
            fstionoI2l2<<setw(8)<<"h:m:s";
            fstionoI3l1<<setw(8)<<"h:m:s";
            fstionoI3l2<<setw(8)<<"h:m:s";
            File_sateleangle<<setw(8)<<"h:m:s";
            File_Cycleslip<<setw(8)<<"h:m:s";

            for(int j=1;j<=32;j++)
            {
                if(j<10)
                {   fstionoI2l1<<setw(12)<<"PRN0"<<j;
                    fstionoI2l2<<setw(12)<<"PRN0"<<j;
                    fstionoI3l1<<setw(12)<<"PRN0"<<j;
                    fstionoI3l2<<setw(12)<<"PRN0"<<j;
                    tecfile<<setw(12)<<"PRN0"<<j;
                    File_sateleangle<<setw(12)<<"PRN0"<<j;
                    File_Cycleslip<<setw(12)<<"PRN0"<<j;
                }//if
                else
                {   fstionoI2l1<<setw(11)<<"PRN"<<j;
                    fstionoI2l2<<setw(11)<<"PRN"<<j;
                    fstionoI3l1<<setw(11)<<"PRN"<<j;
                    fstionoI3l2<<setw(11)<<"PRN"<<j;
                    tecfile<<setw(11)<<"PRN"<<j;
                    File_sateleangle<<setw(11)<<"PRN"<<j;
                    File_Cycleslip<<setw(11)<<"PRN"<<j;
                }//else
            }//for

            tecfile<<endl;
            fstionoI2l1<<endl;
            fstionoI2l2<<endl;
            fstionoI3l1<<endl;
            fstionoI3l2<<endl;
            File_sateleangle<<endl;
            File_Cycleslip<<endl;

        } //if(cont_epoch == 0)

    }if(save_file)


    cout  <<"RINEX_HO - Program to read the RINEX, compute the TEC"<<endl
          <<"and correct GPS data from 2nd and 3rd order ionospheric efffects" << endl<<endl
          <<"Author: Haroldo Antonio Marques - 2009"<<endl
          <<"        Programa de Pos-Graduacao em Ciencias Cartograficas"<<endl
          <<"        Sao Paulo State University (Universidade Estadual Paulista) - FCT/UNESP"<<endl
          <<"        R. Roberto simonsen, 305 - Presidente Prudente - SP, Brazil"<<endl
          <<"        FAPESP Process: 05/03522-1"<<endl<<endl;

    log <<"RINEX_HO - Program to read the RINEX, compute the TEC"<<endl
         <<"and correct GPS data from 2nd and 3rd order ionospheric efffects" << endl<<endl
         <<"Author: Haroldo Antonio Marques - 2009"<<endl
         <<"        Programa de Pos-Graduacao em Ciencias Cartograficas"<<endl
         <<"        Sao Paulo State University (Universidade Estadual Paulista) - FCT/UNESP"<<endl
         <<"        R. Roberto simonsen, 305 - Presidente Prudente - SP, Brazil"<<endl
         <<"        FAPESP Process: 05/03522-1"<<endl<<endl;

        cout<<"*********************************************************************"<<endl<<endl;
        log<<"**********************************************************************"<<endl<<endl;

    cout<<"Navigation file: "<<filenamenav<<endl;
    log<<"Navigation file: "<<filenamenav<<endl;

    RinexFile *myFileNav = new RinexFile();   //class RinexFile

    if( read_tec != 2) //if not to use GIM, read DCBs
     {
         Calc_Tec.Read_DCB_Sat(arqu_p1c1_sat,arqu_p1p2_sat);

         log<<"DCBs read from "<<arqu_p1c1_sat<<" and "<<arqu_p1p2_sat<<": "<<endl;
         log<<"PRN "<<setw(12)<<"P1C1(ns)"<<setw(12)<<"P1P2(ns)"<<endl;
         for(int ik=0;ik<=32;ik++)
         {
        	 if(Calc_Tec.Get_P1C1(ik)!=0.0)
        	  log<<setw(2)<< ik <<" "
        	     <<setprecision(4)<<setw(12)<<Calc_Tec.Get_P1C1(ik)
        	     << " ";
        	 if(Calc_Tec.Get_P1P2(ik)!=0.0)
        	  log<<setprecision(4)<<setw(12)<<Calc_Tec.Get_P1P2(ik)<<endl;
         }
     }

    //Reading broadcast ephemeris file
    if(!ReadNavFile(currentPRNBlock,myFileNav,filenamenav,cont_efe_nav,log))
    {   getchar();
        exit(-1);
    }//if

    cout<<"Observation file: "<<filenameobs<<endl;
    log<<"Observation file: "<<filenameobs<<endl;

    if(read_tec==2)
    { cout<<"GIM File: "<<GIM_File<<endl;
      log<<"GIM File: "<<GIM_File<<endl;
    }


//------------------------------------------------------------------------------------------------------------
   //Computing average phase bias for each stellite arc.
   //The phase bias will be used to compute STEC from carrier levelled by code
   if(read_tec==3)
     Comp_Bias( filenameobs,Calc_Tec,Phase_Bias,File_Cycleslip);
//------------------------------------------------------------------------------------------------------------

    RinexFile *myFile = new RinexFile();   //object from class RinexFile

    try {  //try to open rinex file
            myFile->setPathFilenameMode( filenameobs, ios::in );
            log<<" - Version "<<myFile->getFormatVersion()<<endl;
            log<<"Satellite Sistem: "<<myFile->getSatSystem()<<endl;
            log<<"Rinex Type: "<<myFile->getRinexFileType()<<endl;
            log<<"Rinex Program: "<<myFile->getRinexProgram()<<endl;
            log<<"Created by Agency: "<<myFile->getCreatedByAgency()<<endl;

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

            log << "Error opening file: " << filenameobs << endl
                << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
                << "Exiting the program... " << endl << endl;

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

        try
        {   //try to open new observation file
            mynewobs.setPathFilenameMode(newobsfile,ios_base::out);
        }//try
        catch( RinexFileException &openExcep )
        {
             cout<< "Error opening file: " << newobsfile << endl
                << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
                << "Exiting the program... " << endl << endl;

             log<< "Error opening file: " << newobsfile << endl
                << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
                << "Exiting the program... " << endl << endl;

            exit (-1);
        }//catch

        try
        {   //Reading Header
            myobs.readHeader();

            if(X0!=0.0 && Y0!=0.0 && Z0!=0.0)
            {
                //Station Coordinates input by user
                xest = X0;
                yest = Y0;
                zest = Z0;
            }
            else
            {
                //Station Coordinates from RINEX header
                xest = myobs.getApproxX();
                yest = myobs.getApproxY();
                zest = myobs.getApproxZ();

            }

            //Cartesian to Geodetic Coordinates
            calc_iono.Cart_Curv(xest,yest,zest,lat,lamb,h);

            //Write out the header of new observations file
            myobs.writeHeaderImage( mynewobs.outputStream );

        }//try
        catch( RequiredRecordMissingException &headerExcep )
        {
            cout << " RequiredRecordMissingException is: " << endl
                 << headerExcep.getMessage() << endl;

            log << " RequiredRecordMissingException is: " << endl
                << headerExcep.getMessage() << endl;
        }//catch


        // read the observation epochs
        try
        {
            cout <<endl<< "Reading observation file and applying correction..." << endl;
            cout << endl;
            log <<endl<< "Reading observation file and applying correction..." << endl;
            log << endl;


            SatObsAtEpoch  tempSatObsAtEpoch[NSAT];
            int num = myobs.getNumObsTypes();

            if(read_tec==2)
            {
                //Reading GIM file
                Calc_Tec.ReadGIM(GIM_File);

                if(Calc_Tec.Get_Gim_Flag() == 0)
                {
                    cout<<"Problems Reading GIM File!"<<endl;
                    return 0;
                }//if
                else
                  cout<<"GIM File read succesfully!"<<endl;
            }


            //Decimal degrees
            double lat_deg = lat*180.0/calc_iono.Get_PI();
            double long_deg = lamb*180.0/calc_iono.Get_PI();


            while( myobs.readEpoch( currentObsEpoch ) != 0 )   //read epochs
            {
                //print number of epochs in the screen
                //For linux
                // printf("\033[20;4f");

                //for Windows
                //gotoxy(10,20);  //command not found in linux (#include <conio.h>)
                cout<<"Epoch: "<<cont_epoch<<"..."<<endl;

                if(cont_epoch==0)
                  epoca_ini_cs = true;
                else
                  epoca_ini_cs = false;

                num = myobs.getNumObsTypes();
                calc_iono.Set_year(0);
                calc_iono.Set_dayofyear(0);
                calc_iono.Set_hour(0);
                calc_iono.Set_min(0);
                calc_iono.Set_sec(0.0);

                calc_iono.Set_year(currentObsEpoch.getEpochTime().GetYMDHMS().year);            //year
                calc_iono.Set_month(currentObsEpoch.getEpochTime().GetYMDHMS().month);          //month
                calc_iono.Set_day(currentObsEpoch.getEpochTime().GetYMDHMS().day);              //day
                calc_iono.Set_dayofyear(currentObsEpoch.getEpochTime().GetYDOYHMS().dayOfYear); //day of year
                calc_iono.Set_hour(currentObsEpoch.getEpochTime().GetYMDHMS().hour);            //hours
                calc_iono.Set_min(currentObsEpoch.getEpochTime().GetYMDHMS().min);              //minutes
                calc_iono.Set_sec(currentObsEpoch.getEpochTime().GetYMDHMS().sec);              //seconds
                segweek = currentObsEpoch.getEpochTime().GetGPSTime().secsOfWeek;               //seconds of week

                //GPS Reception time in MJD
                double MJD = currentObsEpoch.getEpochTime().GetMJD().mjd,
                       fracOfDay = currentObsEpoch.getEpochTime().GetMJD().fracOfDay;

                //Start some arrays
                zera_vetor(TECM,NSAT);
                zera_vetor(i2m_l1,NSAT);
                zera_vetor(i2m_l2,NSAT);
                zera_vetor(i3m_l1,NSAT);
                zera_vetor(i3m_l2,NSAT);
                zera_vetor(EleAngSat,NSAT);
                zera_vetor(AzAngSat,NSAT);

                if(read_tec==2)
                {
                   //Check if GIM time match with Obs file
                   if(cont_epoch==0)
                    if( fabs( (MJD+fracOfDay) - Calc_Tec.Get_Epoch_First_MAP() ) > (2.0/24.0) ) //||
                        //fabs( (MJD+fracOfDay) - Calc_Tec.Get_Epoch_Last_Map()) > (2.0/24.0) )
                    {
                        cout<<"GIM file doesn't match with RINEX file!"<<endl;
                        return 0;
                    }//if

/*  SH, move this to a loop over satellite (need one IPP, and VTEC, for each satellite!)
                    VTEC = 0.0;

                    //Interpolate TEC for MJD+fracOfDay
                    if(! Calc_Tec.CalcTecEpoch(MJD,fracOfDay,lat_deg,long_deg) )
                    {
                         log<<"Without TEC information from GIM "
                            <<"in the epoch: "<<calc_iono.Get_hour()<<" "
                            <<calc_iono.Get_min()<<" "
                            <<calc_iono.Get_sec()<<endl;

                         return 0;

                     }//if
                     else
                     {
                         VTEC = Calc_Tec.GetVTEC();
                     }//else
*/

                }//if

                //Reception time in seconds of week
                tr =  segweek;

                //start obs arrays
                zera_vetor(ca,NSAT);
                zera_vetor(p1,NSAT);
                zera_vetor(p2,NSAT);
                zera_vetor(c2,NSAT);
                zera_vetor(c5,NSAT);
                zera_vetor(ph1,NSAT);
                zera_vetor(ph2,NSAT);
                zera_vetor(ph5,NSAT);

                zera_vetor(Isat,NSAT);
                zera_vetor(Pos_GPS_CurrObs,NSAT);

                int numobstype = myobs.getNumObsTypes();

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
                       Pos_GPS_CurrObs[pos]=i;
                       //Store observations in arrays

                       StoreObsVector( numobstype,
                                       currentObsEpoch,
                                       ph1,  //phase L1
                                       ph2,  //phase L2
                                       ca,   // C/A code
                                       p1,   //P1 code
                                       p2,   //P2 code
                                       c2,   //L2C code - not being used
                                       c5,   //C5 code  - not being used
                                       ph5,  //phase L5 - not being used
                                       pos);

                    }//if
                    else
                    {
                        Ns_temp--;
                        count_sat_out++;
                        //...
                        //if you want to store GLONASS for example insert your code here

                     }//else     */
                }//for

                NS = Ns_temp; //number of GPS satelite

                //-------DCB corrections for P1-----------------------------------
                //Added by MHA (may of 2013)
                //Applying DCB P1-C1 to compatibilize CA with P1
                //DCBs p1c1 are not more applied inside the subroutine Calc_STEC_PR
                for(int i=0;i<NS;i++){
                   int PRN = Isat[i];
                   if(PRN && p1[i]==0.0 && ca[i]!=0.0)
                   {
                       double p1c1 = Calc_Tec.Get_P1C1(PRN);
                       p1[i] = ca[i]+c*p1c1*1.0e-9;
                   }
                }//for

                /* We can do the same for L2C and P2, but p2c2 must be read from CODE files for example
                 for(int i=0;i<NS;i++){
                    int PRN = isat[i];
                    if(p2[i]==0.0 && c2[i]!=0.0)
                      p2[i] = c2[i]+c*p2c2[PRN]*1.0e-9;
                 }//for
                */


                //---------//check-out for cycle slips-------------------
                if(read_tec==1) //smothed pseudorange by phase
                {
                     for (int i = 0; i < NS; i++ )
                     {   //int PRN = currentObsEpoch.getSatListElement(i).satNum;
                         int PRN = Isat[i];

                         if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         {

                              wlslip = nlslip = FALSE;
                              //log << "------------------ " << PRN << " "
                              //    << currentObsEpoch.getEpochTime().GetYMDHMS().hour << ":"
                              //    << currentObsEpoch.getEpochTime().GetYMDHMS().min << ":"
                              //    << (int) currentObsEpoch.getEpochTime().GetYMDHMS().sec << endl;
                              //log << "== step1: " << endl;
                              //log << "wlbias: " << CycleSlip[PRN].Get_wlbias() << endl;
                              //log << "sig_wlbias: " << CycleSlip[PRN].Get_sig_wlbias() << endl;
                              //log << "diff_wl: " << CycleSlip[PRN].Get_diff_wl() << endl;
                              //log << "flag: " << CycleSlip[PRN].Get_flag() << endl;
                              //log << "count: " << CycleSlip[PRN].Get_count() << endl;

                              CycleSlip[PRN].Wide_Cycle_Slip( ph1[i],  //PHL1[PRN],
                                                              ph2[i],  //PHL2[PRN],
                                                              p1[i],   //p1[PRN],
                                                              p2[i],   //P2[PRN],
                                                              epoca_ini_cs);


                              //log << "== step2: " << endl;
                              //log << "wlbias: " << CycleSlip[PRN].Get_wlbias() << endl;
                              //log << "sig_wlbias: " << CycleSlip[PRN].Get_sig_wlbias() << endl;
                              //log << "diff_wl: " << CycleSlip[PRN].Get_diff_wl() << endl;
                              //log << "flag: " << CycleSlip[PRN].Get_flag_wl() << endl;
                              //wlslip = CycleSlip[PRN].Get_flag_wl();
                             // log << "count: " << CycleSlip[PRN].Get_count() << endl;

                            /*  if( CycleSlip[PRN].Get_flag_wl() )
                              {
                                   log << "::: PRN:" << setw(2) << setfill('0') << PRN << "  HMS: "
                                       << setw(2) << currentObsEpoch.getEpochTime().GetYMDHMS().hour << ":"
                                       << setw(2) << currentObsEpoch.getEpochTime().GetYMDHMS().min << ":"
                                       << setw(2) << (int) currentObsEpoch.getEpochTime().GetYMDHMS().sec << "  difWL: "
                                       << setw(9) << setprecision(2) << setfill(' ') << CycleSlip[PRN].Get_diff_wl() << endl;
                              }//if
                             */

                              CycleSlip[PRN].WN_Cycle_Slip(ph1[i],  //PHL1[PRN],
                                                           ph2[i],  //PHL2[PRN],
                                                           p1[i],   //p1[PRN],
                                                           p2[i],   //P2[PRN],
                                                           epoca_ini_cs);
                             /* log << "== step3: " << endl;
                              log << "wlbias: " << CycleSlip[PRN].Get_wlbias() << endl;
                              log << "sig_wlbias: " << CycleSlip[PRN].Get_sig_wlbias() << endl;
                              log << "diff_wl: " << CycleSlip[PRN].Get_diff_wl() << endl;
                              log << "flag: " << CycleSlip[PRN].Get_flag() << endl;
                              nlslip = CycleSlip[PRN].Get_flag();
                              log << "count: " << CycleSlip[PRN].Get_count() << endl;


                              if( wlslip == FALSE  &&  nlslip == TRUE )
                              {
                                   log << "::> PRN:" << setw(2) << setfill('0') << PRN << "  HMS: "
                                       << setw(2) << currentObsEpoch.getEpochTime().GetYMDHMS().hour << ":"
                                       << setw(2) << currentObsEpoch.getEpochTime().GetYMDHMS().min << ":"
                                       << setw(2) << (int) currentObsEpoch.getEpochTime().GetYMDHMS().sec << "  difWL: "
                                       << setw(9) << setprecision(2) << setfill(' ') << CycleSlip[PRN].Get_diff_wl()
                                       << "  NL flag set!" << endl;
                              }//if
                             */
                         }//if
                     }//for

                     File_Cycleslip<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                                    <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                                    <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

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

                }
                //------------------------------------------------------------------------------

                //----------//Tec from pseudorange smothed  by phase----------------------------
                if(read_tec==1)
                {
                     for(int ij=0;ij<NSAT;ij++)
                     {   smoothed_code[ij][0] = 0.0;
                         smoothed_code[ij][1] = 0.0;
                     }

                     int iono = 1; // 1 to use smooth P1 and P2 and 2 for ionosphere free

                     double UTTAG = double(currentObsEpoch.getEpochTime().GetGPSTime().GPSWeek*604800.0 +
                                           currentObsEpoch.getEpochTime().GetGPSTime().secsOfWeek);

                     for (int i = 0; i < NS; i++ )
                     {

                         int PRN = Isat[i];

                         //Start cycle slip
                         if(CycleSlip[PRN].Get_flag_wl()==true || CycleSlip[PRN].Get_flag_nl()==true)
                         {
                              CycleSlip[PRN].Start_CycleSlip();
                              if(PRN>0 && PRN <=32)
                                  FilterObs[PRN].Set_cont_epoch(0);  //it means that smoothing filter will be re-initialized in the next epoch


                         }//if
                         else if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         {
                            if( epoca_ini_cs ||
                                 (CycleSlip[PRN].Get_flag_wl()==true || CycleSlip[PRN].Get_flag_nl()==true) ||
                                 FilterObs[PRN].Get_cont_epoch() == 0)
                             {
                                 //Start Filter code
                                 FilterObs[PRN].Start_Smooth_Code( UTTAG,
                                                                     p1[i],   //p1[PRN],
                                                                     p2[i],   //P2[PRN],
                                                                     ph1[i],  //PHL1[PRN],
                                                                     ph2[i],  //PHL2[PRN],
                                                                     var_ca,
                                                                     var_p2,
                                                                     iono);
                                 //Start cycle slip
                                 if(CycleSlip[PRN].Get_flag_wl()==true || CycleSlip[PRN].Get_flag_nl()==true)
                                     CycleSlip[PRN].Start_CycleSlip();

                             }//if
                             else
                             {
                                 FilterObs[PRN].Filter_Obs( UTTAG,
                                                            p1[i],   //p1[PRN],
                                                            p2[i],   //P2[PRN],
                                                            ph1[i],  //PHL1[PRN],
                                                            ph2[i],  //PHL2[PRN],
                                                            var_ca,
                                                            var_p2,
                                                            var_ph1,
                                                            var_ph2,
                                                            iono );
                             }//else

                             if( FilterObs[PRN].Get_quality() )
                             {
                                 if(iono==2) //ion-free
                                 {
                                     smoothed_code[PRN][0] = FilterObs[PRN].Get_p1_filt();
                                 }
                                 else
                                 {   smoothed_code[PRN][0] = FilterObs[PRN].Get_p1_filt();
                                     smoothed_code[PRN][1] = FilterObs[PRN].Get_p2_filt();
                                 }//else
                             }//if
                             else
                             {
                                 FilterObs[PRN].Start_Smooth_Code( UTTAG,
                                                                     p1[i],   //p1[PRN],
                                                                     p2[i],   //P2[PRN],
                                                                     ph1[i],  //PHL1[PRN],
                                                                     ph2[i],  //PHL2[PRN],
                                                                     var_ca,
                                                                     var_p2,
                                                                     iono);
                                 if(iono==2) //ion-free
                                 {
                                     smoothed_code[PRN][0] = FilterObs[PRN].Get_p1_filt();
                                 }
                                 else //CA and P2
                                 {  smoothed_code[PRN][0] = FilterObs[PRN].Get_p1_filt();
                                    smoothed_code[PRN][1] = FilterObs[PRN].Get_p2_filt();
                                 }
                             }//else
                         }//if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         else
                         {  //Enter here if there are no data available to smothing filter
                            if(PRN>0 && PRN <=32){
                               // FilterObs[PRN].Set_cont_epoch(0);

                                smoothed_code[PRN][0] = Old_smoothed_code[PRN][0];
                                smoothed_code[PRN][1] = Old_smoothed_code[PRN][1];
                            }
                         }//else
                     }//for

                      //Store current code smoothed by phase
                     for(int ij=0;ij<NSAT;ij++)
                     {
                         Old_smoothed_code[ij][0] = smoothed_code[ij][0];
                         Old_smoothed_code[ij][1] = smoothed_code[ij][1];
                     }//for
                }// if(read_tec==1) //Tec from pseudorange smothed by phase
                //----------------------------------------------------------------------------


                for (int i = 0; i < NS; i++ )
                {
                     calc_iono.Set_I1_L1(0.0);
                     calc_iono.Set_I1_L2(0.0);
                     //calc_iono.Set_I1_L5(0.0);

                     calc_iono.Set_I2_L1(0.0);
                     calc_iono.Set_I2_L2(0.0);
                     //calc_iono.Set_I2_L5(0.0);

                     calc_iono.Set_I3_L1(0.0);
                     calc_iono.Set_I3_L2(0.0);
                     //calc_iono.Set_I3_L5(0.0);

                    if(ca[i]!=0.0) //if there are data for CA
                    {
                        prn = Isat[i];

                        int k =0;

                        xs=0.0,ys=0.0,zs=0.0,dxs=0.0,dys=0.0,dzs=0.0;

                        while(k<cont_efe_nav)
                        {
                            if ( (prn == currentPRNBlock[k].getSatellitePRN()) &&
                                 (fabs(tr - currentPRNBlock[k].getToe())< 7200.0) )
                            {
                                //store ephemeris to compute orbit
                                double toe = currentPRNBlock[k].getToe(),       //TOE --> TIME ORIGIN EPHEMERIS:
                                       ra = currentPRNBlock[k].getSqrtA(),      //A   --> SQUARE ROOT SEMI MAJOR AXIS AT TOE:
                                       i0 = currentPRNBlock[k].getIo(),         //RI0 --> INCLINATION AT TOE:
                                       idot = currentPRNBlock[k].getIdot(),     //DRI --> RATE  OF INCLINATION:
                                       dn = currentPRNBlock[k].getDeltan(),     //DN  --> CORRECTION  OF MEAN MOTION:
                                       m0 = currentPRNBlock[k].getMo(),         //CM0 --> MEAN ANOMALY AT TOE:
                                       e = currentPRNBlock[k].getEccen(),       //E   --> ECCENTRICITY:
                                       omega = currentPRNBlock[k].getLilOmega(),//W   --> ARGUMENT OF PERIGEE:
                                       cus = currentPRNBlock[k].getCus(),       //Cus --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE ARGUMENT OF LATITUDE:
                                       cuc = currentPRNBlock[k].getCuc(),       //Cuc --> IBID COSINE:
                                       crs = currentPRNBlock[k].getCrs(),       //Crs --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE ORBITS RADIUS:
                                       crc = currentPRNBlock[k].getCrc(),       //Crc --> IBID COSINE:
                                       cis = currentPRNBlock[k].getCis(),       //Cis --> AMPLITUDE OF THE SINE HARMONIC CORRECTION TERM TO THE INCLINATION ANGLE:
                                       cic = currentPRNBlock[k].getCic(),       //Cic --> IBID COSINE:
                                       omega0 = currentPRNBlock[k].getBigOmega(),       //W0  --> RIGHT ASCENSION AT REFERENCE TIME:
                                       omegadot = currentPRNBlock[k].getBigOmegaDot(),  //WD  --> RATE OF RIGHT ASCENCION:
                                       dt = 0.0;

/* ---- this method to compute tt may be affected by pseudorange multipath
   ---- and it does not include the Sagnac effect!

                                double tau = ca[i]/c;

                                double tt = tr - tau;   //~ transmission time

                                //Satellite coordinates and velocities
                                sat_pos_vel(tt, toe,dt,ra,i0,idot,dn,m0, e,omega,cus,cuc,crs,crc,cis,
                                      cic,omega0,omegadot,xs,ys,zs,dxs,dys,dzs);
*/

//------
                               //Get initial approx Satellite coordinates and velocities (at tr)
                               sat_pos_vel( tr, toe,dt,ra,i0,idot,dn,m0, e,omega,cus,cuc,crs,crc,cis,
                                            cic,omega0,omegadot,xs,ys,zs,dxs,dys,dzs);


                               // Use Gerry's combined light loop and Sagnac code from Mergedb
                               //omega_e = 7.2921158553e-05;
                               //we is a const variable defined in the constant.h file
                               dx = xs - xest;
                               dy = ys - yest;
                               dz = zs - zest;
                               dsv_sq = ( (dx*dx) + (dy*dy) + (dz*dz) );
                               sv_range = sqrt( dsv_sq );

                               for( k = 0; k < 4; k++ )
                               {
                                   propagationDelay = (-1.0 * sv_range)/c;
	                               dtom = propagationDelay * we;       //tau*we   //we is a const variable defined in the constant.h file
                                   usx = xs + dxs*propagationDelay;
                                   usy = ys + dys*propagationDelay;
                                   svx = usx*cos(dtom) - usy*sin(dtom);
                                   svy = usy*cos(dtom) + usx*sin(dtom);

                                   svz = zs + dzs*propagationDelay;

                                   dx = svx - xest;
                                   dy = svy - yest;
                                   dz = svz - zest;
                                   dsv_sq = ( (dx*dx) + (dy*dy) + (dz*dz) );
                                   sv_range = sqrt( dsv_sq );
                                }

                                //double tt= tr + propagationDelay;
                                //coordinates at transmission time and corrected from sagnac effects
                                xs = svx;
                                ys = svy;
                                zs = svz;

                              /*  cout << "PRN, LL_SVxyz: " << prn << endl
                                     << setw(30) << setprecision(5) << xs << endl
                                     << setw(30) << setprecision(5) << ys << endl
                                     << setw(30) << setprecision(5) << zs << endl;

                                getchar();
                                */

//--------
                                //Elevation angle and azimuth of satellite
                                Sat_Angle( lat,lamb,xs,ys,zs,
                                           xest,yest,zest,
                                           EleAngSat[prn],AzAngSat[prn]);

                                double Ang_Dec = EleAngSat[prn]*180.0/PI;

                                if(Ang_Dec > Mask_Ele)  //if elevation angle greater than elevation mask
                                {
                                   //One possibility is to use TGD from broadcast ephemeris
                                   //double TGD = currentPRNBlock[k].getTgd();
                                   //double Bp1p2 = c*(TGD/(-1.0*calc_tec.Get_m2()));

                                   switch(read_tec){
                                       case 0: //TEC from raw pseudorange
                                               //One possibility to correct from DCBs is to use TGD from broadcast ephemeris
                                               //double TGD = currentPRNBlock[k].getTgd();
                                               //double Bp1p2 = c*(TGD/(-1.0*calc_tec.Get_m2()));
                                               STEC = 0.0;
                                               if(p1[i]!=0.0 && p2[i]!=0.0)
                                               {
                                                  //Compute slant TEC from Pseudorange
                                                  STEC = Calc_Tec.Calc_STEC_PR(p1[i],p2[i],prn);
                                               }

                                               if(STEC<0){
                                                  //Sometimes, mainly for low ionospheric activities period, we get negative TEC from pseudorange
                                                  //Then, warning user
                            	                  log<<"Warning: "<<"Satellite "<<setw(3)<<prn<<" - Negative STEC ( "<<setprecision(4)<<setw(10)<<STEC<<" ) from Pseudorange in the epoch : "<<cont_epoch<<"..."<<endl;
                            	                  log<<"Difference between CA and P2: "<<p1[i]-p2[i]<<endl;
                                               }//if

                                               if( STEC != 0.0) //Compute Second and third order ionospheric effects
                                                  calc_iono.Calc_Iono(STEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                      xest,yest,zest,lat,lamb,h,
                                                                      Geomag_Type,Name_IGRF_Model);

                                               break;
                                       case 1: //TEC from pseudorange smoothed by phase
                                               STEC = 0.0;

                                               if(smoothed_code[prn][0]!=0.0 && smoothed_code[prn][1]!=0.0){
                                                     //Wait almost 2 epochs until the smoothing filter get good values for smoothed pseudorange
                                                     //when there is cycle slip sometimes we will see gaps in the higher order iono time series due to STEC missig
                                                      if(FilterObs[prn].Get_cont_epoch() > 2){
                                                          //Compute slant TEC from Smoothed Pseudorange
                                                          STEC = (Calc_Tec.Calc_STEC_PR(smoothed_code[prn][0],smoothed_code[prn][1],prn));
                                                      }

                                                      if(STEC<0){
                                                        //Sometimes, mainly for low ionospheric activies periods, we got negative TEC from pseudorange
                                                        //warning user
                                	                log<<"Warning: "<<"Satellite "<<setw(3)<<prn<<" - Negative STEC ( "<<setprecision(4)<<setw(10)<<STEC<<" ) from Pseudorange in the epoch : "<<cont_epoch<<"..."<<endl;
                                	                log<<"Difference between CA and P2: "<<p1[i]-p2[i]<<endl;
                                                     }//if

                                                     if(STEC != 0.0) //Compute Second and third order ionospheric effects
                                                        calc_iono.Calc_Iono(STEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                            xest,yest,zest,lat,lamb,h,
                                                                            Geomag_Type,Name_IGRF_Model);
                                               }//if

                                               break;
                                       case 2: //TEC from GIM  (now get a VTEC for each IPP, each satellite)

                                              VTEC = 0.0;

                                              calc_iono.azelIPP( xest,yest,zest,xs,ys,zs,6371000.0,450000.0,
                                                                 &sh_az,&sh_elv,&sh_dist,&sh_pplat,&sh_pplon,&sh_zp);

                                             // calc_iono.Pierce_Point_Coord(xest,yest,zest,xs,ys,zs,lat,lamb,h);

                                              //Marques, H. A.
                                              //lat. and long. of IPP in degrees
                                               lat_IPP_deg = sh_pplat*180.0/PI;
                                               lon_IPP_deg = sh_pplon*180.0/PI;

                                              //Interpolate TEC for MJD+fracOfDay
                                              //**SH** if(! Calc_Tec.CalcTecEpoch(MJD,fracOfDay,lat_deg,long_deg) )
                                              //if(!Calc_Tec.CalcTecEpoch(MJD,fracOfDay,lat_deg,long_deg) )
                                              if(!Calc_Tec.CalcTecEpoch(MJD,fracOfDay,lat_IPP_deg,lon_IPP_deg) )
                                              {
                                                   log<<"Without TEC information from GIM "
                                                      <<"in the epoch: "<<calc_iono.Get_hour()<<" "
                                                      <<calc_iono.Get_min()<<" "
                                                      <<calc_iono.Get_sec()<<endl;

                                                   return 0;

                                               }//if
                                               else
                                               {
                                                   VTEC = Calc_Tec.GetVTEC();

                                                  // STEC = ( (VTEC)/cos(calc_iono.zl) )*1.0e16;

                                               }//else

                                               //Compute Second and third order ionospheric effects
                                               calc_iono.Calc_Iono(VTEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                   xest,yest,zest,lat,lamb,h,
                                                                  Geomag_Type,Name_IGRF_Model);
                                               break;
                                         case 3: //STEC from Phase Levelled by Code
                                                STEC=0.0;
                                                if(p1[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0){

                                                    if( CycleSlip[prn].Get_flag_wl()==true || CycleSlip[prn].Get_flag_nl()==true ) //flag for wide lane or flag for narrow lane
                                                    {  Cycle_Slip_Flag  = true;
                                                       CycleSlip[prn].Start_CycleSlip();
                                                    }
                                                    else
                                                      Cycle_Slip_Flag=false;

                                                  /*   Calc_Tec.Mean_Phase_Code(ph1[i],ph2[i],p1[i],p2[i],
                                                                    Sat_Average_Phase_Code[prn].Mean_Diff_Li_Pi,
                                                                     Sat_Average_Phase_Code[prn].Last_Diff_Li_Pi,
                                                                     Sat_Average_Phase_Code[prn].Count_Mean,
                                                                     cont_epoch,Cycle_Slip_Flag,prn);
                                                    */
                                                   double ModJDate = MJD + fracOfDay;
                                                   double PhBias = Get_Phase_Bias(prn,Phase_Bias,ModJDate);
                                                   if(PhBias != 0.0)
                                                      STEC = Calc_Tec.Calc_STEC_PH_PR( ph1[i],ph2[i],p1[i],p2[i],PhBias,prn);

                                                   //Get_Phase_Bias(Phase_Bias, PhBias);

                                                   //  if(Sat_Average_Phase_Code[prn].Count_Mean > 2)
                                                   //     STEC = Calc_Tec.Calc_STEC_PH_PR( ph1[i],ph2[i],p1[i],p2[i],
                                                   //                                      Sat_Average_Phase_Code[prn].Mean_Diff_Li_Pi,prn);
                                                 }//if

                                                 if(STEC<0){
                                                 //Sometimes, mainly for low ionospheric activies periods, we got negative TEC from pseudorange
                                                 //warning user
                                	             log<<"Warning: "<<"Satellite "<<setw(3)<<prn<<" - Negative STEC ( "<<setprecision(4)<<setw(10)<<STEC<<" ) from Pseudorange in the epoch : "<<cont_epoch<<"..."<<endl;
                                	             //log<<"Difference between CA and P2: "<<p1[i]-p2[i]<<endl;
                                                 }//if

                                                 if(STEC != 0.0) //Compute Second and third order ionospheric effects
                                                   calc_iono.Calc_Iono( STEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                        xest,yest,zest,lat,lamb,h,
                                                                        Geomag_Type,Name_IGRF_Model);

                                                 break;

                                         default: break;
                                    };//switch

                                    if(save_file) //if true, store data to be saved posteriori
                                    {
                                        i2m_l1[prn] = calc_iono.Get_I2_L1(); //Second order - L1
                                        i2m_l2[prn] = calc_iono.Get_I2_L2(); //Second order - L2
                                                      //calc_iono.Get_I2_L5(); //Second order - L5
                                        i3m_l1[prn] = calc_iono.Get_I3_L1(); //Third order - L1
                                        i3m_l2[prn] = calc_iono.Get_I3_L2(); //Third order - L2
                                                      //calc_iono.Get_I3_L5(); //Third order - L5

                                        if(read_tec==2)
                                        {   //Store Slant TEC from GIM to be printed posteriori
                                            TECM[prn] = static_cast< float > (VTEC);
                                        }
                                        else //if(read_tec==0 || read_tec ==1 ||read_tec==3)
			                            {   //Store Slant TEC computed from GNSS data to be printed posteriori
                                            TECM[prn] = static_cast< float > (STEC)*cos(calc_iono.zl);
		                                }
                                    }//if save_file
                                }//if(Ang_Dec > Mask_Ele)  //if elevation angle is greater than elevation mask
                                else

                                {
                                    //MHA: When we use this line, the satellite still are below the elevation mask
                                    //and the smoothed pseudorange filter will start computations
                                    //only after satellite reach the input elevation mask
                                    //commenting this line means that smoothed filtering starts as soon as satellite data are available

                                     FilterObs[prn].Set_cont_epoch(0);
                                    // Sat_Average_Phase_Code[prn].Count_Mean=0;
                                     CycleSlip[prn].Start_CycleSlip();

                                }//else

                                break;

                            }//if ( (prn == currentPRNBlock[k].getSatellitePRN()) && (fabs(tr - currentPRNBlock[k].getToe())< 7200.0) )

                            k++;
                        }//while(k<cont_efe_nav)
                    }//if (ca[i]!=0.0)

                    if(xs == 0.0 && ys == 0.0 && zs == 0.0 && dxs == 0.0 && dys == 0.0 && dzs == 0.0)
                    {   //warning user
                        log<<"Without ephemeris for the satellite: "<<currentObsEpoch.getSatListElement(i).satNum
                           <<" in the ephoc "<<calc_iono.Get_hour()<<" "
                           <<calc_iono.Get_min()<<" "
                           <<calc_iono.Get_sec()<<endl;

                    }//if

                    //Correct GPS observables from 2nd and 3rd order ionospheric effects - L2C and L5 not being used
                    CorrectIono( numobstype,currentObsEpoch,tempSatObsAtEpoch,
                                 calc_iono.Get_I1_L1(),calc_iono.Get_I1_L2(),calc_iono.Get_I1_L5(), //I1 not being used
                                 calc_iono.Get_I2_L1(),calc_iono.Get_I2_L2(),calc_iono.Get_I2_L5(),
                                 calc_iono.Get_I3_L1(),calc_iono.Get_I3_L2(),calc_iono.Get_I3_L5(),
                                 ph1,ph2,ph5,       //phase
                                 ca,p1,p2,c2,c5,    //code
                                 Pos_GPS_CurrObs[i],i);


                }//for (int i = 0; i < NS; i++ )

                //Updating currentObsEpoch object with corrected GPS observable
                for (int i = 0; i < NS; i++)
                 currentObsEpoch.setSatListElement(tempSatObsAtEpoch[i],num,Pos_GPS_CurrObs[i]);


                //Write corrected observations in the new RINEX obs file
                myobs.writeEpoch(mynewobs.outputStream, currentObsEpoch);

                if(save_file)
                {


                    //Print time
                    //tecfile<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                    //                <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                    //                <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                    tecfile<<setw(2)<<setfill('0')<<calc_iono.Get_hour()+
                                                    calc_iono.Get_min()/60.0+
                                                    calc_iono.Get_sec()/3600.0;

                   // fstionoI2l1<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                   //            <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                   //            <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                   fstionoI2l1<<setw(2)<<setfill('0')<<calc_iono.Get_hour()+
                                                    calc_iono.Get_min()/60.0+
                                                    calc_iono.Get_sec()/3600.0;

                    fstionoI2l2<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                               <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                               <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                    fstionoI3l1<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                               <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                               <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                    fstionoI3l2<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                               <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                               <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                    File_sateleangle<<setw(2)<<setfill('0')<<calc_iono.Get_hour()<<":"
                                    <<setw(2)<<setfill('0')<<calc_iono.Get_min()<<":"
                                    <<setw(4)<<setfill('0')<<setprecision(1)<<calc_iono.Get_sec();

                    tecfile<<setfill(' ');
                    fstionoI2l1<<setfill(' ');
                    fstionoI2l2<<setfill(' ');
                    fstionoI3l1<<setfill(' ');
                    fstionoI3l2<<setfill(' ');
                    File_sateleangle<<setfill(' ');

                    for(int j=1;j<=32;j++)
                    {
                        if(read_tec==2)
                           STEC = TECM[j];
                        else //if(read_tec==0 || read_tec ==1 || read_tec==3)
                           STEC = TECM[j]*1.0e-16;

                        if(STEC!=0.0)
                           tecfile<<setw(13)<<setprecision(6)<<STEC;
                        else {
                           tecfile<<setw(13)<<setprecision(6)<<9999999;
                           //tecfile<<setw(13)<<"NaN";
                        }

                        if(i2m_l1[j] != 0)
                           fstionoI2l1<<setw(13)<<setprecision(6)<<i2m_l1[j];
                        else{
                           fstionoI2l1<<setw(13)<<setprecision(6)<<9999999;
                           //fstionoI2l1<<setw(13)<<"NaN";
                        }

                        if(i2m_l2[j] != 0)
                          fstionoI2l2<<setw(13)<<setprecision(6)<<i2m_l2[j];
                        else
                           fstionoI2l2<<setw(13)<<setprecision(6)<<9999999;

                        if(i3m_l1[j] != 0)
                          fstionoI3l1<<setw(13)<<setprecision(6)<<i3m_l1[j];
                        else
                          fstionoI3l1<<setw(13)<<setprecision(6)<<9999999;

                        if(i3m_l2[j] != 0)
                          fstionoI3l2<<setw(13)<<setprecision(6)<<i3m_l2[j];
                        else
                          fstionoI3l2<<setw(13)<<setprecision(6)<<9999999;


                        if (EleAngSat[j] != 0.0)
                          File_sateleangle<<setw(13)<<setprecision(6)<< EleAngSat[j]*180.0/PI;
                        else
                          File_sateleangle<<setw(13)<<setprecision(6)<<9999999;

                    }//for

                    tecfile<<endl;
                    fstionoI2l1<<endl;
                    fstionoI2l2<<endl;
                    fstionoI3l1<<endl;
                    fstionoI3l2<<endl;
                    File_sateleangle<<endl;
                }//if save_file

                cont_epoch++; //counter of epochs

                currentObsEpoch.initializeData();

           }//while ( myobs.readEpoch( currentObsEpoch ) != 0 )
        }//try
        catch( RinexReadingException &readingExcep )
        {
            cout << " RinexReadingException is: " << endl
                 << readingExcep.getMessage() << endl;

            log << " RinexReadingException is: " << endl
                <<readingExcep.getMessage() << endl;
        }//catch

        log << "*** Error Message in the RINEX file:" << endl;
        log << myobs.getErrorMessages() << endl;
        log << "*** Messages with WARNINGS in the RINEX file:" << endl;
        log << myobs.getWarningMessages() << endl;
    }//if theObsFileType[0] == 'O'

  cout<<endl<<"Process Finished Normally..."<<endl<<endl;
  log<<endl<<"Process Finished Normally..."<<endl<<endl;

  cout<<endl<<"Please, check log files to see results and warnings..."<<endl<<endl;

  log.close();
  tecfile.close();
  fstionoI2l1.close();
  fstionoI2l2.close();
  fstionoI3l1.close();
  fstionoI3l2.close();
  File_sateleangle.close();
  File_Cycleslip.close();

  return 0;

}
