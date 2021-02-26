
/*-------------------------------------------------------------------
    Purpose:  Program to read a RINEX file and apply the second and
              third order ionospheric effects corrections
    -------------------------------------------------------------------
    Authors: Haroldo Antonio Marques
             haroldoh2o@gmail.com
             Programa de Pos-Graduacao em Ciencias Cartograficas
             FCT/UNESP - Presidente Prudente - SP, Brazil
             FAPESP PROCESS: 05/03522-1

    Advisor: Joao Francisco Galera Monico
             galera@fct.unesp.br
             Departamento de Cartografia
             FCT/UNESP - Presidente Prudente - SP, Brazil

    Date: July of 2010
    -------------------------------------------------------------------
    Observation: Use C++ classes to read and save RINEX file:
                 http://www.ngs.noaa.gov/gps-toolbox/rinex.htm
                 The authors would like to thanks CAPES, FAPESP by financial support
    -------------------------------------------------------------------*/

#include "rinex_ho.h"

int main(int argc, char* argv[]) {

    #include "Include/Variable.inc"

    system("clear");  //for linux
    //clrscr();           //for windows

    //check-out input parameters
    if(argc == 1)
      InpFileName = "rinex_ha.inp";
    else
       InpFileName = argv[1];

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
        system("pause");
        exit (-1);
    }//if

    string temp = outfile+"_tec.txt";
    ofstream tecfile(temp.c_str());
    if( tecfile.fail() )
    {
        cout << "It was not possible to open the file: tec.txt!" << endl;
        system("pause");
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
        system("pause");
        exit (-1);
    }//if

    temp = outfile+"_I2L2.txt";
    ofstream fstionoI2l2(temp.c_str());
    if( fstionoI2l2.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        system("pause");
        exit (-1);
    }//if

    temp = outfile+"_I3L1.txt";
    ofstream fstionoI3l1(temp.c_str());
    if( fstionoI3l1.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        system("pause");
        exit (-1);
    }//if

    temp = outfile+"_I3L2.txt";
    ofstream fstionoI3l2(temp.c_str());
    if( fstionoI3l2.fail())
    {
        cout << "It was not possible to open the file: " << temp<< endl;
        system("pause");
        exit (-1);
    }//if

    //---------------Configure output files-------------------------------------
    fstionoI2l1.setf(ios::fixed);       fstionoI2l2.setf(ios::fixed);
    fstionoI2l1.setf(ios::showpoint);   fstionoI2l2.setf(ios::showpoint);
    fstionoI2l1.setf(ios::right);       fstionoI2l2.setf(ios::right);

    fstionoI3l1.setf(ios::fixed);       fstionoI3l2.setf(ios::fixed);
    fstionoI3l1.setf(ios::showpoint);   fstionoI3l2.setf(ios::showpoint);
    fstionoI3l1.setf(ios::right);       fstionoI3l2.setf(ios::right);


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

    if( read_tec == 0 ||  read_tec == 1) //if not to use GIM, read DCBs
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
    {   system("pause");
        exit(-1);
    }//if

    cout<<"Observation file: "<<filenameobs<<endl;
    log<<"Observation file: "<<filenameobs<<endl;

    if(read_tec==2)
    { cout<<"GIM File: "<<GIM_File<<endl;
      log<<"GIM File: "<<GIM_File<<endl;
    }

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
            }
                                //Decimal degrees
            double lat_deg = lat*180.0/calc_iono.Get_PI();
            double long_deg = lamb*180.0/calc_iono.Get_PI();

            while( myobs.readEpoch( currentObsEpoch ) != 0 )
            {
                //print number of epochs in the screen
                //For linux
                printf("\033[20;4f");

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
                }//if

                //Reception time in seconds of week
                tr =  segweek;

                //start obs arrays
                for(int ij=0;ij<NSAT;ij++) ph1[ij] = ph2[ij] = ca[ij] = p1[ij] = p2[ij] = 0.0;

                for(int ij=0;ij<NSAT;ij++)
                {   smoothed_code[ij][0] = 0.0;
                    smoothed_code[ij][1] = 0.0;
                }

                int numobstype = myobs.getNumObsTypes();

                for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )
                {  //Store observations in arrays
                    StoreObsVector(numobstype,currentObsEpoch,ph1,ph2,ca,p1,p2,i);
                }

                //----------//Tec from pseudorange smothed  by phase----------------------------
                if(read_tec==1)
                {
                     //check-out for cycle slips:
                     for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )
                     {   int PRN = currentObsEpoch.getSatListElement(i).satNum;
                         if(ca[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         {
                             CycleSlip[PRN].Wide_Cycle_Slip( ph1[i],  //PHL1[PRN],
                                                             ph2[i],  //PHL2[PRN],
                                                             ca[i],   //CA[PRN],
                                                             p2[i],   //P2[PRN],
                                                             epoca_ini_cs);
                         }//if
                       /*  else
                         {
                            CycleSlip[PRN].Start_CycleSlip();
                         }
                         */
                     }//for

                     double UTTAG = currentObsEpoch.getEpochTime().GetGPSTime().GPSWeek*604800.0 +
                                    currentObsEpoch.getEpochTime().GetGPSTime().secsOfWeek;

                     int iono =1;

                     for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )
                     {
                         int PRN = currentObsEpoch.getSatListElement(i).satNum;

                         if(ca[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         {
                             if( epoca_ini_cs ||
                                 CycleSlip[PRN].Get_flag()==true ||
                                 FilterObs[PRN].Get_cont_epoch() == 0)
                             {   //Start Filter code

                                 FilterObs[PRN].Start_Smooth_Code( UTTAG,
                                                                     ca[i],   //CA[PRN],
                                                                     p2[i],   //P2[PRN],
                                                                     ph1[i],  //PHL1[PRN],
                                                                     ph2[i],  //PHL2[PRN],
                                                                     var_ca,
                                                                     var_p2,
                                                                     iono);
                                 //Start cycle slip
                                 if(CycleSlip[PRN].Get_flag()==true)
                                     CycleSlip[PRN].Start_CycleSlip();

                             }//if
                             else
                             {   //TECM[PRN] = 0.0;
                                 FilterObs[PRN].Filter_Obs( UTTAG,
                                                            ca[i],   //CA[PRN],
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
                                                                     ca[i],   //CA[PRN],
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
                         }//if(ca[i] != 0.0 && p2[i] != 0.0 && ph1[i] !=0.0 && ph2[i] != 0.0)
                         else
                         {

                             FilterObs[PRN].Set_cont_epoch(0);

                             smoothed_code[PRN][0] = ca[i];
                             smoothed_code[PRN][1] = p2[i];

                             //TECM[PRN] = 0.0;

                         }//else

                     }//for
                }// if(read_tec==1) //Tec from pseudorange smothed by phase
                //----------------------------------------------------------------------------

                for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )
                {

                    int PRN = currentObsEpoch.getSatListElement(i).satNum;

                    if(ca[i]!=0.0 && p2[i]!=0.0) //if there are data for CA and P2
                    {
                        prn = currentObsEpoch.getSatListElement(i).satNum;

                        int k =0;

                        calc_iono.Set_I1_L1(0.0);
                        calc_iono.Set_I1_L2(0.0);
                        calc_iono.Set_I2_L1(0.0);
                        calc_iono.Set_I2_L2(0.0);
                        calc_iono.Set_I3_L1(0.0);
                        calc_iono.Set_I3_L2(0.0);

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

                                double tau = ca[i]/c;

                                double tt = tr - tau;

                                //Satellite coordinates and velocities
                                sat_pos_vel(tt, toe,dt,ra,i0,idot,dn,m0, e,omega,cus,cuc,crs,crc,cis,
	                                     cic,omega0,omegadot,xs,ys,zs,dxs,dys,dzs);

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
                                               //double TGD = currentPRNBlock[k].getTgd();
                                               //double Bp1p2 = c*(TGD/(-1.0*calc_tec.Get_m2()));
                                               TEC = 0.0;
                                               //Compute slant TEC from Pseudorange
                                               TEC = Calc_Tec.Calc_STEC_PR(ca[i],p2[i],prn);

                                               if(TEC<0){
                                                  //Sometimes, mainly for low ionospheric activities period, we got negative TEC from pseudorange
                                                  //Then, warning user
                            	                  log<<"Warning: "<<"Satellite "<<setw(3)<<prn<<" - Negative STEC ( "<<setprecision(4)<<setw(10)<<TEC<<" ) from Pseudorange in the epoch : "<<cont_epoch<<"..."<<endl;
                            	                  log<<"Difference between CA and P2: "<<ca[i]-p2[i]<<endl;
                                               }//if

                                               //Compute Second and third order ionospheric effects
                                               calc_iono.Calc_Iono(TEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                   xest,yest,zest,lat,lamb,h,
                                                                   Geomag_Type,Name_IGRF_Model);
                                               break;
                                       case 1: //TEC from pseudorange smoothed by phase
                                               TEC = 0.0;
                                               //Wait almost 2 epochs until the filter get good values for smoothed pseudorange
                                               if(FilterObs[PRN].Get_cont_epoch() > 2){
                                                  //Compute slant TEC from smoothed Smoothed Pseudorange
                                                  TEC = (Calc_Tec.Calc_STEC_PR(smoothed_code[PRN][0],smoothed_code[PRN][1],prn));
                                               }

                                               if(TEC<0){
                                                 //Sometimes, mainly for low ionospheric activies periods, we got negative TEC from pseudorange
                                                 //warning user
                                	             log<<"Warning: "<<"Satellite "<<setw(3)<<prn<<" - Negative STEC ( "<<setprecision(4)<<setw(10)<<TEC<<" ) from Pseudorange in the epoch : "<<cont_epoch<<"..."<<endl;
                                	             log<<"Difference between CA and P2: "<<ca[i]-p2[i]<<endl;
                                               }//if

                                               if(TEC != 0.0) //Compute Second and third order ionospheric effects
                                                  calc_iono.Calc_Iono(TEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                      xest,yest,zest,lat,lamb,h,
                                                                      Geomag_Type,Name_IGRF_Model);

                                               break;
                                       case 2: //TEC from GIM

                                               //Compute Second and third order ionospheric effects
                                               calc_iono.Calc_Iono(VTEC,read_tec,MJD,fracOfDay,xs,ys,zs,
                                                                   xest,yest,zest,lat,lamb,h,
                                                                  Geomag_Type,Name_IGRF_Model);
                                               break;
                                    };//switch
                                }//if(Ang_Dec > Mask_Ele)  //if elevation angle greater than elevation mask

                                if(save_file) //if true, store data to be saved posteriori
                                {
                                    i2m_l1[prn] = calc_iono.Get_I2_L1(); //Second order - L1
                                    i2m_l2[prn] = calc_iono.Get_I2_L2(); //Second order - L2
                                    i3m_l1[prn] = calc_iono.Get_I3_L1(); //Third order - L1
                                    i3m_l2[prn] = calc_iono.Get_I3_L2(); //Third order - L2

                                    if(read_tec==0 || read_tec ==1)
			                        {   //Store Slant TEC from pseudorange to be printed posteriori
                                        TECM[prn] = static_cast< float > (TEC);
		                            }
                                    else if(read_tec==2)
                                    {   //Store Slant TEC from GIM to be printed posteriori
                                        TECM[prn] = static_cast< float > (VTEC)/cos(calc_iono.zl);
                                    }
                                }//if save_file

                                break;

                            }//if ( (prn == currentPRNBlock[k].getSatellitePRN()) && (fabs(tr - currentPRNBlock[k].getToe())< 7200.0) )

                            k++;
                        }//while(k<cont_efe_nav)
                    }//if (ca[i]!=0.0&&p2[i]!=0.0)

                    if(xs == 0.0 && ys == 0.0 && zs == 0.0 && dxs == 0.0 && dys == 0.0 && dzs == 0.0)
                    {   //warning user
                        log<<"Without ephemeris for the satellite: "<<currentObsEpoch.getSatListElement(i).satNum
                           <<" in the ephoc "<<calc_iono.Get_hour()<<" "
                           <<calc_iono.Get_min()<<" "
                           <<calc_iono.Get_sec()<<endl;

                    }//if

                    //Correct GPS observables from 2nd and 3rd order ionospheric effects
                    CorrectIono( numobstype,currentObsEpoch,tempSatObsAtEpoch,
                                 calc_iono.Get_I1_L1(),calc_iono.Get_I1_L2(),
                                 calc_iono.Get_I2_L1(),calc_iono.Get_I2_L2(),
                                 calc_iono.Get_I3_L1(),calc_iono.Get_I3_L2(),
                                 ph1,ph2,ca,p1,p2,i);

                }//for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )


                for (int i = 0; i < (unsigned short) currentObsEpoch.getNumSat(); i++ )
                    currentObsEpoch.setSatListElement(tempSatObsAtEpoch[i],num,i);

                //Write corrected observations in the new RINEX obs file
                myobs.writeEpoch(mynewobs.outputStream, currentObsEpoch);

                if(save_file)
                {
                    if(cont_epoch == 0)
                    {   //In the first epoch print the output file headers
                        string nome = filenameobs.substr(0,8);

                        tecfile<<setw(8)<<"h:m:s";
                        fstionoI2l1<<setw(8)<<"h:m:s"; 
                        fstionoI2l2<<setw(8)<<"h:m:s";
                        fstionoI3l1<<setw(8)<<"h:m:s"; 
                        fstionoI3l2<<setw(8)<<"h:m:s";

                        for(int j=1;j<=32;j++)
                        {
                            if(j<10)
                            {   fstionoI2l1<<setw(12)<<"PRN0"<<j;  
                                fstionoI2l2<<setw(12)<<"PRN0"<<j;
                                fstionoI3l1<<setw(12)<<"PRN0"<<j;  
                                fstionoI3l2<<setw(12)<<"PRN0"<<j;
                                tecfile<<setw(12)<<"PRN0"<<j;
                            }//if
                            else
                            {   fstionoI2l1<<setw(11)<<"PRN"<<j;  
                                fstionoI2l2<<setw(11)<<"PRN"<<j;
                                fstionoI3l1<<setw(11)<<"PRN"<<j;  
                                fstionoI3l2<<setw(11)<<"PRN"<<j;
                                tecfile<<setw(11)<<"PRN"<<j;
                            }//else
                        }//for

                        tecfile<<endl;
                        fstionoI2l1<<endl;  
                        fstionoI2l2<<endl;
                        fstionoI3l1<<endl;  
                        fstionoI3l2<<endl;
                    } //if(cont_epoch == 0)

                    //Print time
                    tecfile<<calc_iono.Get_hour()<<":"<<calc_iono.Get_min()<<":"<<setprecision(1)<<calc_iono.Get_sec();
                    fstionoI2l1<<calc_iono.Get_hour()<<":"<<calc_iono.Get_min()<<":"<<setprecision(1)<<calc_iono.Get_sec();
                    fstionoI2l2<<calc_iono.Get_hour()<<":"<<calc_iono.Get_min()<<":"<<setprecision(1)<<calc_iono.Get_sec();
                    fstionoI3l1<<calc_iono.Get_hour()<<":"<<calc_iono.Get_min()<<":"<<setprecision(1)<<calc_iono.Get_sec();
                    fstionoI3l2<<calc_iono.Get_hour()<<":"<<calc_iono.Get_min()<<":"<<setprecision(1)<<calc_iono.Get_sec();


                    for(int j=1;j<=32;j++)
                    {   if(read_tec==0 || read_tec ==1) TEC = TECM[j]*1.0e-16;
                        else if(read_tec==2) TEC = TECM[j];

                        if(TEC!=0.0)
                           tecfile<<setw(13)<<setprecision(6)<<TEC;
                        else
                           tecfile<<setw(13)<<setprecision(6)<<9999999;

                        if(i2m_l1[j] != 0)
                           fstionoI2l1<<setw(13)<<setprecision(6)<<i2m_l1[j];
                        else
                           fstionoI2l1<<setw(13)<<setprecision(6)<<9999999;

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
                    }//for

                    tecfile<<endl;
                    fstionoI2l1<<endl;  
                    fstionoI2l2<<endl;
                    fstionoI3l1<<endl;  
                    fstionoI3l2<<endl;           
                }//if save_file

                cont_epoch++; //counter of epochs
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

  return 0;

}
