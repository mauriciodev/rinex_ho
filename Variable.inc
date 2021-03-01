  //MHA: may of 2013 - Change TEC variable by STEC

 //Variable for RINEX_HO

 double segweek,tr,
          ph1[MAXPRNID]={0.0},    //Phase - L1
          ph2[MAXPRNID]={0.0},    //Phase - L2
          ca[MAXPRNID]={0.0},     //code - c1
          p1[MAXPRNID]={0.0},     //code - p1
          p2[MAXPRNID]={0.0},     //code - p2
          c2[MAXPRNID]={0.0},     //code - L2C
          c5[MAXPRNID]={0.0},     //code C5 - L5
          ph5[MAXPRNID]={0.0},    //Phase - L5
         // dopp1[MAXPRNID]={0.0},
         // dopp2[MAXPRNID]={0.0},
         // snr1[MAXPRNID]={0.0},
         // snr2[MAXPRNID]={0.0},
          smoothed_code[MAXPRNID][2]={0.0},
          Old_smoothed_code[MAXPRNID][2]={0.0},
          STEC(0.0);

   double X0(0.0),Y0(0.0),Z0(0.0),xest(0.0),yest(0.0),zest(0.0),lat(0.0),lamb(0.0),h(0.0), //Station coordinates (read in the Header of RINEX file)
          xs(0.0),ys(0.0),zs(0.0),dxs(0.0),dys(0.0),dzs(0.0),  //satellite coordinates and velocities
          sigca(0.0),sigp2(0.0),sigph1(0.0),sigph2(0.0), //standard deviations (used to do covariance propagation)
          h_ion(0.0),  //Ionospheric layer height
          R_e(0.0),    //Earth equatorial axis
          B_eq(0.0),   //geomagnetic induction magnitude in the Equator (Tesla)
          Ne_max(0.0); //Maximun density of eletrons

    double EleAngSat[MAXPRNID], //Elevation angles
           AzAngSat[MAXPRNID];  //Azimuth angles

    double Mask_Ele = 10.0;     //Elevation mask - default = 10.0 degrees

    double VTEC(0.0), //Vertical TEC
           dcbr(0.0); //Receiver hardware delay

    int cont_efe_nav(0),cont_epoch(0),prn(0);  //counters
    bool epoca_ini_cs = false;

    string filenameobs,filenamenav,newobsfile,outfile;  //name of files

    string theObsFileType,  //type of observation
           theNavFileType;

    int save_file(0);

    char GIM_File[300]="/0"; //GIM file name

    char arqu_p1c1_sat[300],arqu_p1p2_sat[300]; //DCB files names

    char Name_IGRF_Model[300]; //IGRF model file name

    int read_tec(0);
    int Geomag_Type = 2; //0 = dipolar model; 1 = geomg coord. transformation from PIM; 2 = IGRF model

    //Class_GIM R1; //object of Class Read_Tec

    Class_Cycle_Slip CycleSlip[NSAT];                   //object of class CycleSlip
    FilterCode FilterObs[NSAT];

    string InpFileName; //Input file: Default name is "rinex_ha.inp'

    PRNBlock  currentPRNBlock[500];           //object to store broadcast ephemeris
