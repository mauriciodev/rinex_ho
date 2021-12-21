#include "test.h"



// add necessary includes here



rinex_ho_Tests::rinex_ho_Tests()
{

}

rinex_ho_Tests::~rinex_ho_Tests()
{

}


void rinex_ho_Tests::test_RINEX2() {
    QBENCHMARK {
        NGSrinex::RinexNavFile mynav;
        std::string filenamenav="RECF/recf0101.16n";
        try {
            mynav.setPathFilenameMode(filenamenav,ios_base::in);
        } catch( RinexFileException &openExcep ) {
            cout << "Error opening file: " << filenamenav << endl;
            cout << " RinexFileException is: " << endl <<
                    openExcep.getMessage() << endl;
        }//catch
        try {
            mynav.readHeader();
        } catch( RequiredRecordMissingException &headerExcep ) {
            cout << " RequiredRecordMissingException is: " << endl <<
                    headerExcep.getMessage() << endl;
        }

        //read the PRN Blocks
        PRNBlock currentPRNBlock;
        try
        {   int cont_efe_nav=0; //number of ephemeris block

            while( mynav.readPRNBlock( currentPRNBlock ) != 0 )
            {
                cont_efe_nav++;
            }//while
        }//try
        catch( RinexReadingException &readingExcep )
        {
            cout << " RinexReadingException is: " << endl <<
                    readingExcep.getMessage() << endl;
        }//catch

        if (mynav.getNumberWarnings()>0) {
            qDebug() << "*** Messages WARNINGS:" << endl;
            qDebug() << mynav.getWarningMessages().c_str() << endl;
        }
        if (mynav.getNumberErrors()>0) {
            qDebug() << "*** Messages ERRORS:" << endl;
            qDebug() << mynav.getErrorMessages().c_str() << endl;
        }
    }
}

void rinex_ho_Tests::test_RINEX3nav() {
    QBENCHMARK {
        NGSrinex::RinexNavFile mynav;
        std::string filenamenav="RECF/BRDC00WRD_S_20210680000_01D_MN.rnx";
        try {
            mynav.setPathFilenameMode(filenamenav,ios_base::in);
        } catch( RinexFileException &openExcep ) {
            cout << "Error opening file: " << filenamenav << endl;
            cout << " RinexFileException is: " << endl <<
                    openExcep.getMessage() << endl;
        }//catch
        try {
            mynav.readHeader();
        } catch( RequiredRecordMissingException &headerExcep ) {
            cout << " RequiredRecordMissingException is: " << endl <<
                    headerExcep.getMessage() << endl;
        }

        //read the PRN Blocks
        PRNBlock currentPRNBlock;
        try
        {   int cont_efe_nav=0; //number of ephemeris block
            if (mynav.getFormatVersion()>=3.0) {
                while( mynav.readPRNBlockV3( currentPRNBlock ) != 0 )
                {
                    cont_efe_nav++;
                }//while
            } else {
                while( mynav.readPRNBlock( currentPRNBlock ) != 0 )
                {
                    cont_efe_nav++;
                }//while
            }
        }//try
        catch( RinexReadingException &readingExcep )
        {
            cout << " RinexReadingException is: " << endl <<
                    readingExcep.getMessage() << endl;
        }//catch
        if (mynav.getNumberWarnings()>0) {
            qDebug() << "*** Messages WARNINGS:" << endl;
            qDebug() << mynav.getWarningMessages().c_str() << endl;
        }
        if (mynav.getNumberErrors()>0) {
            qDebug() << "*** Messages ERRORS:" << endl;
            qDebug() << mynav.getErrorMessages().c_str() << endl;
        }
    }
}

void rinex_ho_Tests::test_RINEX3obs() {
    std::string filenameobs="RECF/RJNI00BRA_R_20210680000_01D_15S_MO.rnx";
    RinexFile *myFile = new RinexFile();   //object from class RinexFile
    double X0(0.0),Y0(0.0),Z0(0.0),xest(0.0),yest(0.0),zest(0.0),lat(0.0),lamb(0.0),h(0.0); //Station coordinates (read in the Header of RINEX file)
    string theObsFileType;
    try {  //try to open rinex file
            myFile->setPathFilenameMode( filenameobs, ios::in );
            cout<<" - Version "<<myFile->getFormatVersion()<<endl;
            cout<<"Satellite Sistem: "<<myFile->getSatSystem()<<endl;
            cout<<"Rinex Type: "<<myFile->getRinexFileType()<<endl;
            cout<<"Rinex Program: "<<myFile->getRinexProgram()<<endl;
            cout<<"Created by Agency: "<<myFile->getCreatedByAgency()<<endl;

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
        }//try
        catch( RequiredRecordMissingException &headerExcep )
        {
            cout << " RequiredRecordMissingException is: " << endl
                 << headerExcep.getMessage() << endl;
        }//catch
        while( myobs.readEpoch( currentObsEpoch ) != 0 ){   //read epochs
            cout<< "test"<<endl;
        }
    }
}

void rinex_ho_Tests::test_RINEX2obs() {
    std::string filenameobs="RECF/rjni2361.21o";//"RECF/recf010c.16o";
    RinexFile *myFile = new RinexFile();   //object from class RinexFile
    double X0(0.0),Y0(0.0),Z0(0.0),xest(0.0),yest(0.0),zest(0.0),lat(0.0),lamb(0.0),h(0.0); //Station coordinates (read in the Header of RINEX file)
    string theObsFileType;
    try {  //try to open rinex file
            myFile->setPathFilenameMode( filenameobs, ios::in );
            cout<<" - Version "<<myFile->getFormatVersion()<<endl;
            cout<<"Satellite Sistem: "<<myFile->getSatSystem()<<endl;
            cout<<"Rinex Type: "<<myFile->getRinexFileType()<<endl;
            cout<<"Rinex Program: "<<myFile->getRinexProgram()<<endl;
            cout<<"Created by Agency: "<<myFile->getCreatedByAgency()<<endl;

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
        }//try
        catch( RequiredRecordMissingException &headerExcep )
        {
            cout << " RequiredRecordMissingException is: " << endl
                 << headerExcep.getMessage() << endl;
        }//catch
        while( myobs.readEpoch( currentObsEpoch ) != 0 ){   //read epochs

            cout<< "Sats:" << currentObsEpoch.getNumSat()<<endl;
        }
    }
}


void rinex_ho_Tests::test_RINEX2_glonass() {
    QBENCHMARK {
        ephemeris test;
        vector<GlonassEphemEpoch> glonassData;
        RinexFile* mynav= new RinexFile();
        int nephs;
        fstream log;
        test.ReadGlonassNavFile(glonassData,mynav,"RECF/rjni2361.21g",nephs,log);
        //test.readGlonassNav("RECF/rjni2361.21g");
        Vector3d newPos=test.getSatPos(172800+30*60-1, 1, glonassData);
        std::cout << std::setprecision(3);
        std::cout << std::fixed;
        //cout<<"SP3 (0min GPS?): PR01 -23925811.092   8776865.852    922070.055 "<<endl;
        cout<<"Previous pos: "<< -2.357955810547e+07 <<", "<< 8.771120605469e+06 <<", "<< 4.170250976563e+06 <<endl;
        cout<<"Integrated : PR01 "<< newPos[0] <<" "<< newPos[1] <<" "<< newPos[2]<<endl;
        //cout<<"SP3 (15min): PR01 -23589901.990   8773088.620   4107160.190"<<endl;
        //test.getSatPos(175000, 1, glonassData);
        cout<<"Next pos: "<<-2.190240869141e+07 <<", "<< 8.114974609375e+06 <<", "<<1.023487744141e+07<<endl;
    }
}

QTEST_APPLESS_MAIN(rinex_ho_Tests)
#include "test.moc"

