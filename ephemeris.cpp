#include "ephemeris.h"

ephemeris::ephemeris()
{

}

int ephemeris::readGlonassNav(string navFile)
{

}

int ephemeris::ReadGlonassNavFile(vector<GlonassEphemEpoch> & currentPRNBlock, RinexFile *myFileNav, string filenamenav,
                int &cont_efe_nav, fstream &log)
{/*-------------------------------------------------------------------
    Purpose: Read broadcast orbit file
    -------------------------------------------------------------------
    Input:

    Output:
    -------------------------------------------------------------------
    Authors: Maurício C. M. de Paulo

    Date: Nov 2021
    -------------------------------------------------------------------
    Observation:
    -------------------------------------------------------------------*/
    string theNavFileType;

    try
    {   //opening file
        myFileNav->setPathFilenameMode( filenamenav, ios::in );

        //type of file
        theNavFileType = myFileNav->getRinexFileType();
        delete myFileNav;
        myFileNav = 0; //TODO: This should be replaced by a real representation of the rinex file.
    }//try
    catch( RinexFileException &openExcep )
    {
        cout << "Error opening file: " << filenamenav << endl
             << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
             << "Exiting program... " << endl << endl;

        log << "Error opening file: " << filenamenav << endl
            << "Rinex File\"Exception:\" " << endl << openExcep.getMessage() << endl
            << "Exiting program..." << endl << endl;

        exit (-1);
    }//catch


    //If the file is an NAV file
    if( theNavFileType[0] == 'G' )
    {
        GlonassNavFile mynav;

        try
        {
            mynav.setPathFilenameMode(filenamenav,ios_base::in);
        }
        catch( RinexFileException &openExcep )
        {
            cout << "Error opening file: " << filenamenav << endl;
            cout << " RinexFileException is: " << endl <<
            openExcep.getMessage() << endl;
        }//catch

        try
        {
            mynav.readHeader();

        }
        catch( RequiredRecordMissingException &headerExcep )
        {
            cout << " RequiredRecordMissingException is: " << endl <<
            headerExcep.getMessage() << endl;
        }

        //read the PRN Blocks
        try
        {   cont_efe_nav=0; //number of ephemeris block
            GlonassEphemEpoch tempBlock;
            while( mynav.readEphemEpoch( tempBlock) != 0 )
            {
                currentPRNBlock.push_back(tempBlock);
                cont_efe_nav++;
            }//while
        }//try
        catch( RinexReadingException &readingExcep )
        {
            cout << " RinexReadingException is: " << endl <<
            readingExcep.getMessage() << endl;
        }//catch

        log << endl << endl;
        log << "Number of OBS with WARNINGS = " << mynav.getNumberWarnings() << endl;
        log << endl << endl;
        log << "*** Messages ERRORS:" << endl;
        log << mynav.getErrorMessages() << endl;
        log << "*** Messages WARNINGS:" << endl;
        log << mynav.getWarningMessages() << endl;


    }//if

return (1);

}

double ephemeris::getGlonassF1(GlonassEphemEpoch &epoch) {
    return 1602 + epoch.getFreqNumber() * 9./16.;
}

double ephemeris::getGlonassF2(GlonassEphemEpoch &epoch) {
    return 1246 + epoch.getFreqNumber() * 7./16.;
}

double ephemeris::norm(vector<double> v) {
    double res=0;
    for (auto it=v.begin();it!=v.end();++it) {
        res+=(*it)*(*it);
    }
    return sqrt(res);
}

VectorXd ephemeris::glo_deriv(double tt, VectorXd vx, Vector3d acc) {
    //(Statella et al, 2013)   Calculo dos vetores de posicao e velocidade dos satelites GLONASS a partir das efemerides transmitidas e aspectos relacionados a sua integracao com o GPS

    static const double gmWGS = 398.60044e12;
    static const double AE    = 6378136.0;
    static const double OMEGA = 7292115.e-11;
    static const double C20   = -1082.6257e-6;
    VectorXd deriv(6);
    Vector3d pos=vx.segment(0,3);

    double rho = pos.norm();
    double t1  = -gmWGS/(rho*rho*rho);
    double t2  = 3.0/2.0 * C20 * (gmWGS*AE*AE) / (rho*rho*rho*rho*rho);
    double t3  = OMEGA * OMEGA;
    double t4  = 2.0 * OMEGA;
    double z2  = pos[2] * pos[2];


    deriv[0] = vx[3];
    deriv[1] = vx[4];
    deriv[2] = vx[5];
    deriv[3] = (t1 + t2*(1.0-5.0*z2/(rho*rho)) + t3) * vx[0] + t4*vx[4] + acc[0];
    deriv[4] = (t1 + t2*(1.0-5.0*z2/(rho*rho)) + t3) * vx[1] - t4*vx[3] + acc[1];
    deriv[5] = (t1 + t2*(3.0-5.0*z2/(rho*rho))     ) * vx[2]            + acc[2];

    return deriv;
}

VectorXd ephemeris::rk4_integrate(double xi,
const VectorXd yi, // vector of the initial y-values
double dx,              // the step size for the integration
Vector3d acc            // additional acceleration
) {
//return yi;

VectorXd k1 = glo_deriv(xi       , yi       , acc) * dx;
VectorXd k2 = glo_deriv(xi+dx/2.0, yi+k1/2.0, acc) * dx;
VectorXd k3 = glo_deriv(xi+dx/2.0, yi+k2/2.0, acc) * dx;
VectorXd k4 = glo_deriv(xi+dx    , yi+k3    , acc) * dx;

VectorXd yf = yi + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

return yf;
}
Vector3d ephemeris::integrate(double dt, Vector3d &pos, Vector3d &velocity, Vector3d &acc, double stepSecs)
{

    int steps  = int(fabs(dt) / stepSecs) + 1;
    double step = dt / steps, tt=0.;
    VectorXd currPos(6);
    currPos.segment(0,3)=pos;
    currPos.segment(3,3)=velocity;
    for (int ii = 1; ii <= steps; ii++) {
      currPos = rk4_integrate(tt,currPos, step, acc);
      tt = tt + step;
    }



    return currPos.segment(0,3);
}

Vector3d ephemeris::getSatPos(double secsOfWeek, int SV, vector<GlonassEphemEpoch> navEpochs)
{
    Vector3d Pos,Velocity,Acc,deriv, newPos;
    //find  latest epoch
    auto epoch=navEpochs.begin();
    for (epoch=navEpochs.begin(); epoch!=navEpochs.end();++epoch){
        cout<<epoch->getMessageFrameTime()<<endl;
        if ((epoch->getSatelliteAlmanacNumber()==SV
             && secsOfWeek >= epoch->getMessageFrameTime())
             && secsOfWeek < (epoch+1)->getMessageFrameTime() ){
            break;
        }

    }
    if (epoch!=navEpochs.end()){

        //double dt=secsOfWeek-epoch->getMessageFrameTime()-19;
        double dt=secsOfWeek-epoch->getMessageFrameTime();
        double GammaN=epoch->getSvRelFreqBias();
        double TauN=epoch->getSvClockBias();
        double TauC=0;
        //double TauC=mynav.getTimeScaleCorr();
        double Tb=0;
        double TUTC = secsOfWeek + TauN - GammaN*(dt) + TauC;
        double t_glo=dt;
        //double leapSeconds=mynav.getLeapSec();


        Pos[0]=epoch->getPosX()*1e3;
        Pos[1]=epoch->getPosY()*1e3;
        Pos[2]=epoch->getPosZ()*1e3;
        Velocity[0]=epoch->getVelX()*1e3;
        Velocity[1]=epoch->getVelY()*1e3;
        Velocity[2]=epoch->getVelZ()*1e3;
        Acc[0]=epoch->getAccX()*1e3;
        Acc[1]=epoch->getAccY()*1e3;
        Acc[2]=epoch->getAccZ()*1e3;
        newPos=integrate(t_glo,Pos,Velocity,Acc);
        return newPos;
    }
    /*


    double newPos[3];

    double X[3]={navEpoch.getPosX(), navEpoch.getPosY(), navEpoch.getPosZ()};
    double Vel[3]={navEpoch.getVelX(), navEpoch.getVelY(), navEpoch.getVelZ()};
    double Acc[3]={navEpoch.getAccX(), navEpoch.getAccY(), navEpoch.getAccZ()};*/

}
