#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <string>
#include <vector>
#include "rinex.h"
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
using namespace NGSrinex;

class ephemeris
{
public:
    NGSrinex::GlonassNavFile glonassNav;
    ephemeris();
    int readGlonassNav(string navFile);
    int ReadGlonassNavFile(vector<GlonassEphemEpoch> &currentPRNBlock, RinexFile *myFileNav, string filenamenav, int &cont_efe_nav, fstream &log);
    Vector3d getSatPos(double secsOfWeek, int SV, vector<GlonassEphemEpoch> navEpochs);
    double getGlonassF1(GlonassEphemEpoch &epoch);
    double getGlonassF2(GlonassEphemEpoch &epoch);
    double norm(vector<double> v);
    VectorXd glo_deriv(double tt, VectorXd vx, Vector3d acc);
    Vector3d integrate(double dt, Vector3d &pos, Vector3d &velocity, Vector3d &acc, double stepSecs=1);

    VectorXd rk4_integrate(double xi, const VectorXd yi, double dx, Vector3d acc);
    //vector<double> integrate(double dt, vector<double>& pos, vector<double> &velocity, vector<double> &acc, double stepSecs=1);
};

#endif // EPHEMERIS_H
