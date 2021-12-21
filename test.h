#ifndef TEST_H
#define TEST_H
#include <QtTest>
#include <iostream>
#include "ephemeris.h"
#include <iomanip>      // std::setprecision
using namespace std;
#include "rinex.h"
using namespace NGSrinex;

class rinex_ho_Tests : public QObject
{
    Q_OBJECT

public:
    rinex_ho_Tests();
    ~rinex_ho_Tests();

private slots:
    void test_RINEX2();
    void test_RINEX3nav();
    void test_RINEX3obs();
    void test_RINEX2obs();
    void test_RINEX2_glonass();
};

#endif
