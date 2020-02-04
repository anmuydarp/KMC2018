#ifndef GLOBAL_H
#define GLOBAL_H

#include "Device.h"

using namespace std;

/******************************************************************
                    Function in Initialize.cpp
******************************************************************/
void Initialize_Defects(Defect& Def);
void Initialize_Carriers(CarrierH& Carr);

/******************************************************************
                    Function in CreateDefects.cpp
******************************************************************/
void Initialize_RealSpace_Defects(Defect& Def);
void Initialize_Energy_Defects(Defect& Def);
//void Initialize_Occupation(Defect& Def);
void Initialize_Occupation_2(Defect& Def);
void Initilize_IntDef_Energy(Defect& Def);
void Initilize_IntDef_Occupation(Defect& Def);
double BandEdge(double,double,int);
double QFL(double);

/******************************************************************
                    Function in TransitionTable.cpp
******************************************************************/
void CalculateOccupation(Defect& Def);
void iROC(Defect& Def,RateTable& TRT);
void HopTable(int,Defect& Def,RateTable& TRT,CarrierH& Carr,int);
void test_function(Defect& Def);
void InterfaceRecombination(Defect& Def,RateTable& TRT);

/******************************************************************
                    Function in TransitionRateCalc.cpp
******************************************************************/
double MPROC(int,int,Defect& Def,int);
void NormalizeRCTable(RateTable& TRT,Defect& Def);
void NormalizeS1Table(RateTable& TRT,int);
double EfSum(Defect& Def,RateTable& TRT,int);
void NormalizeS2Table(RateTable& TRT,int);
int ValidateFinalState(double,double);
int BinarySearch(RateTable& TRT,int,double,int);
void NormalizeIRTable(RateTable& TRT,Defect& Def);
double IRSum(RateTable& TRT,int,int);
void iROCSum(RateTable& TRT);

/******************************************************************
              Function in TransitionRateFunctions.cpp
******************************************************************/
double InelasticTransition(double,double,double,int,int);
double ElasticTransition(double,double,double,int);
double PF(int);
double D2D(double,double,int,Defect& Def);
double PhD2D(int,int,int,double,Defect& Def);
double TCTest(double,double,double);
int QMTCReflection(double,double);
double RecombinationRate(int,double,double);
//double QMTCReflection(double,double);

/******************************************************************
              Function in Lifetime.cpp
******************************************************************/
void Calculate_SRH_Interface_Lifetime(Defect& Def,CarrierH& Carr,int);
void Calcualte_SRH_Bulk_Lifetime(Defect& Def);

/******************************************************************
                    Function in KMCLoop.cpp
******************************************************************/
void KMCLoop(Defect& Def,RateTable& TRT,CarrierH& Carr);

/******************************************************************
                    Function in KMCEvents.cpp
******************************************************************/
void fROC(int,Defect& Def,RateTable& TRT,CarrierH& Carr);
int MechSelect(RateTable& TRT,int);
void FinalState(int,int,Defect& Def,RateTable& TRT,CarrierH& Carr);
void IRSelection(int,Defect& Def,RateTable& TRT,CarrierH& Carr);
void InjSelect(RateTable& TRT,CarrierH& Carr,int);

/******************************************************************
                    Function in Write_Calc.cpp
******************************************************************/
void CalcAverage(CarrierH& Carr);
void WriteFinalPos(CarrierH& Carr);

/******************************************************************
                    Function in Math.cpp
******************************************************************/
double trapzd(double (*)(double),double,double,int);
double qsimp(double (*)(double),double,double);
float bessi0(float);
float bessi1(float);
float bessi(int,float);
int factorial(int);
double sinc(const double);
double parallel(double,double);

#endif
