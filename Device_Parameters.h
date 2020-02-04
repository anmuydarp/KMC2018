/********************************************************************/
/*  Device Parameters for Kinetic Monte Carlo  (i-aSi)              */
/********************************************************************/

#ifndef DPARAMETERS_H
#define DPARAMETERS_H

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double Device_X=10*1E-9; //in m (10 nm)
const double Device_Y=100*1E-9; //in m (10 um)
const double Device_Z=100*1E-9; //in m (10 um)

// Band-Edge Description
const double BandSlope=-0.0317;
const double BandIntercept=0.4470; // This is the band offset
const double QFLSlope=0.0;
const double QFLIntercept=-0.01419;
const double aSi_EField=3.15e7; // V/m

const double VBOffset=BandIntercept; // in eV
const double VBDelta=0.7;
const double ithickness=Device_X; // in m
 
const double Volume=1E-9*Device_Y*Device_Z; // Volume of slice

const double dx=1E-9; // in m
const double dy=1E-6; // in m
const double dz=1E-6; // in m

const int xpoint=Device_X/dx;
const int ypoint=Device_Y/dy;
const int zpoint=Device_Z/dz;

const double Nvt=2E21; //in per cm3
const double Evt=0.045; //in eV (I think this value may be 45 meV but need to find a reference for it)

const int CarrierMax=10000;
const double EnergyStep=1E-3;

// Amorphous Silicon Parameters
const double Eg_aSi=1.72; // eV
const double m_aSi=1; // Effective Hole Mass a-Si Double Check value
const double m_Si=0.81; // Effective Hole Mass Si
const double homega=0.060; // Phonon coupling constant for a-Si
const double Ud=8.51E-49; // Deformation Potential in J^2
const double S=0.5; // Huang-Rhys Factor
const double LatticeVibration=1E13; // per second
const double eopt_aSi=5.6; //
const double Nc_aSi=2.5E20;
const double Nv_aSi=2.5E20;
const double Ni_aSi=1E6;
const double ElecMob_aSi=5; // cm2/V.sec
const double HoleMob_aSi=1; // cm2/V.sec

// Crystalline Silicon Parameters
const double Eg_Si=1.12; // eV
const double Nc_Si=1.05E19;
const double Nv_Si=3.92E18;
const double ni_Si=1E10;

extern double EGlobal,ETrans,initialx,GammaMax,GammaMaxIR,Grdi,Grdj,GRi,GRj,Gx,Gd,Gradi,Gradj,Galpha,finalx;
extern int TotalLocalDefects,TotalMGDefects,StuckCarrierCount;
extern int nDefect,PathCount;
extern int PFCount,D2ECount,ThermCount,TotPh;
extern int ph0,ph1,ph2,ph3,ph4,ph5,Ab,Em,El;
extern double AvgE,AvgTau,AvgEpen,Avgx1;
extern int Refl,Refl_Sval,Refl_Occ,Refl_QMTC;
extern double GlobF;

// Global variables for bulk recombination
extern double MG_gamma,MGE,MG2E,MGH,tauB_n0,tauB_p0;
extern int MGHop,HEm,ECap;

// Global variables for interface recombination
extern int IntTotal,selectIR,IRcount;
extern double IntElec,IntHole,tauI_n0,tauI_p0;

// Parameters for Mid gap states
const double H=5E21*0; // Hydrogen concentration in the aSi - 10%
const double Efermi=0.3; // 1.05 eV for bulk i-aSi
const double Efermi_Dev=0.0; //Device Fermi Level
const double sigma=0.190; // Defect Pool width
const double U=0.2; // eV - Correlation energy
const double delta=0.44; // eV
const double Eprob=Efermi+0.5*delta; // eV
const double NSiSi=2E23;
const double Vth=2E7; // cm/sec
const double sigma0=1E-16; // cm^-2 // Capture cross section area for neutral states D0
const double sigma_plus=2.5E-15; // Capture cross section area for positively charged states D+
const double sigma_minus=5E-14; // Capture cross section area for negatively charges states D-

// Parameters for Interface Defect Density
const double Nint=1E12; // cm^-2, Interface defect density
const double sigma_Int=0.1;
const double Int_Min=-0.36; // Min Interface Defect Energy
const double Int_Max=-0.46; // Max Interface Defect Energy
const double Int_Mid=-0.41;

#endif
