#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <stdlib.h>
#include <csignal>
#include "Global.h"
#include "Device.h"
#include "Device_Parameters.h"
#include "Physics_Parameters.h"

using namespace std;

//----------Initializing Global Variables-------------------
double EGlobal,ETrans,initialx,GammaMax,GammaMaxIR,Grdi,Grdj,GRi,GRj,Gx,Gd,Gradi,Gradj,Galpha,finalx;
int nDefect,PathCount=0,PFCount=0,D2ECount=0,ThermCount=0,Refl=0,Refl_Sval=0,Refl_Occ=0,Refl_QMTC=0;
int TotPh=0,ph0,ph1,ph2,ph3,ph4,ph5,Ab,Em,El;
double AvgE,AvgTau,AvgEpen,Avgx1,GlobF=0;
int TotalLocalDefects=0,TotalMGDefects=0,StuckCarrierCount=0;

//----------- Global variables for Interface Recombination a-Si/c-Si -----------
int IntTotal=0,selectIR,IRcount=0;
double IntElec,IntHole,tauI_n0,tauI_p0;

//-------------- Global variables for Bulk Recombination in a-Si ---------------
double MG_gamma;
double MGE,MG2E,MGH,tauB_n0,tauB_p0;
int MGHop=0,HEm=0,ECap=0;

int main(void)
{
	// Initialize all classes
	Defect Def;
	CarrierH Carr;
	RateTable TRT;
	
	srand (time(NULL));

	Def.DefectInit();
	Carr.CarrierInit();	

	ofstream myfile;
	
	int per;
	
	// Read Valence Band Profile from SILVACO (Don't require it right now, analyzing a 10 nm slice. Might want to represent slope more accurately)
	Initialize_Defects(Def);

	cout<<"MGE : "<<MGE<<"   "<<"MG2E : "<<MG2E<<"   "<<"MGH : "<<MGH<<'\n';
   
	//CalculateOccupation(Def);		    

	// Initialize Carriers (Non Maxwellian Distribution or Gaussian Distribution)
	Initialize_Carriers(Carr);

	// Create and Normalize Initial Rate of Capture Transition Rate Table
	TRT.TRTInit();
	iROC(Def,TRT);
	iROCSum(TRT);
	NormalizeRCTable(TRT,Def);
	cout<<'\n'<<"IROC COMPILED"<<'\n';
	//InterfaceRecombination(Def,TRT);
	//NormalizeIRTable(TRT,Def);
	//cout<<'\n'<<"INTERFACE MECH COMPILED"<<'\n';	

	// KMC Loop
	KMCLoop(Def,TRT,Carr); // The KMC version of freeflight and scatter	

	//Calculate_SRH_Interface_Lifetime(Def,Carr,0);

	//Calcualte_SRH_Bulk_Lifetime(Def);		
	
	// Write Final Positions
	WriteFinalPos(Carr);

	// Calculating Averages
	CalcAverage(Carr);

	per=Refl/double(CarrierMax);

	myfile.open("SimulationRun.txt");
	{
		myfile<<"PFCount : "<<PFCount<<'\n';
        myfile<<"D2E : "<<D2ECount<<'\n';
        myfile<<"Thermionic Emission : "<<ThermCount<<'\n';
		myfile<<"0 Phonon : "<<ph0<<'\n';
        myfile<<"1 Phonon : "<<ph1<<'\n';
        myfile<<"2 Phonon : "<<ph2<<'\n';
        myfile<<"3 Phonon : "<<ph3<<'\n';
        myfile<<"4 Phonon : "<<ph4<<'\n';
        myfile<<"5 Phonon : "<<ph5<<'\n';
		myfile<<"Total Absorption : "<<Ab<<"   "<<"Total Emission : "<<Em<<'\n';
		myfile<<"Avg Energy : "<<AvgE<<'\n';
        myfile<<"Avg Time : "<<AvgTau<<'\n';
	}

	//Def.DefectDel();
	//Carr.CarrierDel();
	//TRT.TRTDel();

    cout<<"Stuck Carrier Count : "<<StuckCarrierCount<<'\n';

	cout<<"Kinetic Monte Carlo"<<'\n';
}
