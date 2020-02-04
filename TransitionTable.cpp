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

// Calculate Occupation Function
void CalculateOccupation(Defect& Def)
{
	int i,count,j=0;
	double *EnergyRep,*EnergyOccupied,*EnergyUnOcc;

	EnergyRep=new double[451];
	EnergyOccupied=new double[451];
	EnergyUnOcc=new double[451];

	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		EnergyRep[i]=0;
		EnergyOccupied[i]=0;
		EnergyUnOcc[i]=0;		
	}

	// Calculating total number of occupied and unoccupied states, number of states present at the same energy
	for (i=0;i<nDefect;i++)
	{
		count=(Def.Epos[i]/EnergyStep);

		EnergyRep[count]=EnergyRep[count]+1;

		if (Def.FlagOccupation[i]==0)
			EnergyOccupied[count]=EnergyOccupied[count]+1;
		else
			EnergyUnOcc[count]=EnergyUnOcc[count]+1;

	}

	cout<<'\n';
	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		cout<<"Energy : "<<i*EnergyStep<<"    "<<"Reps : "<<EnergyRep[i]<<"   "<<"Occupied : "<<EnergyOccupied[i]<<"    "<<"UnOccupied : "<<EnergyUnOcc[i]<<'\n';
	}

	delete [] EnergyRep;
	delete [] EnergyOccupied;
	delete [] EnergyUnOcc;
}

// Initial Rate of Caprture Function
void iROC(Defect& Def,RateTable& TRT)
{
	int i,j,count=0;

	initialx=0.0;

	ofstream myfile;

    /****************************************************************
    Setting up table for carrier injection from a-Si:H(i)/c-Si 
    heterointerface
    ****************************************************************/

	for (i=0;i<=VBOffset/EnergyStep;i++)	
	{
		for (j=0;j<nDefect;j++)	
		{
			if (Def.MGFLAG1D[count]==1) // Mid Gap states
				TRT.RCTable[i][j]=0;
			else // Localized States
				TRT.RCTable[i][j]=MPROC(i,j,Def,0);        

			count++;
		}
		count=0;
	}
}

void InterfaceRecombination(Defect& Def,RateTable& TRT)
{
	int i,j,count=0;
	double E,temp1,temp2,jointtemp,HoleConc=1E16;
	
	// TRANSITION TO D0 STATES	
	// Multi Phonon Assisted Transition + Electron Capture
	// Multi Phonon Assisted Transition + Hole Emission
	// Hole Capture + Hole Emission

	/**************************************************************************************************************
	RECOMBINATION VIA D0 STATES
	1. HCap via Mph + ECap
	2. HCap via Mph + HEm
	3. HCap + ECap
	4. HCap + HEm

	RECOMBINATION VIA D+ STATES
	5. ECap + HCap via Mph
	6. ECap + HCap
	7. ECap + EEm
	**************************************************************************************************************/

	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		for (j=0;j<IntTotal;j++)
		{			
			TRT.InterfaceMech[i][j][0]=MPROC(i,j,Def,3); // Hole Capture via Mph
			TRT.InterfaceMech[i][j][1]=RecombinationRate(1,Def.ID_xpos[j],Def.ID_Epos[j]); // Hole Capture
			TRT.InterfaceMech[i][j][2]=RecombinationRate(0,Def.ID_xpos[j],Def.ID_Epos[j]); // Hole Emission
			TRT.InterfaceMech[i][j][3]=RecombinationRate(3,Def.ID_xpos[j],Def.ID_Epos[j]); // Electron Capture
			TRT.InterfaceMech[i][j][4]=RecombinationRate(2,Def.ID_xpos[j],Def.ID_Epos[j]); // Electron Emission
		}
		// Max and normalize
		// Calculate joint probability
	}
}

// Dynamic Transition Rate Table for all defect mechanisms 
void HopTable(int i,Defect& Def,RateTable& TRT,CarrierH& Carr,int MG)
{
	int j,Loci,count=0;
	double VBmin,VBmax,E,xi; // xi - Initial position of carrier

	E=Carr.Energy[i];
	Loci=E/EnergyStep;
	xi=Carr.xpos[i];
	initialx=xi;

	// For Defect to Defect Transitions
	for (j=0;j<nDefect;j++)
	{
		if (Def.FlagOccupation[j]==0 || Def.FlagOccupation[j]==3) // 0 = 1 electron , 3 = 2 electrons
		{
			TRT.D2D[j]=D2D(E,xi,j,Def);
			//TRT.D2D[j]=MPROC(Loci,j,Def,2);
		}
		else
		{
			TRT.D2D[j]=0;			
		}
	}

    VBmin=round(BandEdge(Device_X/dx,0,0)*1e3)/1e3;
    VBmax=round(BandEdge(0,0,0)*1e3)/1e3;
    
    /*cout<<"VBmin : "<<VBmin<<"   "<<"VBmax : "<<VBmax<<'\n';
    cout<<"int min : "<<int(VBmin/EnergyStep)<<"   "<<"int max : "<<int(VBmax/EnergyStep)<<'\n';
    cout<<"min : "<<VBmin/EnergyStep<<"   "<<"max : "<<VBmax/EnergyStep<<'\n';*/

	// For Defect to Electrode
	for (j=(VBmin/EnergyStep)+1;j<VBmax/EnergyStep;j++)
	{
		//TRT.D2E[j]=MPROC(Loci,j,Def,1); //Initial Energy,Final Energies,Defect Class
        TRT.D2E[count]=MPROC(Loci,j,Def,1); // MPROC ( Initial Energy, Final Energy, Defect Class, arg to toggle Defect->Electrode )
        count++;
	}   

	// EfSum = Summation of all final states	
	TRT.S1Table[0]=PF(Loci); // Poole Frenkel Emission    

	TRT.S1Table[1]=EfSum(Def,TRT,1); // Defect to Electrode Emission   

	TRT.S1Table[2]=EfSum(Def,TRT,2); // Defect to Defect Mechanism  

	//TRT.S1Table[3]=RecombinationRate(0,xi,E); //  0 - Hole Emission

	//TRT.S1Table[4]=RecombinationRate(5,xi,E); // 1 - Electron Capture

	TRT.TrackRate[0][i]=TRT.S1Table[0];
	TRT.TrackRate[1][i]=TRT.S1Table[1];
	TRT.TrackRate[2][i]=TRT.S1Table[2];

	NormalizeS1Table(TRT,MG);
}
