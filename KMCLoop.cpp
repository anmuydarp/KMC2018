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

void KMCLoop(Defect& Def,RateTable& TRT,CarrierH& Carr)
{
	double Ei;
	int i,j=0,Flag=0,select;	

	// NEED TO KNOWS
	/************************************************************************************************
	- Transitions are only allowed to states which are 'not occupied' by a hole.
	- 0 : Occupied by an electron, 1 : Occupied by a hole.
	- Transitions are only allowed to Def.FlagOccupation[]=0.
	REMEMBER : In the KMCLoop function 'i' is the carrier number but this 'i' has to be translated to
	carrier energy for every other function.
	************************************************************************************************/
	for (i=0;i<CarrierMax;i++) // CarrierMax - lets pretend that there are only 10 carriers
	{
		Ei=Carr.Energy[i];        

		PathCount=0;

		cout<<'\n'<<"CARRIER IN TRANSIT : "<<i<<'\n';

		Flag=0;

		while (Flag==0)
		{			
		// The switch statement should run uptil the time the carrier has been collected, i.e. till Carr.Status==2
		switch (Carr.Status[i])
		{
			case 0: // At the barrier		
			//IRSelection(i,Def,TRT,Carr); // For Interface Recombination
			//InjSelect(TRT,Carr,i); // For Interface recombination
			fROC(i,Def,TRT,Carr); // Picking the inital transition for a carrier at the barrier (finalROC - final state after Rate of Capture)           
			break;

			case 1: // Hopping via Defects
			HopTable(i,Def,TRT,Carr,0);
			// i - Carrier Number, Def - Defect Class, TRT - Transition Rate Table, Carr - Carrier Class, integer - Add mid gap states
			// 0 - no MG, 1 - yes MG
			select=MechSelect(TRT,0);
			FinalState(i,select,Def,TRT,Carr);

			switch (select)
			{
				case 0: // PF
				Carr.FinalRate[i]=TRT.TrackRate[select][i];
				Carr.TransType[i]=select;
				break;

				case 1: // D2E
				Carr.FinalRate[i]=TRT.TrackRate[select][i];
				Carr.TransType[i]=select;
				break;
			}            
			break;

			case 2: // Collected
			Carr.TotPath[i]=PathCount;
			cout<<"CARRIER HAS BEEN COLLECTED"<<"   "<<"PathCount : "<<PathCount<<'\n';
			Flag=1;
			break;

			case 3: // Captured at a Mid Gap State
			HopTable(i,Def,TRT,Carr,0);
			select=MechSelect(TRT,0);	
			FinalState(i,select,Def,TRT,Carr);
			break;

			case 4: // Hole Emission and currently 'stuck' carrier
			Carr.TotPath[i]=PathCount;
			Flag=1;
            StuckCarrierCount++;
            cout<<"Stuck Carrier"<<'\n';           
			break;

			case 5: // Electron Capture
			Carr.TotPath[i]=PathCount;
			Flag=1;
			break;			
		}
		}		
	}
	//cout<<'\n'<<"Inside KMC Loop"<<'\n';
}
