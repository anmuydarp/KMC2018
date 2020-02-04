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

void fROC(int i,Defect& Def,RateTable& TRT,CarrierH& Carr)
{
	int j,Loci,occuprev,occunext,Flag=0,find=0,locate,State,QMTC,count=0,Therm_Toggle=0;
	double r,tau,Ein,Ef,xin;
	
	Ein=Carr.Energy[i];	
    Loci=Ein/EnergyStep;
	xin=Carr.xpos[i];

	//---------------Initializing Tracking--------------------
	Carr.xPath[i][PathCount]=xin;
	Carr.EPath[i][PathCount]=Ein;
	Carr.TauPath[i][PathCount]=0;
	//---------------Initializing Tracking--------------------
	
	while (Flag==0)
	{
		r=(double)rand()/RAND_MAX;

        if (Carr.Energy[i]<=VBOffset)
        {
            Therm_Toggle=0;
            State=BinarySearch(TRT,Loci,r,0);
        }
        else
        {            
            Therm_Toggle=1;
        }

		if (State==0)
			raise(SIGSEGV);
        
        switch (Therm_Toggle)
        {
            case 0 : // Electrode -> Defect Emission

	    	//if (Def.FlagOccupation[State]==1) // Occupied by Hole
		    if (Def.FlagOccupation[State]==0) // Occupied by Electron
    		{
	    		Flag=0;
		    	//cout<<"State Rejected Due to Occupation"<<'\n';
			    Refl++;
    			Refl_Occ++;
                count++;

                if (count==15000) 
                {
                    /******************************************************
                    1. Carrier is stuck at the interface
                    2. Unable to find a final state 
                    3. This carrier should be selected for recombination
                    ******************************************************/
				    Carr.xPath[i][PathCount]=0;
    				Carr.EPath[i][PathCount]=0;
                    Carr.TauPath[i][PathCount]=0;
                    Carr.Status[i]=4;
                    Flag=1;
                    cout<<"Carrier Excluded"<<'\n';
                }		
	    	}
    		//else if (Def.FlagOccupation[State]==0 || Def.FlagOccupation[State]==3) // Occupied of Electron
	    	else if (Def.FlagOccupation[State]==1) // Occupied by Hole
		    {
			    Ef=Def.Epos[State];
    			locate=ValidateFinalState(Ein,Ef);

	    		if (locate==0) // INVALID STATE
		    	{
			    	Flag=0;
				    //cout<<"State Rejected Due to Invalid State"<<'\n';
    				Refl++;
	    			Refl_Sval++;		
		    	}
			    else if (locate==1) // VALID STATE
    			{
	    			QMTC=QMTCReflection(Ein,Def.xpos[State]);
		    		QMTC=0;

			    	if (QMTC==0) // Transmission
				    {
					    if (Def.MGFLAG1D[State]==1)
    					{
		    				Carr.Status[i]=3; // Captured by a Mid Gap State
	    				}
			    		else
				    	{
					    	Carr.Status[i]=1; // Captured by defect
                            //Carr.Status[i]=2;
    					}	
					
	    				Carr.Energy[i]=Def.Epos[State];
		    			Ef=Carr.Energy[i];
			    		Carr.xpos[i]=Def.xpos[State];
				    	occuprev=Def.FlagOccupation[State];
					    tau=-(1/TRT.GammaRC[Loci])*log(r);

    					//---------------TRACKING---------------------
		    			PathCount++;
	    				Carr.xPath[i][PathCount]=Carr.xpos[i];
			    		Carr.EPath[i][PathCount]=Carr.Energy[i];
                      	Carr.TauPath[i][PathCount]=tau;
					    //---------------TRACKING---------------------
    					//cout<<"Electrode to Defect Injetion : "<<i<<"   "<<"E : "<<Ein<<" -> "<<Ef<<" "<<"xin : "<<xin<<"   "<<"xf : "<<Carr.xpos[i]<<"  "<<"Gamma : "<<TRT.GammaRC[Loci]<<"  "<<"tau : "<<tau<<"  "<<"State : "<<State<<'\n';	                                  
                    
                        if (isinf(tau)>0)
                            raise(SIGSEGV);
                        else
                            Flag=1;

					    TotPh++;
    				}
	    			else if (QMTC==1) // Reflection
		    		{
					Flag=0;					
					//cout<<"State Rejected Due to QM Reflection"<<'\n';
					Refl++;
					Refl_QMTC++;
				}				
				break;
			}
            break;
    
            case 1: // Thermionic Emission
            Carr.Energy[i]=Ein; // Energy remains the same
            Carr.xpos[i]=Device_X;
        
    		//---------------TRACKING---------------------
		    PathCount++;
	    	Carr.xPath[i][PathCount]=Carr.xpos[i];
			Carr.EPath[i][PathCount]=Carr.Energy[i];
            Carr.TauPath[i][PathCount]=Device_X/(HoleMob_aSi*1e-4*aSi_EField);
			//---------------TRACKING---------------------

            ThermCount++;
            break;            
            }
		}
	}
	
	if (Ef<Ein) // Emission
		Em++;
	else if (Ef>Ein) // Absorption
		Ab++;
	else if (Ef==Ein) // Elastic
		El++;
	//cout<<"Final State for Rate of Capture"<<'\n';	
}

int MechSelect(RateTable& TRT,int MG)
{
	//cout<<"Entering Stage 1 of Selection : Selecting Mechanism"<<'\n';
	double r;
	int i,select,Max;

	//Stage one selection using TRT.S1Table
	r=(double)rand()/RAND_MAX;

	/********************************************************************
	Transiton 0: Poole Frenkel Emission
	          1: Defect to Electrode Emission
	          2: Defect to Defect Emission
	          3: Recombination via Hole Emission
	          4: Recombination via Electron Capture
	********************************************************************/

	switch (MG)
	{
		case 0: // Without mid gap states
		Max=2;
		break;

		case 1: // With mid gap states
		Max=4;
		raise(SIGSEGV);
		break;
	}	          

	for (i=0;i<=Max;i++)
	{		
		if (i==0)
		{
			if (r<=TRT.S1Table[i])
			{
				select=i;
				break;
			}
		}
		else if (i>0)
		{
			if (r>TRT.S1Table[i-1] && r<=TRT.S1Table[i])
			{
				select=i;
				break;
			}
		}	
	}

	// KEEPING COUNT OF EXTRACTION MECHANISMS
	switch(select)
	{
		case 0: // Poole Frenkel Emission
		PFCount++;		
		break;

		case 1: // Defect to Electrode Emission
		D2ECount++;
		break;
	}

	//cout<<"Selected Mech (1st Select) : "<<select<<"    "<<"r : "<<r<<'\n';	

	return select;
	//Stage two selection by using TRT.<mechanism> table
}

void FinalState(int i,int Select,Defect& Def,RateTable& TRT,CarrierH& Carr)
{
	/****************************************************************
	Select 0 - Poole Frenkel Emission
	       1 - Defect to Electrode Emission
	       2 - Defect to Defect Emission
	****************************************************************/

	int j,occuprev,Flag=0,locate=0,State=0,count;
	double r,Ein,Ef,xin,tau;

	Ein=Carr.Energy[i];
	xin=Carr.xpos[i];   

	switch (Select)
	{
/*****************************************************************************
							POOLE FRENKEL EMISSION
*****************************************************************************/		
		case 0: // Poole Frenkel Transition		*******FIX THIS****************
		Carr.Status[i]=2; 		
		Carr.Energy[i]=BandEdge(Carr.xpos[i]/dx,0,0);
        Carr.xpos[i]=Device_X;
        r=(double)rand()/RAND_MAX;
		tau=-(1/GammaMax)*log(r);

		//---------------TRACKING---------------------
		PathCount++;
		Carr.xPath[i][PathCount]=Carr.xpos[i];
		Carr.EPath[i][PathCount]=Carr.Energy[i];       
		Carr.TauPath[i][PathCount]=tau;
		Carr.Epen[i]=Ein;
		//---------------TRACKING---------------------
		//cout<<"PF Emission"<<"  "<<"Ein : "<<Ein<<"   "<<"Ef : "<<Carr.Energy[i]<<"   "<<"tau : "<<tau<<'\n';
		break;

/*****************************************************************************
						DEFECT TO ELECTRODE EMISSION
*****************************************************************************/						

		case 1: // Defect to Electrode Emission

		NormalizeS2Table(TRT,1);

        count=0;

		while (Flag==0)
		{
            r=(double)rand()/RAND_MAX;
			
			State=BinarySearch(TRT,0,r,1);

			Ef=BandEdge(Device_X/dx,0,0)+State*EnergyStep; // The energy has to offset by the minima of the valence band edge

			locate=ValidateFinalState(Ein,Ef);         
         
            if (locate==0)
            {   
                if (count>1e5)
                {
                    Carr.Status[i]=4;
                    Flag=1;
                    cout<<"Unable to find final state"<<'\n';
                    raise(SIGSEGV);
                    // No final state available
                }
                count++;
            }
            else if (locate==1)
			{
                Carr.Energy[i]=(State)*EnergyStep;
    			Carr.xpos[i]=Device_X;
	    		tau=-(1/GammaMax)*log(r);
		    	Carr.Status[i]=2;
			    //---------------TRACKING---------------------
    			PathCount++;
	    		Carr.xPath[i][PathCount]=Carr.xpos[i];
		    	Carr.EPath[i][PathCount]=Carr.Energy[i];
			    Carr.TauPath[i][PathCount]=tau;        
    			Carr.Epen[i]=Ein;
	    		//---------------TRACKING---------------------
		    	//cout<<"D2E : "<<"Ein : "<<Ein<<" "<<"Ef : "<<Carr.Energy[i]<<" "<<"xin : "<<xin<<" "<<"xf : "<<Carr.xpos[i]<<"  "<<"tau : "<<tau<<'\n';
			    TotPh++;
    			Flag=1;
            }
		}		
		break;
/*****************************************************************************
							DEFECT TO DEFECT EMISSION
*****************************************************************************/				

		case 2: // Defect to Defect Emission

		NormalizeS2Table(TRT,2);

        count=0;

		while (Flag==0)
		{
            r=(double)rand()/RAND_MAX;
			
            State=BinarySearch(TRT,0,r,2);
			   
			Ef=Def.Epos[State];

			locate=ValidateFinalState(Ein,Ef);

            if (locate==0)
            {
                if (count>1e5)
                {
                    Carr.Status[i]=4;
                    Flag=1;
                    // No final state available
                }
                count++;
            }
            else if (locate==1)
            {

			//Carrier Status does not change
    		occuprev=Def.FlagOccupation[State];
			Carr.Energy[i]=Def.Epos[State];
			Carr.xpos[i]=Def.xpos[State];
			tau=-(1/GammaMax)*log(r);
			Carr.Status[i]=1;
			//Def.FlagOccupation[State]=1;

			//---------------TRACKING---------------------
			PathCount++;
			Carr.xPath[i][PathCount]=Carr.xpos[i];
			Carr.EPath[i][PathCount]=Carr.Energy[i];
			Carr.TauPath[i][PathCount]=tau;
			//---------------TRACKING---------------------

			if (Def.MGFLAG1D[State]==1)
			{
				Carr.Status[i]=3; // Captured by a Mid Gap State					
			}

			//cout<<"D2D : "<<"Ein : "<<Ein<<" "<<"Ef : "<<Carr.Energy[i]<<" "<<"xi : "<<xin<<"  "<<"xf : "<<Carr.xpos[i]<<"  "<<"tau : "<<tau<<'\n';
			TotPh++;
			Flag=1;
            }
			if (Def.FlagOccupation[State]==1)
			{
				Flag==0;
			}
		}

		break;

		case 3: // Hole Emission for Recombination
		cout<<'\n'<<"HOLE EMISSION"<<'\n';
		PathCount++;
		Carr.Energy[i]=Carr.Energy[i]; // Trying to figure at what energy level does recombination happen
		//Carr.xpos[i]=Carr.xpos[i];
		r=(double)rand()/RAND_MAX;
		tau=-(1/GammaMax)*log(r);
		Carr.TauPath[i][PathCount]=tau;
		Carr.Status[i]=4;
		Carr.tauHEm[HEm]=tau;
		HEm++;
		//raise(SIGSEGV);
		break;

		case 4: // Electron Capture for Recombination
		cout<<'\n'<<"ELECTRON CAPTURE"<<'\n';
		PathCount++;
		Carr.Energy[i]=Carr.Energy[i];
		Carr.xpos[i]=Carr.xpos[i];
		r=(double)rand()/RAND_MAX;
		tau=-(1/GammaMax)*log(r);
		Carr.TauPath[i][PathCount]=tau;
		Carr.Status[i]=5;
		Carr.tauECap[ECap]=tau;
		ECap++;
		//raise(SIGSEGV);
		break;

		case 9: // No Transition was selected
		//cout<<"No Transition was selected"<<'\n';
		break;
	}
	
	if (Ef<Ein) // Emission
		Em++;
	else if (Ef>Ein) // Absorption
		Ab++;
	else if (Ef==Ein) // Elastic
		El++;
	//raise(SIGSEGV);
}

void IRSelection(int i,Defect& Def,RateTable& TRT,CarrierH& Carr)
{
	int loci,j;
	double *IRMechSelect,sum=0,r;

	IRMechSelect=new double[7];

	loci=Carr.Energy[i]/EnergyStep;

	// Select IR Mech

	//TRT.IRMech[0]=IRSum(TRT,0,loci)+IRSum(TRT,3,loci); // HCap via Mph + ECap
	TRT.IRMech[0]=0; // HCap via Mph + ECap
	
	//TRT.IRMech[1]=IRSum(TRT,0,loci)+IRSum(TRT,2,loci); // HCap via Mph + HEm	
	TRT.IRMech[1]=0; // HCap via Mph + HEm	
	
	TRT.IRMech[2]=(RecombinationRate(3,0,0)+RecombinationRate(6,0,0)); // HCap + ECap	
	tauI_p0=1/RecombinationRate(3,0,0);
	tauI_n0=1/RecombinationRate(6,0,0);

	//TRT.IRMech[3]=IRSum(TRT,1,loci)+IRSum(TRT,2,loci); // HCap + HEm
	TRT.IRMech[3]=0; // HCap + HEm

	//TRT.IRMech[4]=parallel(IRSum(TRT,3,loci),IRSum(TRT,0,loci)); // ECap + HCap via Mph
	//cout<<'\n'<<"4 : "<<TRT.IRMech[4]<<'\n';

	//TRT.IRMech[5]=parallel(IRSum(TRT,3,loci),IRSum(TRT,1,loci)); // ECap + HCap	
	//cout<<'\n'<<"5 : "<<TRT.IRMech[5]<<'\n';

	// Normalize Mech table
	for (j=1;j<4;j++)
	{		
		TRT.IRMech[j]=TRT.IRMech[j-1]+TRT.IRMech[j];
	}	
	sum=TRT.IRMech[3];

	//cout<<'\n'<<"sum : "<<sum<<'\n';

	for (j=0;j<4;j++)
	{
		IRMechSelect[j]=TRT.IRMech[j]/sum;
		//cout<<"IRMechSelect : j : "<<j<<"    "<<IRMechSelect[j]<<'\n';
	}	

	GammaMaxIR=sum;	

	// Select Final Mechanism
	r=(double)rand()/RAND_MAX;
	for (j=0;j<4;j++)
	{
		if (j==0)
		{
			if (r<IRMechSelect[j])
			{
				selectIR=0;
				break;
			}			
		}
		else
		{
			if (r>IRMechSelect[j-1] && r<=IRMechSelect[j])
			{
				selectIR=j;
				break;
			}
		}
	}	
	
	delete [] IRMechSelect;
}

void InjSelect(RateTable& TRT,CarrierH& Carr,int i)
{
	double r,max,*ChooseMech;
	int loci,j,select;

	ChooseMech=new double[3];

	loci=Carr.Energy[i]/EnergyStep;

	// Calculate Max
	TRT.IRMech[2]=TRT.IRMech[2];
	max=TRT.RCSum[loci]+TRT.IRMech[2];

	ChooseMech[0]=TRT.RCSum[loci]/max; // Multiphonon Hop
	ChooseMech[1]=(TRT.RCSum[loci]+TRT.IRMech[2])/max; // Interface Recombination

	//cout<<"Injection : "<<ChooseMech[0]<<"    "<<"IR : "<<ChooseMech[1]<<'\n';

	r=(double)rand()/RAND_MAX;	

	if (r<ChooseMech[0])
		select=0;
	else if (r>ChooseMech[0] && r<=ChooseMech[1])
		select=1;

	cout<<'\n'<<"SELECT : "<<select<<'\n';
	if (select==1)
	{		
		Carr.tauI_n0[IRcount]=tauI_n0;
		Carr.tauI_p0[IRcount]=tauI_p0;
		//cout<<"tau_n0 : "<<Carr.tauI_n0[IRcount]<<"   "<<"tau_p0 : "<<Carr.tauI_p0[IRcount]<<'\n';
		IRcount++;
		//cout<<"Lifetime : "<<-log(r)/GammaMaxIR<<'\n';
	}
	delete [] ChooseMech;
}
