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

double MPROC(int i,int j,Defect& Def,int t)
{

	double R,Ei,Ef,x,Difference,Ep,DefectDepth,EpTemp=0,EvTail; // F = V/m, hw = eV, S = Huang Rhys Factor
    int DefOcc;

	//Ei = Initial Energy, Ef = Final Energy
	//1 - For electrode to defect, Ei is the carrier energy and Ef is the final energy
	//2 - For defect to electrode, Ei is the defect energy and Ef is the final energy
	//3 - Defect to Defect
	//4 - Electrode to Interface defect, Ei is the carrier energy, Ef is the final defect energy	

	int p,PhononMode,Choice; // Integer phonon number	

	switch(t)
	{
		case 0: // INJECTION Electrode --> Defect 
		Ei=i*EnergyStep;
		Ef=Def.Epos[j];       
		x=Def.xpos[j];       
        DefOcc=Def.FlagOccupation[j];
        EvTail=round(BandEdge(x/dx,0,0)*1e3)/1e3;;
		DefectDepth=EvTail-Ef; // Defect Depth        
		break;

		case 1: // EMISSION Defect --> Electrode
		Ei=i*EnergyStep; // Energy of Carrier in a defect state
		Ef=j*EnergyStep; // Energy of final state	
		x=initialx;
        EvTail=round(BandEdge(x/dx,0,0)*1e3)/1e3;
		DefectDepth=EvTail-Ei; // Defect Depth
		break;

		case 2: // Phonon Assisted Defect to Defect
		Ei=i*EnergyStep;
		Ef=Def.Epos[j];
		x=Def.xpos[j];
        EvTail=BandEdge(x/dx,0,0);
		DefectDepth=EvTail-Ef; // Defect Depth		
		break;

		case 3: // Electrode to Interface Defect
		Ei=i*EnergyStep;
		Ef=Def.ID_Epos[j];		
		x=Def.ID_xpos[j];		
        EvTail=BandEdge(x/dx,0,0);
		DefectDepth=EvTail-Ef; // Defect Depth		
		break;
	}	
	if (DefectDepth<0 && abs(DefectDepth<1E-4))
	{
		raise(SIGSEGV);
	}	
/**************************************************************************************************
					CALCULATING THE COUPLED STATE USING PHONON MODES
					Choice 0 : Inelastic Transition using phonon mode
					       1 : No Transition available
					       2 : Elastic Transition using phonon mode
**************************************************************************************************/
	switch(t)
    {
        case 0: // Injection from extended states -> localized states
        if (Ef>Ei)
        {
            for (p=-5;p<=-1;p++)
            {                
			    Difference=(Ef+p*homega)-Ei;
                if (abs(Difference)<=0.005 && Difference>=0)
                {
                    Ep=Ef+p*homega;
                    PhononMode=p;
                    Choice=0;                    
                    break;
                }
                else
                    Choice=1;
            }
        }
        else if (Ef<Ei)
        {
            for (p=1;p<5;p++)
            {			    
			    Difference=(Ef+p*homega)-Ei;
                if (abs(Difference)<=0.005)
                {
                    Ep=Ef+p*homega;
                    PhononMode=p;
                    Choice=0;
                    break;
                }
                else
                    Choice=1;
            }
        }
        else if (abs(Ef-Ei)<=0.005)
            Choice=2;

        break;

        case 1: // Emission from localized states -> extended states
        if (Ef>Ei)
        {
            for (p=1;p<=5;p++)
            {                
			    Difference=(Ef-p*homega)-Ei;
                if (abs(Difference)<=0.005 && Difference>=0)
                {
                    Ep=Ef-p*homega;
                    PhononMode=p;
                    Choice=0;                    
                    break;
                }
                else
                    Choice=1;
            }           
        }
        else if (Ef<Ei)
        {
            for (p=-5;p<=-1;p++)
            {			    
			    Difference=(Ef-p*homega)-Ei;                
                if (abs(Difference)<=0.005 && Difference>=0)
                {
                    Ep=Ef-p*homega;
                    PhononMode=p;
                    Choice=0;
                    break;
                }
                else
                    Choice=1;
            }
        }
        else if (abs(Ef-Ei)<=0.005)
            Choice=2;

        break;
    }
/**************************************************************************************************
					CALCULATING THE COUPLED STATE USING PHONON MODES
**************************************************************************************************/
    if (t==0 || t==1) // INJECTION/EMISSION
    {
	    switch (Choice)
	    {
		    case 0: // Inelastic Transition R=N(E)f(E)T(E)c(E,x)dE
			R=InelasticTransition(Ep,DefectDepth,x,PhononMode,t);
		    break;

		    case 1: // No Transition Available R=0
		    R=0;            
		    break;

		    case 2: // Elastic Transition R=C(E)f(E)T(E)
		    R=ElasticTransition(Ei,DefectDepth,x,t);
		    break;
	    }
    }
    else if (t==2) // Defect to Defect Transitions
    {
	    switch (Choice)
	    {
		    case 0:
    		R=PhD2D(i,j,PhononMode,Ep,Def);
	    	break;

		    case 1:
    		R=0;
	    	break;

		    case 2:		
    		R=ElasticTransition(Ei,DefectDepth,x,t);
	    	//R=PhD2D(i,j,PhononMode,Ep,Def);
		    break;
    	}	
    }

    //cout<<"Ei : "<<Ei<<"   "<<"Ef : "<<Ef<<"   "<<"p : "<<p<<"    "<<"Ep : "<<Ep<<"   "<<"x : "<<x<<"   "<<"R : "<<R<<'\n';
    
	return R;
}

void NormalizeRCTable(RateTable& TRT,Defect& Def)
{
	int i,j,dim,count=0;;
	double temp,*max;

	dim=VBOffset/EnergyStep;
	max=new double[dim+1];

	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		for (j=1;j<nDefect;j++)
		{
			TRT.RCTable[i][j]=TRT.RCTable[i][j-1]+TRT.RCTable[i][j];
			temp=TRT.RCTable[i][j];
		}
		max[i]=temp;
		TRT.GammaRC[i]=temp;
		temp=0;
        count++;
	}

	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		for (j=0;j<nDefect;j++)
		{
			TRT.NRCTable[i][j]=TRT.RCTable[i][j]/max[i];	
		}
	}	
	delete [] max;    
}

void NormalizeIRTable(RateTable& TRT,Defect& Def)
{
	int i,j,k,dim;
	double *Max;

	dim=VBOffset/EnergyStep;
	Max=new double[dim+1];

	for (i=1;i<2;i++)
	{		
		for (j=0;j<=VBOffset/EnergyStep;j++)
			Max[j]=0;

		switch (i)
		{
			case 0: // Multiphonon Transition

			for (j=0;j<=VBOffset/EnergyStep;j++)
			{
				for (k=1;k<IntTotal;k++)
				{
					TRT.InterfaceMech[j][k][i]=TRT.InterfaceMech[j][k-1][i]+TRT.InterfaceMech[j][k][i];					
				}
				Max[j]=TRT.InterfaceMech[j][IntTotal-1][i];
			}

			for (j=0;j<=VBOffset/EnergyStep;j++)
			{
				for (k=0;k<IntTotal;k++)
					TRT.InterfaceMech[j][k][i]=TRT.InterfaceMech[j][k][i]/Max[j];
			}

			break;

			case 1: // Hole Capture

			for (j=0;j<=VBOffset/EnergyStep;j++)
			{
				for (k=1;k<IntTotal;k++)
				{
					TRT.InterfaceMech[j][k][i]=TRT.InterfaceMech[j][k-1][i]+TRT.InterfaceMech[j][k][i];					
				}
				Max[j]=TRT.InterfaceMech[j][IntTotal-1][i];
			}

			for (j=0;j<=VBOffset/EnergyStep*0;j++)
			{
				for (k=0;k<IntTotal;k++)
				{
					TRT.InterfaceMech[j][k][i]=TRT.InterfaceMech[j][k][i]/Max[j];					
				}
			}

			break;
		}
	}

	delete [] Max;
}

void iROCSum(RateTable& TRT)
{
	int i,j;
	double temp=0,sum=0;

	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		sum=0;
		for (j=0;j<nDefect;j++)
		{
			sum=sum+TRT.RCTable[i][j];
		}
		TRT.RCSum[i]=sum;
	}
}

double EfSum(Defect& Def,RateTable& TRT,int k)
{
	double sum=0,temp=0,VBmin,VBmax;	
	int j=0,count=0,max;

	switch (k)
	{
		case 1: // Defect to Electrode
        VBmin=BandEdge(Device_X/dx,0,0);
        VBmax=BandEdge(0,0,0);
		max=((VBmax-VBmin)/EnergyStep)+2;
		break;

		case 2: // Defect to Defect
		max=nDefect;		
		break;
	}

    /**************************************************************************
    Efsum calculated the sum of the array which is the 'Gamma Max' of the given
    mechanism
    **************************************************************************/	

	for (j=0;j<max;j++)
	{
		if (k==1)
			temp=TRT.D2E[j];
		else if (k==2)
			temp=TRT.D2D[j];

		sum=sum+temp;
	}
	return sum;
}

double IRSum(RateTable& TRT,int k,int i)
{
	double sum=0,temp=0;
	int j;

	// i - is the initial carrier energy

	for (j=0;j<IntTotal;j++)
	{
		switch (k)
		{
			case 0: // MPROC
			temp=TRT.InterfaceMech[i][j][0];
			break;

			case 1: // Hole Capture
			temp=TRT.InterfaceMech[i][j][1];
			break;

			case 2: // Hole Emission
			temp=TRT.InterfaceMech[i][j][2];
			break;

			case 3: // Electron Capture
			temp=TRT.InterfaceMech[i][j][3];
			break;

			case 4: // Electron Emission
			temp=TRT.InterfaceMech[i][j][4];
			break;
		}
		sum=sum+temp;
	}
	return sum;
}

void NormalizeS1Table(RateTable& TRT,int MG)
{
	int i,Max;
	double tau=0;

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

	// Summation of array
	for (i=1;i<=Max;i++)
	{
		TRT.S1Table[i]=TRT.S1Table[i-1]+TRT.S1Table[i];
	}

	// Calculating the max	
	tau=TRT.S1Table[Max];
	GammaMax=tau;

	// Normalizing the table
	for (i=0;i<=Max;i++)
	{
		TRT.S1Table[i]=TRT.S1Table[i]/tau;     
	}
	
}

void NormalizeS2Table(RateTable& TRT,int k)
{
	int j,end;
	double max=0,VBmin,VBmax;

	switch (k)
	{
		case 1: //Defect to Electrode Emission        
        VBmin=BandEdge(Device_X/dx,0,0);
        VBmax=BandEdge(0,0,0);
		end=((VBmax-VBmin)/EnergyStep)+2;

		//Summation
		for (j=1;j<end;j++)
		{
			TRT.D2E[j]=TRT.D2E[j-1]+TRT.D2E[j];   
		}
		max=TRT.D2E[end-1];        

		//Normalization
		for (j=0;j<end;j++)
        {
			TRT.D2E[j]=TRT.D2E[j]/max;            
        }

		break;

		case 2: //Defect to Defect Emission
		//Summation
		for (j=1;j<nDefect;j++)
		{
			TRT.D2D[j]=TRT.D2D[j-1]+TRT.D2D[j];
		}
		max=TRT.D2D[nDefect-1];

		//Normalization
		for (j=1;j<nDefect;j++)
			TRT.D2D[j]=TRT.D2D[j]/max;

		break;
	}
}

int ValidateFinalState(double Ein,double Ef)
{
    int Flag=0,p;
	double Difference;

	if (abs(Ein-Ef)<0.005)
	{
		p=0;
		Flag=1;
	}
	else if (Ein>Ef)
	{
		for (p=-5;p<=-1;p++)
		{
			Difference=(Ein+p*homega)-Ef;
			if (abs(Difference)<=0.005)
			{
				Flag=1;
				break;
			}
			else
				Flag=0;
		}
	}
	else if (Ein<Ef)
	{
		for (p=1;p<=5;p++)
		{
			Difference=(Ein+p*homega)-Ef;
			if (abs(Difference)<=0.005)
			{
				Flag=1;
				break;
			}
			else
				Flag=0;
		}
	}	

	if (Flag==1)
	{
    	switch (abs(p))
	    {
		    case 0:		
    		ph0++;
	    	break;

    		case 1:
	    	ph1++;
		    break;

    		case 2:
		    ph2++;
	    	break;

    		case 3:
	    	ph3++;
		    break;

    		case 4:
	    	ph4++;
		    break;

    		case 5:
	    	ph5++;
		    break;
    	}
	}
	return Flag;
}

int BinarySearch(RateTable& TRT,int loci,double r,int s)
{
    double VBmin,VBmax;
	int state,flag=0,first,last,middle,t,count=0;
	
	if (s==0 || s==2)
	{
		first=0;
		last=nDefect;
		middle=(first+last)/2;
	}
	else if (s==1)
	{
        VBmin=BandEdge(Device_X/dx,0,0);
        VBmax=BandEdge(0,0,0);		

		first=0;
		last=((VBmax-VBmin)/EnergyStep)+2;
		middle=(first+last)/2;        
	}	

	/**************************************
	0 - Initial ROC
	1 - Defect to Electrode Emission
	2 - Defect to Defect Emission
	**************************************/

	while (flag==0)
	{
		count++;
		switch (s)
		{
			case 0: // Initial ROC, max - nDefect

			if (r>TRT.NRCTable[loci][middle])
				first=middle;
			else if (r<=TRT.NRCTable[loci][middle])
				last=middle;

			t=last-first;

			break;

			case 1: // Defect to Electrode Emission, max - VBmax

			if (r>TRT.D2E[middle])
				first=middle;
			else if (r<=TRT.D2E[middle])
				last=middle;

			t=last-first;			

			break;

			case 2: //  Defect to Defect Emission, max - nDefect

			if (r>TRT.D2D[middle])
				first=middle;
			else if (r<=TRT.D2D[middle])
				last=middle;

			t=last-first;				
			break;
		}		

		if (t==1)
			flag=1;

		middle=(first+last)/2;
	}

	state=last;
	return state;
}
