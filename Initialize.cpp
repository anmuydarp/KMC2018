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
#include "Device_Parameters.h"
#include "Physics_Parameters.h"

using namespace std;

#define FUNC(x) ((*func)(x))

double func(double x)
{
    double ExpTerm,EvTail;  
    
    // Similar triangles approach
    //EvTail=VBOffset-(x*VBOffset/(ithickness));

    //Equation of line 
    // y = mx + c [y = Band edge, x = distance, m = slope, c = y-intercept]    
    EvTail=BandEdge(x/dx,0,0);    
  
    //ExpTerm=exp(-(EGlobal-EvTail)/Evt); // This is how it should be for valence bands where EGlobal > EvTail

    ExpTerm=exp(-(EvTail-EGlobal)/Evt); // Using this term right now as I've inverted the energy band, so everything looks like it's a CB

    return ExpTerm;
}

double funcMG(double x)
{
	double D,gamma,P,En,EvTail,Ep,f0,Ef_rel;

    // Similar triangles approach
	//EvTail=-(VBOffset-(x*VBOffset/(ithickness)));
    
    //Equation of line 
    // y = mx + c [y = Band edge, x = distance, m = slope, c = y-intercept]    
    EvTail=BandEdge(x/dx,0,0);    

	Ef_rel=EvTail+Efermi; // Fermi level is Efermi eV above the valence band

	//f0=2*exp((Efermi-EGlobal)/KT)/(1+2*exp((Efermi-EGlobal)/KT)+exp((2*Efermi-2*EGlobal-U)/KT));

	f0=2*exp((Ef_rel-EGlobal)/KT)/(1+2*exp((Ef_rel-EGlobal)/KT)+exp((2*Ef_rel-2*EGlobal-U)/KT));

	Ep=Ef_rel+0.5*delta;

	En=EGlobal+sigma*sigma/(2*Ev0);

	P=(1/(sigma*sqrt(2*pi)))*exp(-pow(En-Ep,2.0)/(2*sigma*sigma));

	gamma=Nvt*pow(H/NSiSi,KT/(4*Ev0))*(2*Ev0*Ev0/(2*Ev0-KT))*exp(-(0.5/Ev0)*(Ep-EvTail-(sigma*sigma/(4*Ev0))));
	MG_gamma=gamma; // Using global variables

	D=gamma*pow((2/f0),0.5*KT/Ev0)*P;

	return D;
	
}

void Initialize_Defects(Defect& Def)
{
	// Calculate total number of defects in each slice. Each slice is x (1 nm) * y (10 nm) * z (10 nm)
	// 1 - Calculate N(E) for each dx in the slice
	// 2 - Sum N(E) for the slice and multiply by dy and dz
	// 3 - Calculate volume of slice
	// 4 - Calculate total number of discrete defects in the device
	// 5 - Spread the defects uniformly in real space in each slice
	// 6 - Spread the defects according to exponential decay function in energy	

	int i,j,k,dim,count;
	double E,X1_Integral,X2_Integral,Ex_Integral,mg_Integral,xmin,xmax,dmax,temp,temp1,total,totaldefects=0.0,Eb,Ea,y;
	double *NEx,*Nmg,**ETemp;

	// X1_Integral - Band tail states , X2_Integral - Mid gap states

	ofstream myfile;

	// SliceDensity - Total number of discrete defects in each slice
	// NE[] - Density of defects along the x axis

    dim=round(VBDelta/EnergyStep);
   
	NEx=new double[dim];
	Nmg=new double[dim];

    ETemp=new double*[xpoint];
    for (i=0;i<xpoint;i++)
    {
        ETemp[i]=new double[dim];
    }

	for (i=0;i<xpoint;i++) // First if - Loop to change slice
	{
        // Calculating the bandedge   
        y=BandEdge(i,0,0);
   
        count=0;
        for (j=y/EnergyStep-VBDelta/EnergyStep;j<y/EnergyStep;j++)	
		{
            /*******************************************************************************
            1. dmax will be calculated by equation of the line
            2. y = mx + c [y = Band edge, x = distance, m = slope, c = y-intercept]
            3. y = -0.031x + 0.470 --> will depend on device operating conditions
            4. x is in nm, hence all x's have to be multiplied by dx
            5. The band minma on the a-Si:H(i)/c-Si interface is always scaled to '0'
            *******************************************************************************/
			
            E=EnergyStep*j;
			EGlobal=E;

            dmax=BandEdge(0,E,1);

			xmin=i*dx;			

			if (dmax>=(i+1)*dx)
				xmax=(i+1)*dx;
			else 
				xmax=dmax;			

            //Intergration of DOS w.r.t x
			if (xmax<xmin)
				X1_Integral=0.0;
			else
			{
				X1_Integral=qsimp(&func,xmin,xmax);
				X2_Integral=qsimp(&funcMG,xmin,xmax);				
			}

            ETemp[i][count]=E;
			NEx[count]=Nvt*X1_Integral*1E2;
			Nmg[count]=X2_Integral*1E2;

            count++;
		}

	    // Integration of DOS along E axis
		temp=0.0;
        temp1=0.0;

        for (k=0;k<count-1;k++)
		{
		    Eb=ETemp[i][k+1];
            Ea=ETemp[i][k];	
		
			Ex_Integral=(Eb-Ea)*(NEx[k]+NEx[(k+1)])*0.5; // Integral using Trapezoidal rule - Localized States
			mg_Integral=(Eb-Ea)*(Nmg[k]+Nmg[(k+1)])*0.5; // Integral using Trapezoidal rule - Mid Gap States

        	temp=temp+Ex_Integral;
			temp1=temp1+mg_Integral;
		}
		Ex_Integral=temp;
		mg_Integral=temp1;	      

		// Integration of DOS along y and z axis
		total=(Ex_Integral+mg_Integral)*Device_Y*Device_Z*1E4;

		Def.MGFLAG[i]=Ex_Integral*Device_Y*Device_Z*1E4;

		// Total Number of defects in each slice
	    Def.SliceDefect[i]=int(total);

		//cout<<"Total Number of Defects in Slice : "<<i<<"   "<<Def.SliceDefect[i]<<"    "<<"xmin : "<<xmin<<"   "<<"xmax : "<<xmax<<'\n';
		cout<<"Slice : "<<i<<"  "<<"Ex States : "<<Ex_Integral*Device_Y*Device_Z*1E4<<"   "<<"Mg States : "<<mg_Integral*Device_Y*Device_Z*1E4<<'\n';

        TotalLocalDefects=TotalLocalDefects+(Ex_Integral*Device_Y*Device_Z*1e4);
        TotalMGDefects=TotalMGDefects+(mg_Integral*Device_Y*Device_Z*1e4);
        
		totaldefects=totaldefects+Def.SliceDefect[i];
	}

	cout<<"Total Number of  Defects : "<<totaldefects<<'\n';

	delete [] NEx;
	delete [] Nmg;
    delete [] ETemp;
    
	// Calculating total number of interface defects
	IntTotal=Nint*Device_Y*Device_Z*1E4;

	cout<<'\n'<<"Interface Defects : "<<IntTotal<<'\n';	

	Initialize_RealSpace_Defects(Def);

	Initialize_Energy_Defects(Def);
    
	//Initilize_IntDef_Energy(Def);

	//Initilize_IntDef_Occupation(Def);	

	//Initialize_Occupation(Def);
	Initialize_Occupation_2(Def);

	cout<<'\n'<<"=========================================INITIALIZATION OF DEFECTS COMPLETE==========================================="<<'\n';    
}

void Initialize_Carriers(CarrierH& Carr)
{
	double Estep=1E-3,sigma=0.07,mu=0.150;

	int i,j,dim;
	double *PDF,sum=0,r;

	ofstream myfile;

	dim=int(VBOffset/Estep)+1;
	PDF=new double[dim];
    
    /*
    // Creating PDF Table
	for (i=0;i<=VBOffset/Estep;i++)
	{
		PDF[i]=(1/(sigma*sqrt(2*pi)))*exp(-pow(i*Estep-mu,2.0)/(2*sigma*sigma));
		sum=sum+PDF[i];
	}

	// Summation and Normalization PDF Table
	PDF[0]=PDF[0]/sum;
	for (i=1;i<=VBOffset/Estep;i++)
	{
		PDF[i]=PDF[i]/sum+PDF[i-1];		
	}

	// Choosing Energies of carreirs according a Gaussian Distribution
	//myfile.open("IMAPSCarrierEnergy.txt");
	for (i=0;i<CarrierMax;i++)
	{		
		r=(double)rand()/RAND_MAX;

		for (j=0;j<VBOffset/Estep;j++)
		{
			if (r>PDF[j] && r<=PDF[j+1])
			{
				Carr.Energy[i]=(j+1)*Estep;
				Carr.Status[i]=0; // All carriers are initialized at the barrier
				Carr.xpos[i]=0; // Initialized at the barrier
				//myfile<<Carr.Energy[i]<<'\n';				
			}
		}
	}
	//myfile.close();

	/************************************************************************
	MAXWELLIAN INITIALIZATION
	************************************************************************/
	for (i=0;i<CarrierMax;i++)
	{
		r=(double)rand()/RAND_MAX;
		Carr.Energy[i]=-1.5*(kb*300/q)*log(r);
		Carr.Status[i]=0;
		Carr.xpos[i]=0; 
	}
    
	//myfile.close();*/
	delete [] PDF;

	cout<<'\n'<<"=========================================INITIALIZATION OF CARRIERS COMPLETE=========================================="<<'\n';
}
