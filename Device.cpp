#include <stdio.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <csignal>
#include "Device.h"
#include "Global.h"
#include "Device_Parameters.h"

using namespace std;

int Defect::DefectInit(void)
{
    xpos=new double[40000];
	ypos=new double[40000];
	zpos=new double[40000];
	Epos=new double[40000];
	FlagOccupation=new int[40000];
	ID_Epos=new double[100];
	ID_xpos=new double[100];
	ID_Occupation=new int[100];

	CoupledState=new double[40000];
	CSRate=new double[40000];

	SliceDefect=new int[xpoint];
	MGFLAG=new int[xpoint];
	MGFLAG1D=new int[40000];
	MGHopCheck=new int[50];

	return 0;
}

int Defect::DefectDel(void)
{
	delete [] xpos;
	delete [] ypos;
	delete [] zpos;
    delete [] Epos;
	delete [] FlagOccupation;
	delete [] ID_Epos;
	delete [] ID_xpos;
	delete [] ID_Occupation;
	delete [] CoupledState;
	delete [] CSRate;
	delete [] SliceDefect;
	delete [] MGFLAG;
	delete [] MGFLAG1D;
	delete [] MGHopCheck;

	return 0;
}

int CarrierH::CarrierInit(void)
{
    int total=2000;

	Energy=new double[CarrierMax];
	xpos=new double[CarrierMax];
	ypos=new double[CarrierMax];
	zpos=new double[CarrierMax];
	Status=new int[CarrierMax];

	int i;
	xPath=new double*[CarrierMax];
	for (i=0;i<CarrierMax;i++)
		xPath[i]=new double[total];

	EPath=new double*[CarrierMax];
	for (i=0;i<CarrierMax;i++)
		EPath[i]=new double[total];

	TauPath=new double*[CarrierMax];
	for (i=0;i<CarrierMax;i++)
		TauPath[i]=new double[total];

	TotalTime=new double[CarrierMax];

	FinalRate=new double[CarrierMax];
	TransType=new int[CarrierMax];
	TotPath=new int[CarrierMax];

	Epen=new double[CarrierMax];
	tauI_n0=new double[CarrierMax];
	tauI_p0=new double[CarrierMax];
	tauHEm=new double[CarrierMax];
	tauECap=new double[CarrierMax];

	return 0;
}

int CarrierH::CarrierDel(void)
{
	delete [] Energy;
	delete [] xpos;
	delete [] ypos;
	delete [] zpos;
	delete [] Status;	
	delete [] xPath;
	delete [] EPath;
	delete [] TauPath;
	delete [] TotalTime;
	delete [] FinalRate;
	delete [] TransType;
	delete [] TotPath;	
    delete [] Epen;
	delete [] tauI_p0;
	delete [] tauI_n0;
	delete [] tauHEm;
	delete [] tauECap;

	return 0;
}

int RateTable::TRTInit(void)
{
	int i,j,dim1,dim2,dim3;
    double VBmin,VBmax;

    /************************************************************************************************************************
    This is assuming that the band-edge decrease from the a-Si:H(i)/c-Si heterointerface to the a-Si:H(i)/a-Si:H(p) interface
    i.e. VBmax > VBmin
    ************************************************************************************************************************/

    VBmin=BandEdge(Device_X/dx,0,0);
    VBmax=BandEdge(0,0,0);

    dim1=(VBOffset/EnergyStep)+1;
	dim2=(VBmax-VBmin)/EnergyStep;
	dim3=sigma_Int/EnergyStep;
   
	RCTable=new double*[dim1];
	for (i=0;i<dim1;i++)
		RCTable[i]=new double[nDefect];

	RCSum=new double[dim1];

	NRCTable=new double*[dim1];
	for (i=0;i<dim1;i++)
		NRCTable[i]=new double[nDefect];

	D2D=new double[nDefect];
	PF=new double[2];
	D2E=new double[dim2+1];    

	S1Table=new double[8];

	GammaRC=new double[dim1];

	TrackRate=new double*[3];
	for (i=0;i<3;i++)
		TrackRate[i]=new double[10000];

	InterfaceMech=new double**[dim1+1];
	for (i=0;i<=VBOffset/EnergyStep;i++)
	{
		InterfaceMech[i]=new double*[dim3+2];
		for (j=0;j<sigma_Int/EnergyStep;j++)
			InterfaceMech[i][j]=new double[5];
	}

	IRMech=new double[8];
   
	return 0;
}

int RateTable::TRTDel(void)
{
	delete [] RCTable;
    delete [] RCSum;
	delete [] NRCTable;
	delete [] D2D;
	delete [] PF;
	delete [] D2E;
	delete [] S1Table;
	delete [] GammaRC;
	delete [] TrackRate;
    delete [] InterfaceMech;
	delete [] IRMech;

	return 0;
}

int Device::DevInit(void)
{
	return 0;
}

int Device::DevDel(void)
{
	return 0;
}
