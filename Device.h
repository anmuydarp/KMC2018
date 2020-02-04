#ifndef DEVICEVAR_H
#define DEVICEVAR_H

using namespace std;

class Device
{ 
public:
	double *xmesh,*ymesh,*zmesh;
	int DevInit();
	int DevDel();
};

class Defect
{
public:
	double *xpos,*ypos,*zpos,*Epos,*CoupledState,*CSRate,*ID_Epos,*ID_xpos;
	int *SliceDefect,*FlagOccupation,*MGFLAG,*MGFLAG1D,*MGHopCheck,*ID_Occupation;
	/***********************************************************************************************
	FlagOccupation Definition - 0 - Occupied by an electron (hopping is permitted to this state)
	                          - 1 - Occupied by an hole (hopping is NOT permitter to this state)
	                          - 2 - Unassigned state 
	***********************************************************************************************/	
	int DefectInit();
	int DefectDel();
	//Rate Capture Table
	//Transition Rate Table (Defect-Defect, Poole-Frenkel Emission, Tunneling Emission)
};

class CarrierH
{
public:
	double *Energy,*xpos,*ypos,*zpos;
	double **TauPath,**xPath,**EPath,*TotalTime,*FinalRate,*Epen,*tauI_n0,*tauI_p0,*tauHEm,*tauECap;
	int *Status,*TransType,*TotPath;
	/********************************************
	Status Definition - 0 - Incident on barrier
	                  - 1 - Hopping Via Defects
	                  - 2 - Collected 
	********************************************/
	int CarrierInit();
	int CarrierDel();	
};

class RateTable
{
public:
	double **RCTable,*RCSum,*D2D,*PF,*D2E; // RCTable - Rate of Capture Table, D2D - Defect to Defect, PF - Poole Frenkel, D2E - Defect to Electrode Emission, MRT - Main Rate Table
	double *S1Table,**NRCTable; // Stage 1 table to select mechanism
	double *GammaRC; // Gamma Max for iROC
	double **TrackRate;
	double ***InterfaceMech,*IRMech;
	
	int TRTInit();
	int TRTDel();
};
// Transition Rate Table
#endif
