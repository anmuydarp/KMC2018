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

void CalculateN1P1(Defect& Def)
{
	

}

void Calculate_SRH_Interface_Lifetime(Defect& Def,CarrierH& Carr,int k)
{
	int i;
	double n1,p1,temp_n0,temp_p0,Ec=1.12,Ev=0.0,SRH_n,SRH_p;
	double del_n,del_p,n0,p0,n,p,avg_taun0=0,avg_taup0=0;

	switch (k)
	{
		case 0: // Nt = 1E10		
		n0=100;
		p0=4E18;
		n=4E14;
		p=4E18;	
		del_n=n-n0;
		del_p=p-p0;
		break;

		case 1: // Nt = 1E12
		break;
	}	

	for (i=0;i<IntTotal;i++)
	{
		n1=Nc_Si*exp((Def.ID_Epos[i]*(-1)-Ec)/Vt);
		p1=Nv_Si*exp((Ev-Def.ID_Epos[i]*(-1))/Vt);

		//if (Def.ID_Occupation[i]==0) // Occupied by electron
			temp_n0=temp_n0+n1;
		//else if (Def.ID_Occupation[i]==1) // Occupied by a hole
			temp_p0=temp_p0+p1;		
	}

	temp_n0=temp_n0/IntTotal;
	temp_p0=temp_p0/IntTotal;

	for (i=0;i<IRcount;i++)
	{
		avg_taun0=avg_taun0+Carr.tauI_n0[i];
		avg_taup0=avg_taup0+Carr.tauI_p0[i];
	}

	avg_taun0=avg_taun0/IRcount;
	avg_taup0=avg_taup0/IRcount;

	cout<<"tauI_n0 : "<<avg_taun0<<"   "<<"tauI_p0 : "<<avg_taup0<<'\n';

	// Electrons in p-type material
	SRH_n=(del_n*avg_taup0*(n+temp_n0)+del_n*avg_taun0*(p+temp_p0))/(n*p-ni_Si*ni_Si);

	// Electrons in n-type material
	SRH_p=(del_p*avg_taup0*(n+temp_n0)+del_p*avg_taun0*(p+temp_p0))/(n*p-ni_Si*ni_Si);

	cout<<'\n';
	cout<<"n SRH : "<<SRH_n<<"    "<<"p SRH : "<<SRH_p<<'\n';
	cout<<"n1 : "<<temp_n0<<"   "<<"p1 : "<<temp_p0<<'\n';
}

void Calcualte_SRH_Bulk_Lifetime(Defect& Def)
{

}
