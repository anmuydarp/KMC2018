#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <string.h>
#include <cstdlib>
#include <stdlib.h>
#include <csignal>
#include "Global.h"
#include "Device_Parameters.h"
#include "Physics_Parameters.h"

using namespace std;

double BandEdge (double x,double E,int i)
{
    double y;

    switch (i)
    {
        case 0: // Band-Edge
        y=BandSlope*x+BandIntercept;
        break;

        case 1: // x
        y=dx*(E-BandIntercept)/BandSlope;
        break;
    }

    return y;
}

double QFL (double x)
{
    double y;

    y=QFLSlope*x+QFLIntercept;

    return y;
}

void Initialize_RealSpace_Defects(Defect& Def)
{
	int i,j,count=0;
	double r;

	for (i=0;i<xpoint;i++) // Loop to change slice
	{		
		for (j=0;j<Def.SliceDefect[i];j++) // Loop to assign position (x,y,z) for each defect
		{
			r=(double)rand()/RAND_MAX;			
			Def.xpos[count]=(i+r)*dx;
			Def.ypos[count]=Device_Y*r;
			Def.zpos[count]=Device_Z*r;
			//cout<<"i : "<<i<<"  "<<"j : "<<j<<"  "<<"x pos : "<<Def.xpos[count]/dx<<'\n';
			count++;
		} 
		//cout<<"Slice Number : "<<i<<"   "<<"Number of defects in slice : "<<Def.SliceDefect[i]<<'\n';
	}	
	cout<<"Total Count : "<<count<<'\n';
	nDefect=count;

	for (i=0;i<IntTotal;i++)
	{
		Def.ID_xpos[i]=0;
	}
		
	cout<<'\n'<<"=========================================INITIALIZED REAL SPACE DEFECTS==========================================="<<'\n';
	
}

void Initilize_IntDef_Energy(Defect& Def)
{
	double r,temp=0,sigma1=0.1,sigma2=0.08,mu1=0.2,mu2=0.3,Estep=1E-3,Nint_Volume,max;
	int i,j,dim,aEi=0,bEi=0;

	double *GD1,*GD2,*GDTotal;

	dim=Eg_Si/Estep;

	GD1=new double[dim];
	GD2=new double[dim];
	GDTotal=new double[dim];

	ofstream myfile;

	Nint_Volume=Nint/(1E-7); // converting interface density from cm^-2 to cm^-3

	//myfile.open("ID_Epos.txt");
	/*for (i=0;i<IntTotal;i++)
	{
		r=(double)rand()/RAND_MAX;
		Def.ID_Epos[i]=Int_Mid+sigma_Int*r; //-0.41 is the mid gap in the upside down valence band picture
		temp=temp+Def.ID_Epos[i];
		myfile<<Def.ID_Epos[i]<<'\n';
	}*/
	//myfile.close();

	//myfile.open("DoubleGaussian_FD1.txt");
	for (i=0;i<=Eg_Si/Estep;i++)
	{		
		GD1[i]=(1/(sigma1*sqrt(2*pi)))*exp(-pow(i*Estep-mu1,2.0)/(2*sigma1*sigma1));
		GD2[i]=(1/(sigma2*sqrt(2*pi)))*exp(-pow(i*Estep-mu2,2.0)/(2*sigma2*sigma2));
		GDTotal[i]=GD1[i]+GD2[i];
		//myfile<<FD1[i]<<'\n';
	}
	//myfile.close();

	// Summation
	for (i=1;i<=Eg_Si/Estep;i++)
	{
		GDTotal[i]=GDTotal[i-1]+GDTotal[i];
	}
	max=GDTotal[int(Eg_Si/Estep)];

	// Normalization
	for (i=0;i<=Eg_Si/Estep;i++)
	{
		GDTotal[i]=GDTotal[i]/max;
	}

	// Assigning Energy to Interface Defects
	//myfile.open("DoubleGaussian_Energies.txt");
	double en,ep;
	for (j=0;j<IntTotal;j++)
	{
		r=(double)rand()/RAND_MAX;
		for (i=1;i<=Eg_Si/Estep;i++)
		{
			if (r>GDTotal[i-1] && r<=GDTotal[i])
			{
				Def.ID_Epos[j]=i*Estep*(-1);				

				if (Def.ID_Epos[j]*(-1)>0.56)
					aEi++;
				else if (Def.ID_Epos[j]*(-1)<0.56)
					bEi++;

			}
		}
	}
	
	cout<<"aEi : "<<aEi<<"   "<<"bEi : "<<bEi<<'\n';
	//raise(SIGSEGV);
	//myfile.close();

	//cout<<'\n'<<"Average Int Def Energy : "<<temp/IntTotal<<'\n';

	delete [] GD1;
	delete [] GD2;
	delete [] GDTotal;	
}

void Initilize_IntDef_Occupation(Defect& Def)
{
	int i,count=0,Hole=0,Elec=0;
	double FD,FDMax=0,r,Interface_HoleConc,Interface_ElecConc,nbar,pbar,en,ep,Ec,Ev;

	ofstream myfile;

	// Initializing Interface Parameters
	Interface_HoleConc=1E19; // Will change with the interface defect density
	Interface_ElecConc=1E12;
	Ec=1.12;
	Ev=0;

	nbar=Vth*sigma_plus*Interface_ElecConc;
	pbar=Vth*sigma_plus*Interface_HoleConc;

	for (i=0;i<IntTotal;i++)
	{
		Def.ID_Occupation[i]=2; // Unassigned State
	}

	/*for (i=0;i<IntTotal;i++)
	{
		en=Vth*sigma_plus*Nc_aSi*exp((Def.ID_Epos[i]*(-1)-Ec)/Vt);
		ep=Vth*sigma_plus*Nv_aSi*exp((Ev-Def.ID_Epos[i]*(-1))/Vt);

		FD=(nbar+ep)/(en+nbar+pbar+ep);

		//cout<<"E : "<<Def.ID_Epos[i]*(-1)<<"   "<<"en : "<<en<<"   "<<"ep : "<<ep<<"   "<<"FD : "<<FD<<'\n';

		r=(double)rand()/RAND_MAX;

		if (r<FD)
		{
			Def.ID_Occupation[i]=0; // Occupied by electron
			Elec++;
		}
		else if (r>=FD)
		{
			Def.ID_Occupation[i]=1; // Occupied by hole
			Hole++;
		}
	}

	cout<<"Elec : "<<Elec<<"   "<<"Hole : "<<Hole<<'\n';*/

	//raise(SIGSEGV);

	//for (i=0;i<IntTotal;i++)
	//{
		//cout<<"Def.ID_Epos : "<<Def.ID_Epos[i]*(-1)<<"    "<<"Occ : "<<Def.ID_Occupation[i]<<'\n';
	//}

	/*myfile.open("Interface_DefectEnergy_1E13.txt");
	for (i=0;i<IntTotal;i++)
	{
		myfile<<Def.ID_Epos[i]*(-1)<<'\n';
	}
	myfile.close();

	myfile.open("Interface_DefectOcc_1E13.txt");
	for (i=0;i<IntTotal;i++)
	{
		myfile<<Def.ID_Occupation[i]<<'\n';
	}
	myfile.close();*/

	Elec=0;
	Hole=0;
	for (i=0;i<IntTotal;i++)
	{
		r=(double)rand()/RAND_MAX;

		if (r>0.5)
		{
			Def.ID_Occupation[i]=0; // Occupied by electron
			Elec++;
		}
		else 
		{
			Def.ID_Occupation[i]=1; // Occupied by hole
			Hole++;
		}
	}

	// Calculating interface defect densities in a 1 nm slice
	IntElec=Elec/(Device_Y*Device_Z*1E-9*1E6); 
	IntHole=Hole/(Device_Y*Device_Z*1E-9*1E6);

	cout<<'\n'<<"Occupation Test : Hole : "<<Hole<<"   "<<"Electron : "<<Elec<<'\n';
}

void Initialize_Energy_Defects(Defect& Def)
{
	int i,j,k,dimE,dimX,dimMG,count=0;
	double xstep=0.1,EvTail,dmax,E,r,En,Ep,P,f0,Ef_rel;
	double ***DOSTable,***DOSTableMG,**DOSMax,temp=0.0;

	ofstream myfile;

	dimE=round(VBDelta/EnergyStep);
	dimMG=(VBOffset+Eg_aSi)/EnergyStep;
	dimX=1/xstep;
    
   	DOSTable=new double**[xpoint];
	for (i=0;i<xpoint;i++)
	{
		DOSTable[i]=new double*[dimX];
		for (j=0;j<dimX;j++)
			DOSTable[i][j]=new double[dimE];
	}

	DOSTableMG=new double**[xpoint];
	for (i=0;i<xpoint;i++)
	{
		DOSTableMG[i]=new double*[dimX];
		for (j=0;j<dimX;j++)
			DOSTableMG[i][j]=new double[dimMG];
	}

	DOSMax=new double*[xpoint];
	for (i=0;i<xpoint;i++)
		DOSMax[i]=new double[dimX];

/*********************************************************************************
           				CREATING DENSITY OF STATES TABLE (Localized States)
*********************************************************************************/
	for (i=0;i<xpoint;i++) // Loop to change slice, different slices will have different Evtails (i - Slice, j - Slice Interval, k - Energy)
	{
		for (j=0;j<10;j++) // Slice Mesh to calculate different EvTails
		{
			dmax=(i+j*xstep)*dx;

            EvTail=BandEdge(dmax/dx,0,0);

            count=0;
			for (k=(EvTail-VBDelta)/EnergyStep;k<=EvTail/EnergyStep;k++)
			{				
				E=k*EnergyStep;
				if (E>EvTail)
                {
					DOSTable[i][j][count]=0;
                }
				else
                {			
					DOSTable[i][j][count]=Nvt*exp(-(EvTail-E)/Evt);
                }                
				//cout<<"E : "<<E<<"  "<<"EvTail : "<<EvTail<<"   "<<"DOS : "<<DOSTable[i][j][count]<<'\n';
                count++;
			}
		}
	}
    count=0;
    
/*********************************************************************************
           				CREATING DENSITY OF STATES TABLE (Mid Gap States)
*********************************************************************************
	for (i=0;i<xpoint;i++) // Loop to change slice, different slices will have different Evtails (i - Slice, j - Slice Interval, k - Energy)
	{
		for (j=0;j<10;j++) // Slice Mesh to calculate different EvTails
		{
			dmax=(i+j*xstep)*dx;
			EvTail=VBOffset-VBOffset*(dmax/ithickness);
			EvTail=-EvTail; // I've treated everything like a conduction band, but that'll needlessly complicate the VB based formulae of MG's
			//cout<<"dmax : "<<dmax/dx<<"   "<<"Evtail : "<<EvTail<<'\n';

			for (k=0;k<=(VBOffset+Eg_aSi)/EnergyStep;k++)
			{
				E=(k*EnergyStep)-VBOffset;

				if (E<EvTail)
				{
					DOSTableMG[i][j][k]=0;
				}
				else
				{
					Ef_rel=EvTail+Efermi; // (Should probably be EvTail - Efermi because EvTail is a pseudo conduction band)
					//f0=2*exp((Efermi-E)/KT)/(1+2*exp((Efermi-E)/KT)+exp((2*Efermi-2*E-U)/KT));
					f0=2*exp((Ef_rel-E)/KT)/(1+2*exp((Ef_rel-E)/KT)+exp((2*Ef_rel-2*E-U)/KT));
					//Ep=Efermi+0.5*delta;
					Ep=Ef_rel+0.5*delta;
					En=E+sigma*sigma/(2*Ev0);
					P=(1/(sigma*sqrt(2*pi)))*exp(-pow(En-Ep,2.0)/(2*sigma*sigma));
					//DOSTableMG[i][j][k]=MG_gamma*pow((2/f0),0.5*Vt/Ev0)*P;
					DOSTableMG[i][j][k]=MG_gamma*pow((2/f0),0.5*KT/Ev0)*P;
					//cout<<"E : "<<E<<"  "<<"k : "<<k<<"  "<<"DOS : "<<DOSTableMG[i][j][k]<<'\n';
				}
			}
		}
	}

/*********************************************************************************
           				NORMALIZING DENSITY OF STATES TABLE (Localized States)
*********************************************************************************/
 	for (i=0;i<xpoint;i++) // Summation along Energy for each slice interval
 	{
 		for (j=0;j<10;j++)
 		{
 			for (k=1;k<=VBDelta/EnergyStep;k++)
 			{
 				DOSTable[i][j][k]=DOSTable[i][j][k-1]+DOSTable[i][j][k]; 				
 			}
 			DOSMax[i][j]=0.0;
 		}
 	}

 	for (i=0;i<xpoint;i++) // Finding max value in the array 		
 	{
 		for (j=0;j<10;j++)
 		{
 			temp=0.0;
 			for (k=0;k<=VBDelta/EnergyStep;k++)
 			{
 				if (DOSTable[i][j][k]>temp)
 					temp=DOSTable[i][j][k];
 			}
 			DOSMax[i][j]=temp;   
 		} 		
 	}
    
 	for (i=0;i<xpoint;i++)
 	{
 		for (j=0;j<10;j++)
 		{
 			for (k=0;k<=VBDelta/EnergyStep;k++)
 			{
 				DOSTable[i][j][k]=DOSTable[i][j][k]/DOSMax[i][j]; 				
 			}
 		}
 	}
    
/*********************************************************************************
           				NORMALIZING DENSITY OF STATES TABLE (Mid Gap States)
*********************************************************************************
 	for (i=0;i<xpoint;i++) // Summation along Energy for each slice interval
 	{
 		for (j=0;j<10;j++)
 		{
 			for (k=1;k<=(VBOffset+Eg_aSi)/EnergyStep;k++)
 			{
 				DOSTableMG[i][j][k]=DOSTableMG[i][j][k-1]+DOSTableMG[i][j][k];
 			}
 			DOSMax[i][j]=0.0;
 		}
 	}

 	for (i=0;i<xpoint;i++) // Finding max value in the array 		
 	{
 		for (j=0;j<10;j++)
 		{
 			temp=0.0;
 			for (k=0;k<=(VBOffset+Eg_aSi)/EnergyStep;k++)
 			{
 				if (DOSTableMG[i][j][k]>temp)
 					temp=DOSTableMG[i][j][k];
 			}
 			DOSMax[i][j]=temp; 		
 		} 		
 	}

 	for (i=0;i<xpoint;i++)
 	{
 		for (j=0;j<10;j++)
 		{
 			for (k=0;k<=(VBOffset+Eg_aSi)/EnergyStep;k++)
 			{
 				DOSTableMG[i][j][k]=DOSTableMG[i][j][k]/DOSMax[i][j]; 				
 			} 			
 		}
 	} 	

 	delete [] DOSMax;

/*********************************************************************************
						ASSIGNING ENERGY TO EACH DEFECT
*********************************************************************************/
	double r1,*MGStates;
    int xmesh,Flag=0,MGcount=0,ExCount=0,MidGap,count1=0;

    count=0;
    MGStates=new double[3000];
	for (i=0;i<xpoint;i++) // Silce
	{		
		//cout<<"Slice Defect in assigning energy routine : "<<Def.SliceDefect[i]<<'\n';

		for (j=0;j<Def.SliceDefect[i];j++)
		{			
			xmesh=int(Def.xpos[count]*10/dx)-i*10;

			Flag=0;			
			while (Flag==0)
			{	
				if (j<Def.MGFLAG[i])
					MidGap=0; // Localized States
				else
					MidGap=1; // Mid Gap States

				switch (MidGap)
				{
					case 0: // Localized States
					
					Def.MGFLAG1D[count]=0;
					r=(double)rand()/RAND_MAX;
    
                    EvTail=BandEdge(Def.xpos[count]/dx,0,0);
                    count1=0;
					for (k=(EvTail-VBDelta)/EnergyStep;k<=EvTail/EnergyStep;k++) // have to change thins logic
					{	                        			
						if (r>DOSTable[i][xmesh][count1] && r<=DOSTable[i][xmesh][count1+1])
						{
                            EvTail=BandEdge(Def.xpos[count]/dx,0,0);
							Def.Epos[count]=k*EnergyStep;
						
							if (Def.Epos[count]>EvTail)
							{
								r1=(double)rand()/RAND_MAX;
								Def.Epos[count]=Def.Epos[count]-r1*(Def.Epos[count]-EvTail);
                                cout<<"Check Condition"<<'\n';
                                raise(SIGSEGV);
							}
							ExCount++;
							Flag=1;
							break;					
						}
						else
						{
							Flag=0;						
						}
                        count1++;
					}

					break;

					case 1: // Mid Gap States
                                      
					/*Def.MGFLAG1D[count]=1;
					r=(double)rand()/RAND_MAX;

					for (k=0;k<(VBOffset+Eg_aSi)/EnergyStep;k++)
					{
						if (r>DOSTableMG[i][xmesh][k] && r<=DOSTableMG[i][xmesh][k+1])
						{
							EvTail=VBOffset-VBOffset*(Def.xpos[count]/ithickness);
							Def.Epos[count]=((k*EnergyStep)-VBOffset)*(-1);

							if (Def.Epos[count]>EvTail)
							{
								r1=(double)rand()/RAND_MAX;	
								Def.Epos[count]=Def.Epos[count]-(1+r1)*(Def.Epos[count]-EvTail);			
							}
							MGStates[MGcount]=Def.Epos[count];
							MGcount++;
							Flag=1;
							break;
						}
						else
						{
							Flag=0;
						}
					}*/

					break;
				}				
			}
			count++;			
		}        
	}    

	cout<<"MGcount Check : "<<MGcount<<"   "<<"ExCount Check : "<<ExCount<<"   "<<"count : "<<count<<'\n';

	delete [] DOSTable;
	delete [] DOSTableMG;
	delete [] MGStates;	 	

	cout<<'\n'<<"=========================================INITIALIZED ENERGY SPACE DEFECTS========================================="<<'\n';
}

void Initialize_Occupation(Defect& Def)
{
	//Work on initializing the occupation and study about the defects
	int i,j,k,count=0,Max,zerocheck=0,zerocount=0;
	int countEx,countMG,var=0;
	double *FD,*F0MG,*FplusMG,*FminusMG,aSi_BandGap=1.72,HQFL=0.0,FDMax=0.0,FDMaxMG=0.0,temp=0.0,temp1=0.0,temp2=0.0,r;

    // Import HQFL

	ofstream myfile;

	FD=new double[TotalLocalDefects];
	F0MG=new double [TotalMGDefects];
	FplusMG=new double[TotalMGDefects];
	FminusMG=new double[TotalMGDefects];

	for (i=0;i<xpoint;i++)
	{
		for (j=0;j<Def.SliceDefect[i];j++)
		{
			Def.FlagOccupation[count]=2; // Unassigned state
			count++;
		}
	}
	Max=count;
	count=0;
	
	/********************************************************************************
	1. Calculating Fermi - Dirac probability of occupation for each defect site.
	2. Finiding FD Max and normalizing Fermi - Dirac table.
	3. 
	********************************************************************************/
    
	for (k=0;k<Max;k++) // Main outer loop
	{
		countEx=0;
		countMG=0;
		FDMaxMG=0;
		FDMax=0;
		zerocount=0;
		//cout<<"Occupation : "<<k<<'\n';

		for (i=0;i<Max;i++)
		{
			if (Def.MGFLAG1D[i]==0) // Localized States
			{				
				if (Def.FlagOccupation[i]==0) // Occupied by an electron
				{
					FD[countEx]=0;					
				}
				else if (Def.FlagOccupation[i]==2 || Def.FlagOccupation[i]==1)
				{
                    HQFL=QFL(Def.xpos[i]/dx);                    
					//FD[countEx]=1/(1+exp((-Def.Epos[i]-HQFL)/Vt));					
					FD[countEx]=1/(1+exp((Def.Epos[i]-HQFL)/Vt));
                    //cout<<"Count Ex : "<<countEx<<"   "<<FD[countEx]<<'\n';
				}
				FDMax=FDMax+FD[countEx];
				countEx++;
			}
			else if (Def.MGFLAG1D[i]==1) // Mid Gap States
			{				
				if (Def.FlagOccupation[i]==0)
				{					
					F0MG[countMG]=0;
					FplusMG[countMG]=1/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt));
					FminusMG[countMG]=exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt)/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt));
				}
				else if (Def.FlagOccupation[i]==3)
				{
					FminusMG[countMG]=0;
					F0MG[countMG]=2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt)); // The minus of minus is to flip CB interpretation to VB
					FplusMG[countMG]=1/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt));
				}
				else if (Def.FlagOccupation[i]==2)
				{
					F0MG[countMG]=2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt)); // The minus of minus is to flip CB interpretation to VB
					FplusMG[countMG]=1/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt));
					FminusMG[countMG]=exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt)/(1+2*exp((Efermi_Dev-(-Def.Epos[i]))/Vt)+exp((2*Efermi_Dev-2*(-Def.Epos[i])-U)/Vt));
				}
				FDMaxMG=FDMaxMG+F0MG[countMG]+FplusMG[countMG]+FminusMG[countMG];				
				countMG++;
			}			
		}

		temp=0;
		for (i=0;i<countEx;i++)
		{
			temp=temp+FD[i]/FDMax;
			FD[i]=temp;
			//cout<<"i : "<<i<<"  "<<"FD : "<<FD[i]<<'\n';
		}        
				
		temp=0;
		temp1=0;
		temp2=0;

		// Normalize
		for (i=0;i<countMG;i++)
		{
			temp=temp+F0MG[i]/FDMaxMG;
			temp1=temp1+(F0MG[i]+FplusMG[i])/FDMaxMG;
 			temp2=temp2+(F0MG[i]+FplusMG[i]+FminusMG[i])/FDMaxMG;
			F0MG[i]=temp;
			FplusMG[i]=temp1;
			FminusMG[i]=temp2;
			//cout<<"i : "<<i<<"   "<<"0 : "<<temp<<"  "<<"+ : "<<temp1<<"   "<<"- : "<<temp2<<'\n';
		}		
		
		countEx=0;
		countMG=0;	
		for (i=0;i<Max-1;i++)
		{
			switch (Def.MGFLAG1D[i])
			{
				case 0:  // Occupancy for localizes states
				
				r=(double)rand()/RAND_MAX;
		
				if (Def.FlagOccupation[i]==0)
				{
					// do nothing					
				}
				else if (Def.FlagOccupation[i]==2 || Def.FlagOccupation[i]==1)
				{
					if (r>FD[countEx] && r<=FD[countEx+1])
						Def.FlagOccupation[i]=0; // Occupied by electron
					else 
						Def.FlagOccupation[i]=1; // Occupied by hole
				}
				countEx++;
				
				break;

				case 1: // Occupany for Mid Gap States
				r=(double)rand()/RAND_MAX;
				zerocheck=Def.FlagOccupation[i];
				
				if (zerocheck==0 || zerocheck==1 || zerocheck==3)
				{
					// do nothing
				}
				else if (zerocheck==2)
				{
					if (r<=F0MG[countMG])
					{
						Def.FlagOccupation[i]=0; // 1 electron						
					}
					else if (r>F0MG[countMG] && r<=FplusMG[countMG])
					{
						Def.FlagOccupation[i]=1; // 1 hole						
					}
					else if (r>FplusMG[countMG] && r<=FminusMG[countMG])
					{
						Def.FlagOccupation[i]=3; // 2 electrons						
					}
				}				
				countMG++;
				break;
			}		
		}
		//cout<<"Max : "<<Max<<"  "<<"countEx : "<<countEx<<'\n';		
	}		

	cout<<'\n'<<"=========================================COMPLETED OCCUPATION OF STATES==========================================="<<'\n';

	zerocheck=0;

	int ExC1=0,ExC2=0,MGC1=0,MGC2=0,MGC3=0,Excheck=0,MGcheck=0;
	double *Ee,*Eh,*ED0,*EDminus,*EDplus;
	double *Xe,*Xh,*XD0,*XDminus,*XDplus;

	Ee=new double[TotalLocalDefects];
	Eh=new double[TotalLocalDefects];
	ED0=new double[TotalMGDefects];
	EDminus=new double[TotalMGDefects];
	EDplus=new double[TotalMGDefects];

	Xe=new double[TotalLocalDefects];
	Xh=new double[TotalLocalDefects];
	XD0=new double[TotalMGDefects];
	XDminus=new double[TotalMGDefects];
	XDplus=new double[TotalMGDefects];
	
	for (i=0;i<Max;i++)
	{
		zerocheck=Def.FlagOccupation[i];

		if (Def.MGFLAG1D[i]==0) // Localized States
		{
			switch (zerocheck)
			{
				case 0: // Electron
				Ee[ExC1]=Def.Epos[i];
				Xe[ExC1]=Def.xpos[i];
				ExC1++;				
				break;

				case 1: // Hole
				Eh[ExC2]=Def.Epos[i];
				Xh[ExC2]=Def.xpos[i];
				ExC2++;
				break;
			}
			Excheck++;
		}
		else if (Def.MGFLAG1D[i]==1) // Mid Gap States
		{
			switch (zerocheck)
			{
				case 0: // 1 Electron
				ED0[MGC1]=Def.Epos[i];
				XD0[MGC1]=Def.xpos[i];
				MGC1++;
				break;

				case 1: // 1 Hole
				EDplus[MGC2]=Def.Epos[i];
				XDplus[MGC2]=Def.xpos[i];
				MGC2++;
				break;

				case 3: // 2 Electrons
				EDminus[MGC3]=Def.Epos[i];
				XDminus[MGC3]=Def.xpos[i];				
				MGC3++;
				break;
			}
			MGcheck++;
		}
	}	
	
	delete [] FD;
	delete [] F0MG;
	delete [] FplusMG;
	delete [] FminusMG;

	int assigncheck=0;
	for (i=0;i<Max;i++)
	{
		if (Def.FlagOccupation[i]==2)
			assigncheck++;
	}

//-------------- E pos ------------------------
	myfile.open("ExEe.txt");
	for (i=0;i<ExC1;i++)
	{
		myfile<<Ee[i]<<'\n';
	}
	myfile.close();

	myfile.open("ExEh.txt");
	for (i=0;i<ExC2;i++)
	{
		myfile<<Eh[i]<<'\n';
	}
	myfile.close();
    /*
	myfile.open("MGED0.txt");
	for (i=0;i<MGC1;i++)
	{
		myfile<<ED0[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGEDplus.txt");
	for (i=0;i<MGC2;i++)
	{
		myfile<<EDplus[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGEDminus.txt");
	for (i=0;i<MGC3;i++)
	{
		myfile<<EDminus[i]<<'\n';
	}
	myfile.close();*/

//----------------x pos ---------------------------
	myfile.open("ExXe.txt");
	for (i=0;i<ExC1;i++)
	{
		myfile<<Xe[i]<<'\n';
	}
	myfile.close();

	myfile.open("ExXh.txt");
	for (i=0;i<ExC2;i++)
	{
		myfile<<Xh[i]<<'\n';
	}
	myfile.close();

    /*
	myfile.open("MGXD0.txt");
	for (i=0;i<MGC1;i++)
	{
		myfile<<XD0[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGXDplus.txt");
	for (i=0;i<MGC2;i++)
	{
		myfile<<XDplus[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGXDminus.txt");
	for (i=0;i<MGC3;i++)
	{
		myfile<<XDminus[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGEnergy.txt");
	for (i=0;i<Max;i++)
	{
		myfile<<Def.Epos[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGXPos.txt");
	for (i=0;i<Max;i++)
	{
		myfile<<Def.xpos[i]<<'\n';
	}
	myfile.close();

	myfile.open("MGOcc.txt");
	for (i=0;i<Max;i++)
	{
		myfile<<Def.FlagOccupation[i]<<'\n';
	}
	myfile.close();*/

	delete [] Ee;
	delete [] Eh;
	delete [] ED0;
	delete [] EDminus;
	delete [] EDplus;

	delete [] Xe;
	delete [] Xh;
	delete [] XD0;
	delete [] XDminus;
	delete [] XDplus;	

	cout<<'\n'<<"Ex Electrons : "<<ExC1<<"   "<<"Ex Holes : "<<ExC2<<'\n'<<'\n';
	cout<<"MG 1 Electron : "<<MGC1<<"   "<<"MG 2 Electron : "<<MGC3<<"   "<<"MG 1 Hole : "<<MGC2<<'\n';
	cout<<'\n'<<"Unassigned States : "<<assigncheck<<'\n';
	cout<<'\n'<<"Ex States : "<<Excheck<<"   "<<"MG States : "<<MGcheck<<'\n';
	cout<<'\n'<<"Unassigned : "<<assigncheck<<'\n';

	cout<<'\n'<<"D0 : "<<MGC1<<"   "<<"D + : "<<MGC2<<"    "<<"D - : "<<MGC3<<'\n';	

	// Calculating defect density based on occupation
	MGE=MGC1/(Device_X*Device_Y*Device_Z*1E6);
	MGH=MGC2/(Device_X*Device_Y*Device_Z*1E6);
	MG2E=MGC3/(Device_X*Device_Y*Device_Z*1E6);
}

void Initialize_Occupation_2(Defect& Def)
{
    int i,j,countEx=0,Exe=0,Exh=0;
    double r,fE,HQFL;
    double *Ee,*Eh,*Xe,*Xh;

    ofstream myfile;

    Ee=new double[TotalLocalDefects];
    Eh=new double[TotalLocalDefects];
    Xe=new double[TotalLocalDefects];
    Xh=new double[TotalLocalDefects];

    for (i=0;i<xpoint;i++)
    {
        for (j=0;j<Def.SliceDefect[i];j++)
        {
            HQFL=QFL(Def.xpos[countEx]/dx);
            fE=1/(1+exp((Def.Epos[countEx]-HQFL)/Vt));
			r=(double)rand()/RAND_MAX;

            if (r<=fE)
            {
                Def.FlagOccupation[countEx]=0; // occupied by electron
                Ee[Exe]=Def.Epos[countEx];
                Xe[Exe]=Def.xpos[countEx];
                Exe++;
            }
            else
            {
                Def.FlagOccupation[countEx]=1; // unoccupied    
                Eh[Exh]=Def.Epos[countEx];
                Xh[Exh]=Def.xpos[countEx];
                Exh++;
            }
            countEx++;
        }
    }

	cout<<'\n'<<"Ex Electrons : "<<Exe<<"   "<<"Ex Holes : "<<Exh<<'\n'<<'\n';

	myfile.open("ExEe.txt");
	for (i=0;i<Exe;i++)
	{
		myfile<<Ee[i]<<'\n';
	}
	myfile.close();

	myfile.open("ExXe.txt");
	for (i=0;i<Exe;i++)
	{
		myfile<<Xe[i]<<'\n';
	}
	myfile.close();

	myfile.open("ExEh.txt");
	for (i=0;i<Exh;i++)
	{
		myfile<<Eh[i]<<'\n';
	}
	myfile.close();

	myfile.open("ExXh.txt");
	for (i=0;i<Exh;i++)
	{
		myfile<<Xh[i]<<'\n';
	}
	myfile.close();

    delete [] Ee;
    delete [] Eh;
    delete [] Xe;
    delete [] Xh;
}
