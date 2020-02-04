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

#define FUNC(x) ((*func)(x))

double TunCoeff (double x)
{
	double Term,EvTail;

    //Find Ev function based on x

	//EvTail=VBOffset-(x*VBOffset/(ithickness));
    EvTail=BandEdge(x/dx,0,0);

	//EvTail=3.15-pow(x/(10*1E-9),2.0); - Testing for figures from Liu paper

	Term=sqrt((2*m_aSi*m0/pow(hbar,2.0))*(EvTail*q-ETrans*q));

	if (EvTail<ETrans)
		Term=0;

	return Term;
}

double sigVdelta (double x)
{
	double V,del,term1,term2,I1,I2,total,constant;

  	constant=Ud/(sqrt(pow(Gradi,3.0))*sqrt(pow(Gradj,3.0))*4*pi/3);

  	term1=3*Ud*pow(4*pi/3,1/3.);
  	term2=(1/pow(Grdj,2.0))*sinc(x*GRj/2)-(1/pow(Grdi,2.0))*sinc(x*GRi/2);
  	del=term1*term2;
  	del=del/(sqrt(2)*homega*q);
  
  	//cout<<"term1 : "<<term1/(homega*q)<<"  "<<"term2 : "<<term2<<'\n';
  
  	term1=-(4*Galpha/sin(Galpha))*(Gradi*sinc(x*Gradi*sin(Galpha))-Gx*sinc(x*Gx*sin(Galpha)));
  	term2=(4*Galpha/sin(Galpha))*(Gradi*cos(x*Gradi*sin(Galpha))-Gx*cos(x*Gx*sin(Galpha)));  	
  	I1=(term1+term2)*constant;
  
  	term1=-(4*Galpha/sin(Galpha))*(Gradj*sinc(x*Gradj*sin(Galpha))-(Gd-Gx)*sinc(x*(Gd-Gx)*sin(Galpha)));
  	term2=(4*Galpha/sin(Galpha))*(Gradj*cos(x*Gradj*sin(Galpha))-(Gd-Gx)*cos(x*(Gd-Gx)*sin(Galpha)));  	
  	I2=(term1+term2)*constant;
  
  	V=I1+I2;

  	//cout<<"V : "<<V<<'\n';

  	//del=2.449*q; // For constant del case

  	total=V*del;

  	//cout<<"Sig Vdelta : "<<"V : "<<V<<"   "<<"del : "<<del<<"   "<<"total : "<<total<<'\n';

  	return total;  	
}

double sigV (double x)
{
	double V,term1,term2,I1,I2,constant;

  	constant=Ud/(sqrt(pow(Gradi,3.0))*sqrt(pow(Gradj,3.0))*4*pi/3);
  
  	term1=-(4*Galpha/sin(Galpha))*(Gradi*sinc(x*Gradi*sin(Galpha))-Gx*sinc(x*Gx*sin(Galpha)));
  	term2=(4*Galpha/sin(Galpha))*(Gradi*cos(x*Gradi*sin(Galpha))-Gx*cos(x*Gx*sin(Galpha)));  	
  	I1=(term1+term2)*constant;
  
  	term1=-(4*Galpha/sin(Galpha))*(Gradj*sinc(x*Gradj*sin(Galpha))-(Gd-Gx)*sinc(x*(Gd-Gx)*sin(Galpha)));
  	term2=(4*Galpha/sin(Galpha))*(Gradj*cos(x*Gradj*sin(Galpha))-(Gd-Gx)*cos(x*(Gd-Gx)*sin(Galpha)));  	
  	I2=(term1+term2)*constant;
  
  	V=I1+I2;

  	//cout<<"sigV V: "<<V<<'\n';

  	return V;
}

double sigDelta (double x)
{
	double del,term1,term2;

  	term1=3*Ud*pow(4*pi/3,1/3.);
	term2=(1/pow(Grdj,2.0))*sinc(x*GRj/2)-(1/pow(Grdi,2.0))*sinc(x*GRi/2);
  	del=term1*term2;
  	del=del/(sqrt(2)*homega*q);

	return del;
}

double InelasticTransition(double E,double Ed,double x,int p,int t)
{
	// In this function Ed is defect depth
	double R=0,F=aSi_EField,Ec=0,Ef=0; // F = V/m, hw = eV, S = Huang Rhys Factor, Ef = Fermi Level is pinned at 0 (assumption)

	double DOS,fE,TC,c0,hbarE0,Ip; // DOS = per m3 per J (can convert into per cm3), Ip = Modified Bessel Function, hbarE0 = electro - optical energy
	double rd,fBE,z,Lp;
	
	//R=c0*N(Ep)*f(Ep)*T(Ep,x)*Lp(z)       

    rd=hbar/sqrt(2*m_aSi*m0*m_aSi*Ed*q);
    hbarE0=pow(pow(q*hbar*F,2.0)/(2*m_aSi*m0),1/3.);
	c0=pow(hbarE0,3.0)*pow(4*pi,2.0)*pow(rd,3.0)/(hbar*Eg_aSi*q);

    ETrans=E;

	TC=exp(-2*qsimp(&TunCoeff,initialx,x));	  

    fBE=1/(exp(homega/Vt)-1);
	z=2*S*sqrt(fBE*(fBE+1));

	fE=1/(1+exp((E-0)/Vt)); // for t = 0 - f(E) is either 1 or 0 depending on initial occupation    

    if (abs(p)==0)
	    Ip=bessi0(z);
   	else if (abs(p)==1)
    	Ip=bessi1(z);
   	else if (abs(p)>1)
    	Ip=bessi(abs(p),z);

   	Lp=pow((fBE+1)/fBE,p/2.)*exp(-S*(2*fBE+1))*Ip;   

	switch (t)
	{
		case 0: // INJECTION Electrode --> Defect

        Ec=0;

        if (E<Ec)
            DOS=0;
        else if (E>=Ec)
            DOS=((8*pi*sqrt(2))/pow(hbar*2*pi,3.0))*(pow(m_Si*m0,1.5))*sqrt(E-Ec)*sqrt(q); // This is the DOS of the Carriers

		if (p<0) // p<0 - Absorption
        {
			R=c0*DOS*TC*fE*Lp*exp(p*homega/Vt);	
        }
		else if (p>0) // p>0 - Emission
        {
			R=c0*DOS*TC*Lp*fE; // Should be 1-fE but 1-fE is not occupied by electrons, and no transitions to states with holes are allowed	
        }
		break;

		case 1: // EMISSION Defect --> Electrode
        
        Ec=BandEdge(Device_X/dx,0,0);
        
        if (E<Ec)
            DOS=0;
        else if (E>=Ec)
            DOS=((8*pi*sqrt(2))/pow(hbar*2*pi,3.0))*(pow(m_Si*m0,1.5))*sqrt(E-Ec)*sqrt(q); // This is the DOS of the Carriers

		if (p<0) // p<0 - Emission
			R=c0*DOS*TC*Lp*(1-fE);
		else if (p>0) // p>0 - Absorption
			R=c0*DOS*TC*Lp*(1-fE)*exp(-p*homega/Vt); 

		break;
	}

	if (isnan(R)==1)
		raise(SIGSEGV);  

	return R;
}

double ElasticTransition(double E,double Ed,double x,int t)
{
	double R,TC,fE,tempxi,tempxf,VBmin=0;	

	// For electrode to defect : R = C(E)f(E)T(E)

	ETrans=E;	

	if (initialx<x)
	{
		tempxi=initialx;
		tempxf=x;
	}
	else
	{
		tempxi=x;
		tempxf=initialx;
	}

	TC=exp(-2*qsimp(&TunCoeff,tempxi,tempxf));

	fE=1/(1+exp(E/Vt));

	if (t==0) // Electrode to Defect
    {
        // The bottom of the band is 0
		R=pow(m_Si/m_aSi,5/2.)*(8*pow(E*q,1.5)/(3*hbar*sqrt(Ed*q)))*TC*fE;
    }
	else if (t==1) // Defect to Electrode
    {
        // The bottom of the band is VBmin
        VBmin=round(BandEdge(Device_X/dx,0,0)*1e3)/1e3;
		R=pow(m_Si/m_aSi,5/2.)*(8*pow((E-VBmin)*q,1.5)/(3*hbar*sqrt(Ed*q)))*TC*(1-fE);
    }

	if (isnan(R)==1)
    {
        cout<<"t : "<<t<<"    "<<"E : "<<E<<'\n';
        cout<<"VBmin : "<<round(BandEdge(Device_X/dx,0,0)*1e3)/1e3<<'\n';
        cout<<"R : "<<R<<'\n';
		raise(SIGSEGV);
    }

	return R;
}

double PF(int i) // POOLE FRENKEL EMISSION
{
	double R,DefectDepth,E,F=aSi_EField,deltaPF,EvTail; // V/m

	/*****************************************************************************************
	NOTE: 1 - Equation 1 is for PF emission parallel to the oxide or in our case the aSi-cSi
	          Heterointerface, thus this is 1 D PF emission.
	      2 - Equation 2 is for PF emission at any angle w.r.t to the interface. YET TO BE 
	          IMPLEMENTED
	*****************************************************************************************/
	E=i*EnergyStep;	
    EvTail=BandEdge(initialx/dx,0,0);
	DefectDepth=EvTail-E;

	deltaPF=sqrt(pow(q,3.0)*F/(pi*eps*eopt_aSi));
	
	R=LatticeVibration*exp(-(1/(kb*T))*(DefectDepth*q-deltaPF)); // Equation 1

	if (DefectDepth<0)
	{	
		raise(SIGSEGV);
	}

	return R;
}

double D2D(double E,double xi,int j,Defect& Def)
{
	double R,xf,E1,E2,DelE,TC; // E1 is the position of the initial defect, E2 is the position of the final defect
	double tempxi,tempxf;

	E1=E; // Initial Energy
	E2=Def.Epos[j]; // Final Energy
	xf=Def.xpos[j];	

	ETrans=E1; 

	/***************************************************************************
	Elastic Tunneling from defect 'i' to the virtual position 'xj' of defect 'j'
	and subsequent relaxation to energy 'Ej'
	***************************************************************************/  

	if (xi>xf)
	{
		tempxf=xi;
		tempxi=xf;    	
	}
	else 
	{
		tempxi=xi;
		tempxf=xf;    
	}

	TC=exp(-2*qsimp(&TunCoeff,tempxi,tempxf));

	if (E2>E1) // Absorption
	{	
        DelE=E2-E1; // if DelE > 0, it is an emission process, as E2 < E1, if DelE = 0, then it is an elastic process
		R=LatticeVibration*TC*exp(-DelE/Vt);
	}
	else if (E2<=E1)
	{			
		R=LatticeVibration*TC;
	}
	return R;
}

double PhD2D(int i,int j,int p,double Ep,Defect& Def)
{
	double R,R0,Rsum,TC,c,term1,term2,term3,tempxf,tempxi;
	double V,delta,Edi,Edj,Vdelta;
	double xi,xf,fBE,newx,EvTail,newS;	
	double a=5.431E-10;
			
	// i - initial state , j - final state
	// Deltaq = <psij|U|psij> - <psii|U|psii>

	/*cout<<'\n';
	cout<<"Initial Energy : "<<i*EnergyStep<<'\n';
	cout<<"Defect no : "<<j<<'\n';
	cout<<"Final Energy (coupled) : "<<Ep<<'\n';*/	

	newx=int((initialx)*10)/10;
	//EvTail=VBOffset-(newx*dx*VBOffset/ithickness);
    EvTail=BandEdge(initialx/dx,0,0);
	Edi=EvTail-i*EnergyStep;
	newx=int((Def.xpos[j])*10)/10;
	//newx=int(finalx*10)/10;
	EvTail=VBOffset-(newx*dx*VBOffset/ithickness);
	Edj=EvTail-Def.Epos[j];
	//Edj=EvTail-0.13;

	if (Edj<0 || Edi<0)
		raise(SIGSEGV);

	//cout<<"Initial Energy : "<<i*EnergyStep<<"  "<<"Final Energy : "<<Ep<<'\n';

	//cout<<"Edi : "<<Edi<<"   "<<"Edj : "<<Edj<<'\n';
	
	// Localization radius
	Grdi=hbar/sqrt(2*m_aSi*m0*Edi*q);
	Grdj=hbar/sqrt(2*m_aSi*m0*Edj*q);

	//cout<<"rdi : "<<Grdi<<"  "<<"rdj : "<<Grdj<<'\n';
	
	// Radius of Sphere - to be used for delta q terms
	GRi=pow(4*pi*pow(Grdi,3.0)/3,1/3.);
	GRj=pow(4*pi*pow(Grdj,3.0)/3,1/3.);

	//cout<<"Ri : "<<GRi<<"   "<<"Rj : "<<GRj<<'\n';
	
	// Intersection of Sphere - to be used for Vq terms
	Gradi=Grdi;
	Gradj=Grdj;
	//Gd=abs(Def.xpos[j]-initialx);
	Gd=abs(finalx-initialx);
	Gx=(abs(pow(Gradi,2.0)-pow(Gradj,2.0))+pow(Gd,2.0))/(2*Gd);
	Galpha=acos(Gx/Gradi);
	
	//cout<<"d : "<<Gd<<"  "<<"x : "<<Gx<<"   "<<"alpha : "<<Galpha<<'\n';	

	if (Gx>Gradi)
	{
		//cout<<"No coupling achieved"<<'\n';
		return 0;
	}
	
	Vdelta=qsimp(&sigVdelta,0,pi/a);

	delta=qsimp(&sigDelta,0,pi/a);

	V=qsimp(&sigV,0,pi/a);

	newS=pow(abs(delta),2.0);

	//cout<<"Vdelta : "<<Vdelta<<'\n';
	//cout<<"x : "<<finalx<<"   "<<"delta : "<<delta<<"   "<<"S : "<<newS<<'\n';
	//cout<<"V : "<<V<<'\n';

	term1=(pi/(hbar*homega*q));
	//term2=pow(S,abs(p)-1)*exp(-S)/factorial(abs(p)-1);
	term2=pow(newS,abs(p)-1)*exp(-newS)/factorial(abs(p)-1);
	//term3=pow(abs(V),2.0)*0.26+pow(abs(Vdelta),2.0)*(abs(p)-1)*0.18/S;
	term3=pow(abs(V),2.0)*0.26+pow(abs(Vdelta),2.0)*(abs(p)-1)*0.18/newS;

	//cout<<"R0 term1 : "<<term1<<"  "<<"R0 term2 : "<<term2<<"   "<<"R0 term3 : "<<term3<<'\n';
	
	R0=term1*term2*term3;

	//cout<<"R0 : "<<R0<<'\n';
	
	fBE=1/(exp(homega/Vt)-1);

	ETrans=Ep;

	if (initialx>Def.xpos[j])
	{
		tempxi=Def.xpos[j];
		tempxf=initialx;
	}
	else
	{
		tempxi=initialx;
		tempxf=Def.xpos[j];
	}

	TC=exp(-2*qsimp(&TunCoeff,tempxi,tempxf));
	//TC=e(-2*qsimp(&TunCoeff,initialx,finalx));
	
	//cout<<"TC : "<<TC<<'\n';	

	if (p>0)
	{
		//Absorption
		//c=R0*pow(fBE,p)*exp(-2*fBE*S);
		c=R0*pow(fBE,p)*exp(-2*fBE*newS);		
	}
	else if (p<0)
	{
		//Emission
		//c=R0*pow(fBE+1,-p)*exp(-2*fBE*S);
		c=R0*pow(fBE+1,-p)*exp(-2*fBE*newS);
	}
	
	//cout<<"c : "<<c<<"   "<<"Phonon Mode : "<<p<<'\n';
	//cout<<"c term 2 : "<<term2<<"  "<<"c term 3 : "<<term3<<'\n';

	R=TC*c;

	if (isnan(R)!=0)
		raise(SIGSEGV);
		
	//cout<<"PhD2D R : "<<R<<'\n';
	
	return R;
}

double TCTest(double E,double xi,double xf)
{
	double TC,tempxi,tempxf;

	ETrans=E;

	if (xi>xf)
	{
		tempxi=xf;
		tempxf=xi;		
	}
	else
	{
		tempxi=xi;
		tempxf=xf;
	}

	TC=exp(-2*qsimp(&TunCoeff,tempxi,tempxf));

	return TC;
}

int QMTCReflection(double E,double x)
//double QMTCReflection(double E,double x)
{	
	int Q;
	double T,R,r,L1,L2;

	int i;

	ETrans=E;

	T=exp(-2*qsimp(&TunCoeff,0,x));
	
	R=1-T;

	r=(double)rand()/RAND_MAX;

	if (r<T)
		Q=0; // Transmission
	else
		Q=1; // Reflection

	//cout<<"T : "<<T<<"   "<<"R : "<<R<<"    "<<"r : "<<r<<"   "<<"Q : "<<Q<<'\n';

	return Q;
//	return T;
}

double RecombinationRate(int r,double x,double Et)
{
	// This function calculates hole emission coefficient or electron capture co-efficient
	double rate,Ev,Interface_HoleConc=4E18,Interface_ElecConc=4E14,Bulk_ElecConc=1E12,Bulk_HoleConc=1E14,Ec=1.27;

	switch (r)
	{
		//---------------------- RATES FOR HOLES ------------------------------

		//******************** BULK ********************
		case 0: // Hole Emission for BULK RECOMBINATION
		Ev=VBOffset-(x*VBOffset/(ithickness));
		Ev=-Ev;
		Et=-Et;
		rate=Vth*sigma_plus*Nv_aSi*exp((Ev-Et)/Vt);
		//cout<<'\n'<<"Hole Emission Rate : "<<rate<<'\n';
		break;

		case 1: // Hole Capture for BULK RECOMBINATION
		rate=(MGE+MG2E)*Vth*sigma_plus*Bulk_HoleConc;
		break;
		//******************** BULK ********************

		//****************** INTERFACE *****************
		case 2: // Hole Emission for INTERFACE RECOMBINATION
		//Ev=VBOffset-(x*VBOffset/(ithickness));
		Ev=0;		
		Et=-Et;
		rate=Vth*sigma_plus*Nv_aSi*exp((Ev-Et)/Vt);
		break;

		case 3: // Hole Capture from the valence band for INTERFACE RECOMBINATION		
		rate=IntElec*Vth*sigma_plus;
		break;
		//****************** INTERFACE *****************

		//---------------------- RATES FOR HOLES ------------------------------

		//-------------------- RATES FOR ELECTRONS ----------------------------

		//******************** BULK ********************
		case 4: // Electron Emission for BULK RECOMBINATION
		Ev=VBOffset-(x*VBOffset/(ithickness));
		Ev=-Ev;
		Ec=Ev+Eg_aSi;
		rate=Vth*sigma_plus*Nc_aSi*exp((Et-Ec)/Vt);
		break;

		case 5: // Electron Capture for BULK RECOMBINATION		
		rate=MGH*Vth*sigma_plus;
		break;
		//******************** BULK ********************

		//****************** INTERFACE *****************
		case 6: // Electron Capture from the conduction band for INTERFACE RECOMBINATION
		//rate=IntHole*Vth*sigma_plus*Interface_ElecConc/Nc_aSi;
		rate=IntHole*Vth*sigma_plus;		
		break;

		case 7: // Electron Emission for Interface RECOMBINATION
		Ec=1.12;
		rate=Vth*sigma_plus*Nc_aSi*exp((Et-Ec)/Vt);
		break;
		//****************** INTERFACE *****************

		//-------------------- RATES FOR ELECTRONS ----------------------------
	}
	return rate;
}
