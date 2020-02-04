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

void CalcAverage(CarrierH& Carr)
{
	int i,j,count=0;
	double sum=0,sumE=0,sumTau=0,sumEpen=0,sumx1=0;

	ofstream myfile;

	// Average Tau for each carrier and total average time
	for (i=0;i<CarrierMax;i++)
	{
        if (Carr.Status[i]==4) 
        {
            sum=0;
            count++;
        }
        else
        {
		    sum=0;
    		for (j=0;j<=Carr.TotPath[i];j++)
	    	{
		    	sum=sum+Carr.TauPath[i][j];
		    }
		    sumE=sumE+Carr.Energy[i];    
    		sumEpen=sumEpen+Carr.Epen[i];
	    	sumx1=sumx1+Carr.xpos[i];
		    
            Carr.TotalTime[i]=sum;
        	sumTau=sumTau+Carr.TotalTime[i];		
       }
	}
	AvgE=sumE/(CarrierMax-count);
	AvgTau=sumTau/(CarrierMax-count);
	AvgEpen=sumEpen/(CarrierMax-count);
	Avgx1=sumx1/(CarrierMax-count);

	myfile.open("TotalTime.txt");
	for(i=0;i<CarrierMax;i++)
	{
		myfile<<Carr.TotalTime[i]<<'\n';		
	}
	myfile.close();

    myfile.open("FinalEnergy.txt");
	for(i=0;i<CarrierMax;i++)
	{
		myfile<<Carr.Energy[i]<<'\n';
	}
    myfile.close();    
}

void WriteFinalPos(CarrierH& Carr)
{
	int i,j,CarrNum;

	ofstream myfile;

//------------------ Final Positions ------------------
	/*myfile.open("Finalx.txt");
	for (i=0;i<CarrierMax;i++)
	{
		myfile<<Carr.xpos[i]<<'\n';
	}
	myfile.close();*/
//------------------ Final Positions ------------------

//------------------- Path Tracking -------------------
	myfile.open("xPathTrack_Carr.txt");
    for (j=0;j<CarrierMax;j++)
    {
    	for (i=0;i<=Carr.TotPath[j];i++)
	    {
		    myfile<<Carr.xPath[j][i]<<'\n';
    	}
        myfile<<'\n';
    }
	myfile.close();

   
	myfile.open("EPathTrack_Carr.txt");
    for (j=0;j<CarrierMax;j++)
    {    	
    	for (i=0;i<=Carr.TotPath[j];i++)
	    {
		    myfile<<Carr.EPath[j][i]<<'\n';
    	}
        myfile<<'\n';
    }
	myfile.close();
	    
    myfile.open("TauPathTrack_Carr.txt");
    for (j=0;j<CarrierMax;j++)
    {
    	for (i=0;i<=Carr.TotPath[j];i++)
	    {
		    myfile<<Carr.TauPath[j][i]<<'\n';
	    }
        myfile<<'\n';
    }
	myfile.close();

	/*int temp;

	myfile.open("Ef-1.txt");
	for (i=0;i<CarrierMax;i++)
	{
		temp=Carr.TotPath[i];
		myfile<<Carr.EPath[i][temp-1]<<'\n';
	}
	myfile.close();

	myfile.open("xf-1.txt");
	for (i=0;i<CarrierMax;i++)
	{
		temp=Carr.TotPath[i];
		myfile<<Carr.xPath[i][temp-1]<<'\n';		
	}
	myfile.close();

	myfile.open("FinalTrans.txt");
	for (i=0;i<CarrierMax;i++)
	{
		myfile<<Carr.TransType[i]<<'\n';
	}
	myfile.close();

	myfile.open("FinalRate.txt");
	for (i=0;i<CarrierMax;i++)
	{
		myfile<<Carr.FinalRate[i]<<'\n';
	}
	myfile.close();*/
//------------------- Path Tracking -------------------
}
