/********************************************************************/
/*  General Physics Parameters for Kinetic Monte Carlo  	        */
/********************************************************************/

#ifndef GPPARAMETERS_H
#define GPPARAMETERS_H

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double q=1.60219e-19; // Coulomb
const double m0=9.1E-31; // Kg
const double hbar=1.054E-34; //Jsec
const double eps=8.85419e-12; // F/m
const double kb=1.38066e-23; // J/K
const double pi=3.1416;
const double T=300; // K
const double Vt=kb*T/q; // eV

//Parameters for Mid Gap states
const double Tstar=500; // K - Equilibration center
const double KT=kb*Tstar/q; // Vt at equilibration temperature
const double EvT0=0.036; // eV - characteristic decay energy
const double Ev0=sqrt(pow(EvT0,2.0)+pow(KT,2.0)); // Ev0 at equilibration temperature

#endif
