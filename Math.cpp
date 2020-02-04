#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <stdlib.h>
#include <csignal>
#include "Global.h"
#include "nrutil.h"

using namespace std;

#define EPS 1.0e-6
#define JMAX 20

double trapzd(double (*func)(double),double a,double b,int n)
{
double x,tnm,sum,del;
static double s;
int it,j;

if (n == 1) {
return (s=0.5*(b-a)*(func(a)+func(b)));
} else {
for (it=1,j=1;j<n-1;j++) it <<= 1;
tnm=it;
del=(b-a)/tnm;

x=a+0.5*del;

for (sum=0.0,j=1;j<=it;j++,x+=del) sum += func(x);
s=0.5*(s+(b-a)*sum/tnm);

return s;
}
}

double qsimp(double (*func)(double),double a,double b)
{
int j;
double s,st,ost,os;
ost = os = -1.0e30;
for (j=1;j<=JMAX;j++) {
st=trapzd(func,a,b,j);
s=(4.0*st-ost)/3.0;

if (j > 5)

if (fabs(s-os) < EPS*fabs(os) ||
(s == 0.0 && os == 0.0)) return s;
os=s;
ost=st;
}
return 0.0;
}

float bessi0(float x)
{
float ax,ans;
double y;

if ((ax=fabs(x)) < 3.75) {
y=x/3.75;
y*=y;
ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
} else {
y=3.75/ax;
ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
+y*0.392377e-2))))))));
}
return ans;
}

float bessi1(float x)
{
float ax,ans;
double y;

if ((ax=fabs(x)) < 3.75) {

y=x/3.75;
y*=y;
ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
} else {
y=3.75/ax;
ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
-y*0.420059e-2));
ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
ans *= (exp(ax)/sqrt(ax));
}
return x < 0.0 ? -ans : ans;
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

float bessi(int n, float x)
{
float bessi0(float x);
void nrerror(char error_text[]);
int j;
float bi,bim,bip,tox,ans;

//if (n < 2) nrerror("Index n less than 2 in bessi");
if (x == 0.0)
return 0.0;
else {
tox=2.0/fabs(x);
bip=ans=0.0;
bi=1.0;
for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
bim=bip+j*tox*bi;
bip=bi;
bi=bim;
if (fabs(bi) > BIGNO) {
ans *= BIGNI;
bi *= BIGNI;
bip *= BIGNI;
}
if (j == n) ans=bip;
}
ans *= bessi0(x)/bi;

return x < 0.0 && (n & 1) ? -ans : ans;
}
}

int factorial(int x)
{
	int f,i,temp=1;

	i=x;

	if (x>=1)
	{
		while (i>=1)
		{
			temp=temp*i;	
			i--;
		}		
	}
	else if (x==0)
		temp=1;
	
	f=temp;
	return f;
}

double sinc (const double x)
{
	if (x==0)
		return 1;
	
	return sin(x)/x;
}

double parallel (double a, double b)
{
	double c;

	c=(a*b)/(a+b);

	return c;
}