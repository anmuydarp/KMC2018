clc;
clear all;
close all;

%[E]=textread('RecombStats.txt','%f');

%[m,n]=size(E);

%Eavg=mean(E)

VBOffset=0.45;
ithick=10;
Vth=2E7;
sigma=2.5E-15;
Nv_aSi=2.5E20;
Et=-0.45;
ElecConc=1E-5;
Vt=0.0258;

count=0;
count1=0;

for x=0:0.1:10
    count=count+1;
    Ev=VBOffset-(x*VBOffset/ithick);
    count1=0;
    for Et=-0.3:0.01:Ev
        count1=count1+1;
        HEm(count1,count)=Vth*sigma*Nv_aSi*exp(-(Ev-Et)/Vt);
        x1(count1)=Et;        
    end
        ECap(count)=Vth*sigma*ElecConc;        
        x2(count)=x;  
end

figure(1);
surf(x2,x1,HEm);

figure(2);
semilogy(x2,HEm(1,:));