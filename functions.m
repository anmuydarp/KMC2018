clc;
clear all;
close all;

%Testing phonon sensitivity

Vt=0.0258;
p=1;
S=0.5;
count=0;
count1=0;

for p=1:5
    count1=count1+1;
    for h=0.005:0.005:0.4
        count=count+1;
        f=1/(exp(h/Vt)-1);
        z=2*S*sqrt(f*(f+1));
        Lp(count)=((f+1)/f)^(p/2)*exp(-S*(f+1))*besseli(p,z);        
        x(count)=h;
    end
    figure(1);
    semilogy(x,Lp);
    hold on;
    count=0;
end

%Coupled State Check

[Ep]=textread('CoupledStateAnalysis.txt','%f');
[Ei]=textread('InitialDefectState.txt','%f');
[R]=textread('CSRate.txt','%f');

[m,n]=size(Ep);
[m1,n1]=size(Ei);

x=1:m;
x1=1:m1;

figure(2);
plot(x,Ep);
hold on;
plot(x1,Ei);

count=0;
for i=1:m
    if (Ep(i)>=0)
        count=count+1;
        RealEp(count)=Ep(i);
        RealEi(count)=Ei(i);
        RealR(count)=R(i);
        RealX(count)=count;
    end    
end

figure(3);
plot(RealX,RealEp);
hold on;
plot(RealX,RealEi);

figure(4);
semilogy(RealX,RealR);

MEAN_EP=mean(RealEp)

MEAN_Ei=mean(RealEi)

DelMean=MEAN_Ei-MEAN_EP
