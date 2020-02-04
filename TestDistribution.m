clc;
close all;
clear all;

sigma=0.1;
mu=0.25;
pi=3.14;

count=0;

for x=0:0.01:.45
    count=count+1;
    f(count)=(1/(sigma*sqrt(2*pi)))*exp((-(x-mu)^2)/(2*sigma*sigma));
    xaxis(count)=x;
end

figure(1);
plot(xaxis,f);
fsdfsd
pdf=textread('Carriers.txt','%f');
[m,n]=size(pdf);

x=1:m;
x=x*1E-3;

figure(2);
plot(x,pdf);

energy=textread('CarrierEnergy.txt','%f');
[m,n]=size(energy);

x=1:m;

figure(3);
hist(energy,m/100);