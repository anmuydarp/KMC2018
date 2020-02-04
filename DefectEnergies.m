clc;
close all;
clear all;

[DefectEnergy]=textread('DefectEnergy.txt','%f');
[m,n]=size(DefectEnergy);

figure(1);
hist(DefectEnergy,m/1);

%[Eocc]=textread('E_OccupiedStates.txt','%f');
%[m,n]=size(Eocc);

%figure(2);
%hist(Eocc,m/1);

%[Hocc]=textread('H_OccupiedStates.txt','%f');
%[m,n]=size(Hocc);

%figure(3);
%hist(Hocc,m/1);