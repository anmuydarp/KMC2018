clc;
close all;
clear all;

[nD2D]=textread('Norm_D2D.txt','%f');
[nPF]=textread('Norm_PF.txt','%f');
[nD2E]=textread('Norm_D2E.txt','%f');

[m,n]=size(nD2D);

x=1:m;
x=x*1E-3;

figure(1);
semilogy(x,nD2D);
hold on;
semilogy(x,nPF);
hold on;
semilogy(x,nD2E);