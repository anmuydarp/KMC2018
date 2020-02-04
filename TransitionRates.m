clc; 
clear all;
close all;

% [rc50]=textread('NormTableRC50.txt','%f');
% [rc150]=textread('NormTableRC150.txt','%f');
% [rc250]=textread('NormTableRC250.txt','%f');
% [rc350]=textread('NormTableRC350.txt','%f');
% 
% [m,n]=size(rc50);
% 
% x=1:m;
% 
% figure(1);
% plot(x,rc50);
% hold on;
% plot(x,rc150);
% hold on;
% plot(x,rc250);
% hold on;
% plot(x,rc350);

[d1]=textread('D2DRate_150_1_3.txt','%f');
[d2]=textread('D2DRate_150_1_5.txt','%f');
[d3]=textread('D2DRate_150_1_7.txt','%f');

[m,n]=size(d1);

x=1:m;
x=x*1E-3;

figure(1);
semilogy(x,d1);
hold on;
semilogy(x,d2);
hold on;
semilogy(x,d3);

[d1]=textread('D2DRate_080_1_3.txt','%f');
[d2]=textread('D2DRate_080_1_5.txt','%f');
[d3]=textread('D2DRate_080_1_7.txt','%f');

[m,n]=size(d1);

x=1:m;
x=x*1E-3;

figure(2);
semilogy(x,d1);
hold on;
semilogy(x,d2);
hold on;
semilogy(x,d3);
