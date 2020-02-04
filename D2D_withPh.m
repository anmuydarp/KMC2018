clc;
clear all;
close all;

[d2demp1]=textread('15Dec_D2DEmP1.txt','%f');
[d2dabp1]=textread('15Dec_D2DAbP1.txt','%f');

[phd2demp1]=textread('15Dec_PhD2DEm_delP1.txt','%f');
[phd2dabp1]=textread('15Dec_PhD2DAb_delP1.txt','%f');

[phd2dems1]=textread('15Dec_PhD2DEm_S1.txt','%f');
[phd2dabs1]=textread('15Dec_PhD2DAb_S1.txt','%f');

[d2demp2]=textread('15Dec_D2DEmP2.txt','%f');
[d2dabp2]=textread('15Dec_D2DAbP2.txt','%f');

[phd2demp2]=textread('15Dec_PhD2DEm_delP2.txt','%f');
[phd2dabp2]=textread('15Dec_PhD2DAb_delP2.txt','%f');

[phd2dems2]=textread('15Dec_PhD2DEm_S2.txt','%f');
[phd2dabs2]=textread('15Dec_PhD2DAb_S2.txt','%f');


[m,n]=size(d2demp1);

x=1:m;
x=(x-1);
x=x';

figure(1);
semilogy(x,d2demp1);
hold on;
semilogy(x,d2dabp1);
hold on;
semilogy(x,phd2demp1);
hold on;
semilogy(x,phd2dabp1);
hold on;
semilogy(x,phd2dems1);
hold on;
semilogy(x,phd2dabs1);

figure(2);
semilogy(x,d2demp2);
hold on;
semilogy(x,d2dabp2);
hold on;
semilogy(x,phd2demp2);
hold on;
semilogy(x,phd2dabp2);
hold on;
semilogy(x,phd2dems2);
hold on;
semilogy(x,phd2dabs2);