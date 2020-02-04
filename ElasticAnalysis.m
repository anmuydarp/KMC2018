clc;
clear all;
close all;

[e]=textread('ElasticTT_E.txt');
[x]=textread('ElasticTT_x.txt');
[r]=textread('ElasticTT_R.txt');

[m,n]=size(e);

k=0;
j=0;
for i=1:m/3
    i
    k=k+1;
    j=j+1;
    E(j)=e(i);
    X(k)=x(i);
    R(j,k)=r(i);
end

figure(1);
surf(E,X,R);
