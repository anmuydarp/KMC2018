clc;
clear all;
close all;

[E]=textread('DefectEnergy.txt','%f');
[me,n]=size(E);
[Fermi]=textread('DefectFermi.txt','%f');
[mF,n]=size(Fermi);

%[B,I]=sort(E);

%figure(1);
%plot(B,Fermi(I));

x=1:me;
%figure(2);
%hist(E,me/350);

[Ed]=textread('EdvxE_Ed.txt','%f');
[m,n]=size(Ed);

x=1:m;

figure(3);
hist(Ed,m/1);

[xE]=textread('EdvxE_x.txt','%f');
[m,n]=size(xE);

x=1:m;

figure(4);
hist(xE,m/1);

[occ]=textread('DefOcc.txt','%f');

[Bx,I]=sort(xE);

[m,n]=size(Bx);

count=0;
count1=0;
for i=1:m
    if (occ(i)==0) %Hole
        count=count+1;
        By1(count)=Ed(I(i));
        Bx1(count)=Bx(i);
    elseif (occ(i)==1) %Electron
        count1=count1+1;
        By2(count1)=Ed(I(i));
        Bx2(count1)=Bx(i);        
    end    
end

figure(5);
plot(Bx1,By1);
hold on;
plot(Bx2,By2);