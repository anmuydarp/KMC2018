clc;
clear all;
close all;

[Ein]=textread('Ein.txt','%f');
[m,n]=size(Ein);

x=1:m;

figure(1);
hist(Ein,m/50);

[Ef]=textread('Efinal.txt','%f');
%[Ef]=textread('Etest_Ef_Sdel.txt');
%[Ef]=textread('Etest_Ef_Sdel.txt');

figure(2);
hist(Ef,m/25);

%[tau]=textread('TransTime.txt','%f');

%[xf]=textread('Finalx.txt','%f');

%[m,n]=size(tau);

%x=1:m;

%figure(3);
%semilogy(x,tau);

%figure(4);
%plot(x,xf);

EIN_MEAN=mean(Ein)
EF_MEAN=mean(Ef)

%TAU_MEAN=mean(tau)
%xf_MEAN=mean(xf)

%[xpath]=textread('xPathTrack_Carr.txt','%f');
%[m,n]=size(xpath);

%x=1:m;

%figure(5);
%plot(xpath,x);

%[Epath]=textread('EPathTrack_Carr.txt','%f');
%[m,n]=size(Epath);

%x=1:m;

%figure(6);
%plot(Epath,x);

[TimeTrack]=textread('TotalTime.txt','%f');
[m,n]=size(TimeTrack);

x=1:m;

figure(7);
semilogy(x,TimeTrack);

count=0;
for i=1:m
    if (TimeTrack(i)>1E-12)
        count=count+1;
    end
end
count

figure(8);
plot(Ein,TimeTrack);