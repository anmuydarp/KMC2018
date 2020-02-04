clc;
close all;
clear all;

[d2d]=textread('D2DTEST.txt','%f');

[m,n]=size(d2d);

E=1;
x=0;

for i=1:m
    x=x+1;
    r(E,x)=d2d(i);    
    Ce(E)=E-1;
    Cx(x)=x-1;
    if (mod(i,101)==0)
        E=E+1;
        x=1;
    end
end

Ce=Ce';
Cx=Cx';

figure(1);
surf(Cx,Ce,r);

[pf1]=textread('PF1TEST.txt','%f');

[m,n]=size(pf1);

E=1;
x=0;

for i=1:m
    x=x+1;
    r1(E,x)=pf1(i);    
    Ce1(E)=E-1;
    Cx1(x)=x-1;
    if (mod(i,101)==0)
        E=E+1;
        x=1;
    end
end

Ce1=Ce1';
Cx1=Cx1';

figure(2);
surf(Cx1,Ce1,r1);

[pf2]=textread('PF2TEST.txt','%f');

[m,n]=size(pf2);

E=1;
x=0;

for i=1:m
    x=x+1;
    r2(E,x)=pf2(i);    
    Ce2(E)=E-1;
    Cx2(x)=x-1;
    if (mod(i,101)==0)
        E=E+1;
        x=1;
    end
end

Ce2=Ce2';
Cx2=Cx2';

figure(3);
surf(Cx1,Ce1,r1);
hold on;
surf(Cx2,Ce2,r2);
