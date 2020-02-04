clc;
clear all;
close all;


for i=1:100
    test(i)=rand();
    x(i)=i;
end

sum=0;
for i=2:100
    sum=sum+test(i-1);
    temp(i)=sum;    
end
temp(1)=test(1);

temp=temp/temp(100);

figure(1);
plot(x,temp);

rnd=rand()

for i=2:100
    if (temp(i-1)<rnd && temp(i)>=rnd)
        prev=temp(i-1)
        next=temp(i)
        search=i
        break;
    end
end

%Binary Search

target=rnd
first=1;
last=100;
middle=floor((first+last)/2);

flag=0;

while (flag==0)   
    if (target>temp(middle))
        first=middle;        
    elseif (target<temp(middle))
        last=middle;        
    end
    
    t=last-first;
    if (t==1)
        flag=1;
        first
        last
    end    
    middle=floor((first+last)/2);
end