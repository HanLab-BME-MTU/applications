function [c,ceq]=brBothGrowingOptimFctConst45(x,velData1,velData2,cAngle1,cAngle2,delta,v01,v02);
% brBothGrowingOptimFctConst45 generate the obj. fct for the growth-growth case
% with the unload vel.
% INPUT
%     x         :  paramters of pol
%     velData1  :  velocity of the MT1 w>0
%     velData2  :  velocity of MT2 >0
%     delta     :  1/13 of subunit length
%     v01       : unload vel of the first MT
%     v02       : unload vel of the second MT
%     without
% OUPUT 
%     c         :  constraint inquality


ceq=[];

alpha1=x(1);
beta1=x(2);
alpha2=x(3);
beta2=x(4);


omega=[0:0.1:5];
v1=delta*(alpha1*exp(-omega)-beta1);
v2=delta*(alpha2*exp(-omega)-beta2);
% indexTooLarge1=find(velData1>v01);
% indexTooLarge2=find(velData2>v02);
omega1=peskinInv(velData1,x(1:2),delta);
omega2=peskinInv(velData2,x(3:4),delta);
% omega1(indexTooLarge1)=1;
% omega2(indexTooLarge2)=1;

if size(v1)~=size(omega1)
    omega1=omega1';
    omega2=omega2';
end
v0=0.024;
c(1)=x(2)-x(1)+0.0001;
c(2)=x(4)-x(3)+0.0001;
c(3)=v1(end)-5e-2;
c(4)=v2(end)-5e-2 ;
c(5)=0.6*v01-max(v1);
c(6)=max(v1)-1.4*v01;
c(7)=0.6*v02-max(v2);
c(8)=max(v2)-1.4*v02;
% c(9)=x(1)-3*max(velData1)/(delta*(1-exp(-5)));
% c(10)= x(3)-3*max(velData2)/(delta*(1-exp(-5)));

c=[ -v1 -v2 -omega1 -omega2  c    ]';

