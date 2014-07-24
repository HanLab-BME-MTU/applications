function [c,ceq]=brBothGrowingOptimFctConst(x,velData1,velData2,cAngle1,cAngle2,delta);
% brBothGrowingOptimFct`Const generate the obj. fct for the growth-growth case
% INPUT
%     x         :  paramters of pol
%     velData1  :  velocity of the MT1 w>0
%     velData2  :  velocity of MT2 >0
%     delta     :  1/13 of subunit length
%     without
% OUPUT 
%     c         :  inequality constraint

ceq=[];

alpha1=x(1);
beta1=x(2);
alpha2=x(3);
beta2=x(4);


omega=[0:0.1:5];
v1=delta*(alpha1*exp(-omega)-beta1);
v2=delta*(alpha2*exp(-omega)-beta2);

omega1=peskinInv(velData1,x(1:2),delta);
omega2=peskinInv(velData2,x(3:4),delta);

if size(v1)~=size(omega1)
    omega1=omega1';
    omega2=omega2';
end
v01=max(velData1);
v02=max(velData2);
c(1)=x(2)-x(1)+0.0001;
c(2)=x(4)-x(3)+0.0001;
c(3)=v1(end)-5e-2;
c(4)=v2(end)-5e-2 ;
c(5)=0.75*v01-max(v1);
c(6)=max(v1)-1.25*v01;
c(7)=0.75*v02-max(v2);
c(8)=max(v2)-1.25*v02;
% c(9)=x(1)-3*max(velData1)/(delta*(1-exp(-5)));
% c(10)= x(3)-3*max(velData2)/(delta*(1-exp(-5)));

c=[ -v1 -v2 -omega1 -omega2  c    ]';

