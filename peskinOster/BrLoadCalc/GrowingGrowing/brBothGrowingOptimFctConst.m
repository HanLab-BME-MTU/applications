function [c,ceq]=brBothGrowingOptimFctConst(x,velData1,velData2,cAngle1,cAngle2,delta);

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
v0=0.024;
c(1)=x(2)-x(1)+0.0001;
c(2)=x(4)-x(3)+0.0001;
c(3)=v1(end)-5e-2;
c(4)=v2(end)-5e-2 ;
c(5)=0.75*max(velData1)-max(v1);
c(6)=max(v1)-1.25*max(velData1);
c(7)=0.75*max(velData2)-max(v2);
c(8)=max(v2)-1.25*max(velData2);
% c(9)=x(1)-3*max(velData1)/(delta*(1-exp(-5)));
% c(10)= x(3)-3*max(velData2)/(delta*(1-exp(-5)));

c=[ -v1 -v2 -omega1 -omega2  c    ]';

