% Main: growing optimization validation
clear
delta= 8e-3;

load('DataGrow');
v0=0.0241;


opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',300000000,'TolFun',1e-20,'TolX',1e-20,'TolCon',1e-5);

x0=[ 100 0.001 100 0.02];




[xRes2,objVal2]=fminimax(@brBothGrowingOptimFct,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6  ],[1000 1000 1000 1000 ],@brBothGrowingOptimFctConst,opts,velData1,velData2,cAngle1,cAngle2,delta);
