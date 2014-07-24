% main analyze real data OLD VERSION USE BRCOMPMOVIES
clear;

path='/home/yloosli/MatlabSave/';

load([path 'data5SecFull2']);

delta=8e-3;

velMt1=velMt1(~isnan(cAngle1T)&~isnan(cAngle2T));
velMt2=velMt2(~isnan(cAngle1T)&~isnan(cAngle2T));
cAngle1=cAngle1T(~isnan(cAngle1T)&~isnan(cAngle2T));
cAngle2=cAngle2T(~isnan(cAngle1T)&~isnan(cAngle2T));

indexGG=find(velMt1>0 & velMt2>0);

velData1=velMt1(indexGG);
velData2=velMt2(indexGG);

cAngle1=cAngle1(indexGG) ;
cAngle2=cAngle2(indexGG) ;

% velData1=velData1([2:8 10:15 18:19 21:24 ]);
% velData2=velData2([2:8 10:15 18:19 21:24 ]);
% 
% cAngle1=cAngle1([2:8 10:15 18:19 21:24 ]);
% cAngle2=cAngle2([2:8 10:15 18:19 21:24 ]);


opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',250000,'TolFun',1e-5,'TolX',1e-5,'TolCon',1e-5);
v0=0.0198;
x0=[0 0 0 0 ];

%[xResGGT,objValT,exitFlag,outPut]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst45,opts,velData1,velData2,cAngle1,cAngle2,delta/13,v0,v0);

[xRes,objVal,exitFlag,outPut]=fminimax(@brBothGrowingOptimFct,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13);
%[xResGGT,objValT,exitFlag,outPut]=fminimax(@brBothGrowingOptimFct45,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst45,opts,velData1,velData2,cAngle1,cAngle2,delta/13,v0,v0);

load([path 'data5SecFull2']);
velMt1=velMt1(~isnan(cAngle1TdT)&~isnan(cAngle2TdT));
velMt2=velMt2(~isnan(cAngle1TdT)&~isnan(cAngle2TdT));
cAngle1=cAngle1TdT(~isnan(cAngle1TdT)&~isnan(cAngle2TdT));
cAngle2=cAngle2TdT(~isnan(cAngle1TdT)&~isnan(cAngle2TdT));

indexGG=find(velMt1>0 & velMt2>0);

velData1=velMt1(indexGG);
velData2=velMt2(indexGG);

cAngle1=cAngle1(indexGG) ;
cAngle2=cAngle2(indexGG) ;
[xRes2,objVal2,exitFlag2,outPut2]=fminimax(@brBothGrowingOptimFct,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13);



figure('Name','Growth Growth case for t projection');

omega=[0:0.1:5];
plot(omega,peskin(xRes(1:2),omega,delta/13),'r');hold on;
plot(omega,peskin(xRes(3:4),omega,delta/13),'k');
plot((peskinInv(velData1,xRes(1:2),delta/13)),velData1,'dr');
plot(peskinInv(velData2,xRes(3:4),delta/13),velData2,'dk');




figure('Name','Growth Growth case for t+dt projection');

plot(omega,peskin(xRes2(1:2),omega,delta/13),'r');hold on;
plot(omega,peskin(xRes2(3:4),omega,delta/13),'k');


plot(peskinInv(velData1,xRes2(1:2),delta/13),velData1,'dr');
plot(peskinInv(velData2,xRes2(3:4),delta/13),velData2,'dk');
