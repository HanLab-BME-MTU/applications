%  main function

clear;
delta=8e-3;
v0=0.03;

load('velDataGeneral.mat');



opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',150,'TolFun',1e-5,'TolX',1e-5,'TolCon',1e-5);


alpha1=9;
betaG1=0.02;
beta1=1.1 ; 
gamma1=29;

alpha2=12;
betaG2=0.08;
beta2= 0.1;
gamma2=36;
% 
% %xres=fmincon(@brDirectOptim,x0,[],[],[],[],[1e-4 1e-4  1e-4 1e-1 1e-4 1e-4 1e-4  1e-4 1e-4 1e-4 ],[10000 10000  10000 10000 1  10000 10000  10000 10000 1 ],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,'single');
% 
% fct=brDirectOptim(xGS,velData1(indexGS),velData2(indexGS),cAngle1(indexGS),cAngle2(indexGS),delta,vUnload,'GS','multi','fixed');
% fctg=brDirectOptim(xSG,velData1(indexSG),velData2(indexSG),cAngle1(indexSG),cAngle2(indexSG),delta,vUnload,'SG','multi','fixed');
% 
fctf=brDirectOptim(x2,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'All','multi','fixed');


x0=[10 0.1 0.1 10  ];
x0=[10 0.1 0.1 10 [10 0.1 0.1 10 ]*1.09];



xObj=[alpha1 betaG1 beta1 gamma1 alpha2 betaG2 beta2 gamma2];

%[xres,objVal]=fminimax(@brDirectOptim,[x0 (omega2.*cAngle2)'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4 1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1))],[10000 10000  10000 10000  10000 10000  10000 10000   5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,v0,'All','multi','fixed');
%
% fct=brGeneralOptim(x2(1:8),velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
% [c,ceq]=brGeneralOptimConst(x2(1:8),velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
 %[xResOld,objValOLD]=fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% x01=1.2*xObj;
% [xRes1,objVal]=fminimax(@brDirectOptim,[x01 (omega2.*cAngle2)'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4 1e-4 1e-4  1e-4 1e-4 1e-4*ones(1,length(velData1))],[10000 10000  10000 10000  10000 10000  10000 10000  5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,v0,'All','multi','fixed');
% [xResOld1,objValOLD]=fmincon(@brGeneralOptim,x01,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 


x02=[10 0.1 0.1 10 [10 0.1 0.1 10 ]*1.09] ;
[xRes2,objVal]=fminimax(@brDirectOptim,[x02 (omega2.*cAngle2)'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4 1e-4 1e-4  1e-4 1e-4 1e-4*ones(1,length(velData1))],[10000 10000  10000 10000  10000 10000  10000 10000  5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,v0,'All','multi','fixed');
[xResOld2,objValOLD]=fmincon(@brGeneralOptim,x02,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');




x03=[50 0.01 0.01 50 50 0.01 0.01 50 ] ;
[xRes3,objVal]=fminimax(@brDirectOptim,[x03 (omega2.*cAngle2)'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4 1e-4 1e-4  1e-4 1e-4 1e-4*ones(1,length(velData1))],[10000 10000  10000 10000  10000 10000  10000 10000  5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,v0,'All','multi','fixed');
[xResOld3,objValOLD]=fmincon(@brGeneralOptim,x03,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');



x04=[50 0.01  0.1 10 50 0.01  0.1 10 ] ;
[xRes4,objVal]=fminimax(@brDirectOptim,[x04 (omega2.*cAngle2)'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4 1e-4 1e-4  1e-4 1e-4 1e-4*ones(1,length(velData1))],[10000 10000  10000 10000  10000 10000  10000 10000  5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta,v0,'All','multi','fixed');
[xResOld4,objValOLD]=fmincon(@brGeneralOptim,x04,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');

plotFct(xRes4,xResOld4,xObj,x04,'free','new',v0,velData1,velData2,cAngle1,cAngle2,'4');

plotFct(xRes1,xResOld1,xObj,x01,'free','new',v0,velData1,velData2,cAngle1,cAngle2,'1');
plotFct(xRes3,xResOld3,xObj,x03,'free','new',v0,velData1,velData2,cAngle1,cAngle2,'3');
plotFct(xRes2,xResOld2,xObj,x02,'free','new',v0,velData1,velData2,cAngle1,cAngle2,'2');