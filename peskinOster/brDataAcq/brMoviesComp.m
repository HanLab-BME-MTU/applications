function [oldVelMt1,oldVelMt2,velMt1,velMt2,oldCAngle1T,oldCAngle1TdT,oldCAngle2T,oldCAngle2TdT,cAngle1T,cAngle2T,cAngle1TdT,cAngle2TdT,sigmaMt1,sigmaMt2,stateMt1,stateMt2]=brMoviesComp
%BRCOMPMOVIE load and format the movies data

%SYNOPSIS
%[oldVelMt1,oldVelMt2,velMt1,velMt2,oldCAngle1T,oldCAngle1TdT,oldCAngle2T,oldCAngle2TdT,cAngle1T,cAngle2T,cAngle1TdT,cAngle2TdT,sigmaMt1,sigmaMt2,stateMt1,stateMt2]=brMoviesComp
%    
%
%INPUT 
%       NO
%
%OUTPUT
%       oldVelMt1         :  vector with velocity of first calculated without joans
%                            soft
%       oldVelMt2         :  vector with velocity of first calculated without joans
%                            soft
%       velMt1            :  vector with velocitywith joans soft.
%       velMt2            :  vector with velocitywith joans soft.
%       cAngle            :  vectors with cosine between Mt anf the
%                            centromere axis -old (without jonas soft) -T
%                            (projected at T) -TdT (projected at T+dT)
%       sigmaMTi          :  standart deviation for -MT1 and MT2
%       stateMTi          :  state of the MTi  

%path='H:\MatlabSave\'
path='/home/yloosli/MatlabSave/';

load([path 'WT_30C_Metaphase_Fast3_corr-data-05-Dec-2003-10-24-16.mat']);
idlisttrack3=idlisttrack_L;
dataProperties3=dataProperties;
lastResult3=lastResult;


load([path 'WT_30C_Metaphase_Fast5_corr-data-08-Dec-2003-11-36-23.mat']);
idlisttrack5=idlisttrack_L;
dataProperties5=dataProperties;
lastResult5=lastResult;

load([path 'WT_30C_Metaphase_Fast9_corr-data-08-Dec-2003-12-46-46.mat']);
idlisttrack9=idlisttrack_L;
dataProperties9=dataProperties;
lastResult9=lastResult;

load([path 'WT_30C_Metaphase_Fast11-data-05-Dec-2003-10-49-37.mat']);
idlisttrack11=idlisttrack_L;
dataProperties11=dataProperties;
lastResult11=lastResult;

load([path 'WT_30C_Metaphase_Fast18_corr-data-05-Dec-2003-11-45-46.mat']);
idlisttrack18=idlisttrack_L;
dataProperties18=dataProperties;
lastResult18=lastResult;

load([path 'RC203_2-5-04_30C_13_corr-data-17-Feb-2004-14-31-01.mat']);
idlisttrackRC13=idlisttrack_L;
dataPropertiesRC13=dataProperties;
lastResultRC13=lastResult;

load([path 'RC203_2-5-04_30C_14_corr-data-17-Feb-2004-14-45-09.mat']);
idlisttrackRC14=idlisttrack_L;
dataPropertiesRC14=dataProperties;
lastResultRC14=lastResult;

load([path 'RC203_2-4-04_30C_07_corr-data-12-Feb-2004-12-32-19.mat']);
idlisttrackRC07=idlisttrack_L;
dataPropertiesRC07=dataProperties;
lastResultRC07=lastResult;


[mt1State3,mt2State3,mt1Velocity3,mt2Velocity3,cAngle13,cAngle23,trajMatMt13,trajMatMt23,dataSpb3]=brGetAnaData(idlisttrack3 ,dataProperties3 ,lastResult3);
[mt1State5,mt2State5,mt1Velocity5,mt2Velocity5,cAngle15,cAngle25,trajMatMt15,trajMatMt25,dataSpb5]=brGetAnaData(idlisttrack5 ,dataProperties5 ,lastResult5);
[mt1State9,mt2State9,mt1Velocity9,mt2Velocity9,cAngle19,cAngle29,trajMatMt19,trajMatMt29,dataSpb9]=brGetAnaData(idlisttrack9 ,dataProperties9 ,lastResult9);
[mt1State11,mt2State11,mt1Velocity11,mt2Velocity11,cAngle111,cAngle211,trajMatMt111,trajMatMt211,dataSpb11]=brGetAnaData(idlisttrack11 ,dataProperties11 ,lastResult11);
[mt1State18,mt2State18,mt1Velocity18,mt2Velocity18,cAngle118,cAngle218,trajMatMt118,trajMatMt218,dataSpb18]=brGetAnaData(idlisttrack18 ,dataProperties18 ,lastResult18);
[mt1StateRC13,mt2StateRC13,mt1VelocityRC13,mt2VelocityRC13,cAngle1RC13,cAngle2RC13,trajMatMt1RC13,trajMatMt2RC13,dataSpbRC13]=brGetAnaData(idlisttrackRC13 ,dataPropertiesRC13 ,lastResultRC13);
[mt1StateRC14,mt2StateRC14,mt1VelocityRC14,mt2VelocityRC14,cAngle1RC14,cAngle2RC14,trajMatMt1RC14,trajMatMt2RC14,dataSpbRC14]=brGetAnaData(idlisttrackRC14 ,dataPropertiesRC14 ,lastResultRC14);
[mt1StateRC07,mt2StateRC07,mt1VelocityRC07,mt2VelocityRC07,cAngle1RC07,cAngle2RC07,trajMatMt1RC07,trajMatMt2RC07,dataSpbRC07]=brGetAnaData(idlisttrackRC07 ,dataPropertiesRC07 ,lastResultRC07);

close all;

oldVelMt1=[mt1Velocity3 ;mt1Velocity5 ;mt1Velocity9 ;mt1Velocity11 ;mt1Velocity18 ;mt1VelocityRC13 ;mt1VelocityRC14 ;mt1VelocityRC07];
oldVelMt2=[mt2Velocity3 ;mt2Velocity5 ;mt2Velocity9 ;mt2Velocity11 ;mt2Velocity18 ;mt2VelocityRC13 ;mt2VelocityRC14 ;mt2VelocityRC07];

oldCAngle1T   =[cAngle13(1:end-1); cAngle15(1:end-1) ;cAngle19(1:end-1) ;cAngle111(1:end-1) ;cAngle118(1:end-1) ;cAngle1RC13(1:end-1) ;cAngle1RC14(1:end-1) ;cAngle1RC07(1:end-1)];
oldCAngle1TdT =[cAngle13(2:end)  ; cAngle15(2:end)   ;cAngle19(2:end)   ;cAngle111(2:end)   ;cAngle118(2:end)   ;cAngle1RC13(2:end)   ;cAngle1RC14(2:end)  ;cAngle1RC07(2:end)  ];

oldCAngle2T   =[cAngle23(1:end-1) ;cAngle25(1:end-1) ;cAngle29(1:end-1) ;cAngle211(1:end-1) ;cAngle218(1:end-1) ;cAngle2RC13(1:end-1) ;cAngle2RC14(1:end-1) ;cAngle2RC07(1:end-1)];
oldCAngle2TdT =[cAngle23(2:end)   ;cAngle25(2:end)   ;cAngle29(2:end)   ;cAngle211(2:end)   ;cAngle218(2:end)   ;cAngle2RC13(2:end)   ;cAngle2RC14(2:end) ;cAngle2RC07(2:end)];



velMt1=[trajMatMt13(:,4);trajMatMt15(:,4);trajMatMt19(:,4);trajMatMt111(:,4);trajMatMt118(:,4);trajMatMt1RC13(:,4);trajMatMt1RC14(:,4);trajMatMt1RC07(:,4)]; 
velMt2=[trajMatMt23(:,4);trajMatMt25(:,4);trajMatMt29(:,4);trajMatMt211(:,4);trajMatMt218(:,4);trajMatMt2RC13(:,4);trajMatMt2RC14(:,4);trajMatMt2RC07(:,4)]; 

cAngle1T=[cAngle13(trajMatMt13(:,1)) ; cAngle15(trajMatMt15(:,1)) ; cAngle19(trajMatMt19(:,1)) ; cAngle111(trajMatMt111(:,1)) ; cAngle118(trajMatMt118(:,1)) ; cAngle1RC13(trajMatMt1RC13(:,1))  ; cAngle1RC14(trajMatMt1RC14(:,1)) ; cAngle1RC07(trajMatMt1RC07(:,1))];
cAngle2T=[cAngle23(trajMatMt23(:,1)) ; cAngle25(trajMatMt25(:,1)) ; cAngle29(trajMatMt29(:,1)) ; cAngle211(trajMatMt211(:,1)) ; cAngle218(trajMatMt218(:,1)) ; cAngle2RC13(trajMatMt2RC13(:,1))  ; cAngle2RC14(trajMatMt2RC14(:,1))  ; cAngle2RC07(trajMatMt2RC07(:,1))];

cAngle1TdT=[cAngle13(trajMatMt13(:,2)) ; cAngle15(trajMatMt15(:,2)) ; cAngle19(trajMatMt19(:,2)) ; cAngle111(trajMatMt111(:,2)) ; cAngle118(trajMatMt118(:,2)) ; cAngle1RC13(trajMatMt1RC13(:,2)) ; cAngle1RC14(trajMatMt1RC14(:,2)) ; cAngle1RC07(trajMatMt1RC07(:,2))];
cAngle2TdT=[cAngle23(trajMatMt23(:,2)) ; cAngle25(trajMatMt25(:,2)) ; cAngle29(trajMatMt29(:,2)) ; cAngle211(trajMatMt211(:,2)) ; cAngle218(trajMatMt218(:,2)) ; cAngle2RC13(trajMatMt2RC13(:,2)) ; cAngle2RC14(trajMatMt2RC14(:,2)) ; cAngle2RC07(trajMatMt2RC07(:,2))];

sigmaMt1=[trajMatMt13(:,5) ; trajMatMt15(:,5) ; trajMatMt19(:,5) ; trajMatMt111(:,5);trajMatMt118(:,5) ; trajMatMt1RC13(:,5);trajMatMt1RC14(:,5);trajMatMt1RC07(:,5)];
sigmaMt2=[trajMatMt23(:,5) ; trajMatMt25(:,5) ; trajMatMt29(:,5) ; trajMatMt211(:,5);trajMatMt218(:,5) ; trajMatMt2RC13(:,5);trajMatMt2RC14(:,5);trajMatMt2RC07(:,5)];

stateMt1=[trajMatMt13(:,3) ; trajMatMt15(:,3) ; trajMatMt19(:,3) ; trajMatMt111(:,3);trajMatMt118(:,3) ; trajMatMt1RC13(:,3);trajMatMt1RC14(:,3);trajMatMt1RC07(:,3)];
stateMt2=[trajMatMt23(:,3) ; trajMatMt25(:,3) ; trajMatMt29(:,3) ; trajMatMt211(:,3);trajMatMt218(:,3) ; trajMatMt2RC13(:,3);trajMatMt2RC14(:,3);trajMatMt2RC07(:,3)];

timeMt1=[trajMatMt13(:,6) ; trajMatMt15(:,6) ; trajMatMt19(:,6) ; trajMatMt111(:,6);trajMatMt118(:,6) ; trajMatMt1RC13(:,6);trajMatMt1RC14(:,6);trajMatMt1RC07(:,6)];
timeMt2=[trajMatMt23(:,6) ; trajMatMt25(:,6) ; trajMatMt29(:,6) ; trajMatMt211(:,6);trajMatMt218(:,6) ; trajMatMt2RC13(:,6);trajMatMt2RC14(:,6);trajMatMt2RC07(:,6)];

save([path 'data5SecFull2'],'oldVelMt1','oldVelMt2','velMt1','velMt2','oldCAngle1T','oldCAngle1TdT','oldCAngle2T','oldCAngle2TdT','cAngle1T','cAngle2T','cAngle1TdT','cAngle2TdT','sigmaMt1','sigmaMt2','stateMt1','stateMt2');

sigma1Global=sigmaMt1./(timeMt1.^(0.5));
sigma2Global=sigmaMt2./(timeMt2.^(0.5));

[weightedMean,weightedStd]=brVelocityAnalyze(abs(velMt1),sigma1Global,sigmaMt1,1e-3,'Mt1 all cases','k');

[weightedMean,weightedStd]=brVelocityAnalyze(abs(velMt2),sigma2Global,sigmaMt2,1e-3,'Mt2 all cases','k');


qualityRatio1=abs(sigmaMt1./velMt1);
qualityRatio2=abs(sigmaMt2./velMt2);
indexQualityRatioOk1=find(qualityRatio1<0.45);
indexQualityRatioOk2=find(qualityRatio2<0.45);

indexOK=intersect(indexQualityRatioOk1 , indexQualityRatioOk2);

[weightedMeanMt1,weightedStdMt1]=brVelocityAnalyze(velMt1(indexOK),sigma1Global(indexOK),sigmaMt1(indexOK),1e-3,'Mt1 all cases with discarding the \sigma/v > 0.45','k');

[weightedMeanMt2,weightedStdMt2]=brVelocityAnalyze(velMt2(indexOK),sigma2Global(indexOK),sigmaMt2(indexOK),1e-3,'Mt2 all cases with discarding the \sigma/v > 0.45','k');

% both growing

indexGG=intersect(find(velMt1>0&velMt2>0),indexOK);
[weightedMeanMt2,weightedStdMt2]=brVelocityAnalyze(velMt2(indexGG),sigma2Global(indexGG),sigmaMt2(indexGG),1e-3,'Mt2 GG with discarding the \sigma/v > 0.45','k');
[weightedMeanMt1,weightedStdMt1]=brVelocityAnalyze(velMt1(indexGG),sigma1Global(indexGG),sigmaMt1(indexGG),1e-3,'Mt1 GG with discarding the \sigma/v > 0.45','k');


save([path 'data5SecFull2'],'oldVelMt1','oldVelMt2','velMt1','velMt2','oldCAngle1T','oldCAngle1TdT','oldCAngle2T','oldCAngle2TdT','cAngle1T','cAngle2T','cAngle1TdT','cAngle2TdT','sigmaMt1','sigmaMt2','stateMt1','stateMt2','indexOK');
