% main G1 Velocity analysis


list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_01\WT_30C_G1_noBen_01_corr-data-23-Jul-2003-16-49-09.mat');
[trajectoryDescription1,data1,mtVelocity1,frameNb1,mtVelCons1,frameCons1]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_02\WT_30C_G1_noBen_02_corr-data-23-Jul-2003-16-56-24.mat');
[trajectoryDescription2,data2,mtVelocity2,frameNb2,mtVelCons2,frameCons2]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_03\WT_30C_G1_noBen_03_corr-data-23-Jul-2003-17-04-36.mat');
[trajectoryDescription3,data3,mtVelocity3,frameNb3,mtVelCons3,frameCons3]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_04\WT_30C_G1_noBen_04_corr-data-23-Jul-2003-17-11-51.mat');
[trajectoryDescription4,data4,mtVelocity4,frameNb4,mtVelCons4,frameCons4]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_05\WT_30C_G1_noBen_05_corr-data-23-Jul-2003-17-19-32.mat');
[trajectoryDescription5,data5,mtVelocity5,frameNb5,mtVelCons5,frameCons5]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_06\WT_30C_G1_noBen_06_corr-data-23-Jul-2003-17-27-06.mat');
[trajectoryDescription6,data6,mtVelocity6,frameNb6,mtVelCons6,frameCons6]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_07\WT_30C_G1_noBen_07_corr-data-23-Jul-2003-17-34-33.mat');
[trajectoryDescription7,data7,mtVelocity7,frameNb7,mtVelCons7,frameCons7]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_08\WT_30C_G1_noBen_08_corr-data-23-Jul-2003-17-42-26.mat');
[trajectoryDescription8,data8,mtVelocity8,frameNb8,mtVelCons8,frameCons8]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

list=open('M:\yeastData\Wildtype30C\5secMovies\WT_30C_G1_noBen_09\WT_30C_G1_noBen_09_corr-data-23-Jul-2003-17-50-06.mat');
[trajectoryDescription9,data9,mtVelocity9,frameNb9,mtVelCons9,frameCons9]=brGetAnaDataG1(list.idlisttrack_L,list.dataProperties,list.lastResult);

% Extracting the needed value

state=[trajectoryDescription1.individualStatistics.dataListGroup(:,3);trajectoryDescription2.individualStatistics.dataListGroup(:,3);trajectoryDescription3.individualStatistics.dataListGroup(:,3);trajectoryDescription4.individualStatistics.dataListGroup(:,3);trajectoryDescription5.individualStatistics.dataListGroup(:,3)];
state=[state;trajectoryDescription6.individualStatistics.dataListGroup(:,3);trajectoryDescription7.individualStatistics.dataListGroup(:,3);trajectoryDescription8.individualStatistics.dataListGroup(:,3);trajectoryDescription9.individualStatistics.dataListGroup(:,3)];

velocity=[trajectoryDescription1.individualStatistics.dataListGroup(:,4);trajectoryDescription2.individualStatistics.dataListGroup(:,4);trajectoryDescription3.individualStatistics.dataListGroup(:,4);trajectoryDescription4.individualStatistics.dataListGroup(:,4);trajectoryDescription5.individualStatistics.dataListGroup(:,4)];
velocity=[velocity;trajectoryDescription6.individualStatistics.dataListGroup(:,4);trajectoryDescription7.individualStatistics.dataListGroup(:,4);trajectoryDescription8.individualStatistics.dataListGroup(:,4);trajectoryDescription9.individualStatistics.dataListGroup(:,4)];

sigma=[trajectoryDescription1.individualStatistics.dataListGroup(:,5);trajectoryDescription2.individualStatistics.dataListGroup(:,5);trajectoryDescription3.individualStatistics.dataListGroup(:,5);trajectoryDescription4.individualStatistics.dataListGroup(:,5);trajectoryDescription5.individualStatistics.dataListGroup(:,5)];
sigma=[sigma;trajectoryDescription6.individualStatistics.dataListGroup(:,5);trajectoryDescription7.individualStatistics.dataListGroup(:,5);trajectoryDescription8.individualStatistics.dataListGroup(:,5);trajectoryDescription9.individualStatistics.dataListGroup(:,5)];

time=[trajectoryDescription1.individualStatistics.dataListGroup(:,7);trajectoryDescription2.individualStatistics.dataListGroup(:,7);trajectoryDescription3.individualStatistics.dataListGroup(:,7);trajectoryDescription4.individualStatistics.dataListGroup(:,7);trajectoryDescription5.individualStatistics.dataListGroup(:,7)];
time=[time;trajectoryDescription6.individualStatistics.dataListGroup(:,7);trajectoryDescription7.individualStatistics.dataListGroup(:,7);trajectoryDescription8.individualStatistics.dataListGroup(:,7);trajectoryDescription9.individualStatistics.dataListGroup(:,7)];


close all;
% Looking for the shrinkage

stateSIndex=find(state==2);
stateGIndex=find(state==1);

velS   = abs(velocity(stateSIndex));
sigmaS = sigma(stateSIndex);
timeS=time(stateSIndex);

sigmaGS=sigmaS./(timeS.^(0.5));

dataVel=[velS,0.001*ones(size(velS))];

[velMeanS,weightedStd]=brVelocityAnalyze(velS,sigmaGS,sigmaS,1e-3,'G1 Shrinkage case','o');

velG = velocity(stateGIndex);
sigmaG =sigma(stateGIndex);
timeG=time(stateGIndex);

sigmaGG=sigmaG./(timeG.^(0.5));

[velMeanG,weightedStd]=brVelocityAnalyze(velG,sigmaGG,sigmaG,1e-3,'G1 growing case','o');
 






