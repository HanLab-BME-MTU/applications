% DataGeneration fro a validation
clear;
delta=8e-3;

v0=0.0241;

omega1=5*rand(80,1);
omega2=5*rand(80,1);

cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(80,1));
cAngle2=(omega1.*cAngle1)./omega2;

alpha1=6;
alpha2=9;

xObj=[alpha1 0.02 alpha2 0.01];

velData1=peskin(xObj(1:2),omega1,delta);
velData2=peskin(xObj(3:4),omega2,delta);

save('DataGrow','velData1','velData2','cAngle1','cAngle2','omega1','omega2');
