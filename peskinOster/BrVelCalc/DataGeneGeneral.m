% Data generation general case

clear;

delta=8e-3;
v0=0.03;


alpha1=9;
betaG1=0.02;
beta1=1.1 ; 
gamma1=29;

alpha2=12;
betaG2=0.08;
beta2= 0.1;
gamma2=36;

p2=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;
p1=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;

xObj=[alpha1 betaG1 beta1 gamma1 p1 alpha2 betaG2 beta2 gamma2 p2];

omega1=5*rand(300,1);
omega2=5*rand(300,1);


cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(300,1));
cAngle2=(omega1.*cAngle1)./omega2;
% cAngle1=ones(300,1);
% cAngle2=ones(300,1);
% omega2=omega1;


nbData=300;

index=randperm(nbData);


indexGG=index(1:75);
indexSS=index(76:150);
indexSG=index(151:225);
indexGS=index(226:300);


velData1(indexGG)=peskin(xObj(1:2),omega1(indexGG),delta);
velData2(indexGG)=peskin(xObj(6:7),omega2(indexGG),delta);

velData1(indexSS)=-brShrinkDirect(xObj(3:5),omega1(indexSS),delta);
velData2(indexSS)=-brShrinkDirect(xObj(8:10),omega2(indexSS),delta);

velData1(indexGS)=peskin(xObj(1:2),omega1(indexGS),delta);
velData2(indexGS)=-brShrinkDirect(xObj(8:10),omega2(indexGS),delta);

velData1(indexSG)=-brShrinkDirect(xObj(3:5),omega1(indexSG),delta);
velData2(indexSG)=peskin(xObj(6:7),omega2(indexSG),delta);

vUnload=v0;

x=[xObj';omega1.*cAngle1];
x2=[xObj(1:4)'; xObj(6:9)' ;omega1.*cAngle1];
xGS=[xObj(1:2)'; xObj(8:9)'; omega1(indexGS).*cAngle1(indexGS)];
xSG=[xObj(3:4)'; xObj(6:7)'; omega1(indexSG).*cAngle1(indexSG)];

save('velDataGeneral','velData1','velData2','omega1','omega2','indexSS','indexGG','indexSG','indexGS','cAngle1','cAngle2','vUnload','x','x2','xGS','xSG','delta','xObj');


