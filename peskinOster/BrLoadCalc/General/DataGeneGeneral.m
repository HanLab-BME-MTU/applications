% Data generation general case

clear;

delta=8e-3;
v0=0.02;


alpha1=9;
betaG1=0.02;
beta1=0.1 ; 
gamma1=12;

alpha2=7;
betaG2=0.012;
beta2= 0.6;
gamma2=9.08;


p2=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;
p1=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;

xObj=[alpha1 betaG1 beta1 gamma1 p1 alpha2 betaG2 beta2 gamma2 p2];

omega1=0.1+(4-0.1)*rand(100,1);
omega2=0.1+(4-0.1)*rand(100,1);


cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(100,1));
cAngle2=(omega1.*cAngle1)./omega2;
% cAngle1=ones(300,1);
% cAngle2=ones(300,1);
% omega2=omega1;


nbData=100;

index=randperm(nbData);


indexGG=index(1:12);
indexSS=index(13:26);
indexSG=index(27:65);
indexGS=index(66:100);


velData1(indexGG)=peskin(xObj(1:2),omega1(indexGG),delta);
velData2(indexGG)=peskin(xObj(6:7),omega2(indexGG),delta);

velData1(indexSS)=-brShrinkDirect(xObj(3:5),omega1(indexSS),delta);
velData2(indexSS)=-brShrinkDirect(xObj(8:10),omega2(indexSS),delta);

velData1(indexGS)=peskin(xObj(1:2),omega1(indexGS),delta);
velData2(indexGS)=-brShrinkDirect(xObj(8:10),omega2(indexGS),delta);

velData1(indexSG)=-brShrinkDirect(xObj(3:5),omega1(indexSG),delta);
velData2(indexSG)=peskin(xObj(6:7),omega2(indexSG),delta);

vUnload=v0;
save('velDataGeneral','velData1','velData2','omega1','omega2','indexSS','indexGG','indexSG','indexGS','cAngle1','cAngle2','vUnload','xObj');


