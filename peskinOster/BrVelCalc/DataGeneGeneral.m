% Data generation general case

%clear;

delta=8e-3;


v0=0.9;

alpha1=80;
betaG1=0.02;
beta1=30; 
gamma1=291;

alpha2=382;
betaG2=0.9;
beta2= 20;
gamma2=428;

p2=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;
p1=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;

xObj=[alpha1 betaG1 beta1 gamma1 p1 alpha2 betaG2 beta2 gamma2 p2];

omega1=5*rand(250,1);
omega1(find(omega1==min(omega1)))=0.1;
omega2=5*rand(250,1);


cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(250,1));
cAngle2=(omega1.*cAngle1)./omega2;
% cAngle1=ones(300,1);
% cAngle2=ones(300,1);
% omega2=omega1;


nbData=250;

index=randperm(nbData);


indexGG=index(1:20);
indexSS=index(60:80);
indexSG=index(121:190);
indexGS=index(191:250);


velData1(indexGG)=peskin(xObj(1:2),omega1(indexGG),delta/13);
velData2(indexGG)=peskin(xObj(6:7),omega2(indexGG),delta/13);

velData1(indexSS)=-brShrinkDirect(xObj(3:5),omega1(indexSS),delta);
velData2(indexSS)=-brShrinkDirect(xObj(8:10),omega2(indexSS),delta);

velData1(indexGS)=peskin(xObj(1:2),omega1(indexGS),delta/13);
velData2(indexGS)=-brShrinkDirect(xObj(8:10),omega2(indexGS),delta);

velData1(indexSG)=-brShrinkDirect(xObj(3:5),omega1(indexSG),delta);
velData2(indexSG)=peskin(xObj(6:7),omega2(indexSG),delta/13);

vUnload=v0;

x=[xObj';omega1.*cAngle1];
x2=[xObj(1:4)'; xObj(6:9)' ;omega1.*cAngle1];
xGS=[xObj(1:2)'; xObj(8:9)'; omega1(indexGS).*cAngle1(indexGS)];
xSG=[xObj(3:4)'; xObj(6:7)'; omega1(indexSG).*cAngle1(indexSG)];
figure
plot(omega1(indexSS),-velData1(indexSS),'d');hold on;
plot([0:0.01:5],0,'r');

plot(omega2(indexSS),-velData2(indexSS),'dk');
plot([0:0.01:5],0,'r');hold off;

noise=0.1;

velDataNoise1=velData1.*(0.55+(1.45-0.55)*rand(size(velData1)));
velDataNoise2=velData2.*(0.55+(1.45-0.55)*rand(size(velData2)));

save('H:\MatlabSave\velDataGeneralGeneral','velData1','velData2','velDataNoise1','velDataNoise2','omega1','omega2','indexSS','indexGG','indexSG','indexGS','cAngle1','cAngle2','vUnload','x','x2','xGS','xSG','delta','xObj','v0');


