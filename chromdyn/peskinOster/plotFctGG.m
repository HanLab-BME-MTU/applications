function plotFctGG(xRes,xObj,velData1,velData2,cAngle1,cAngle2,delta,name)


omegaRes1=xRes(5:end)./cAngle1';
omegaRes2=xRes(5:end)./cAngle2';

alpha1=xRes(1);
betaG1=xRes(2);


alpha2=xRes(3);
betaG2=xRes(4);



xRes= [alpha1 betaG1  alpha2 betaG2 ];

alpha1=xObj(1);
betaG1=xObj(2);


alpha2=xObj(3);
betaG2=xObj(4);




xObj= [alpha1 betaG1   alpha2 betaG2];




omega=[0:0.01:8];

 nom=['Mt1, ',name];
figure('Name',nom)
plot(omega,peskin(xObj(1:2),omega,delta),'--r');hold on;
plot(omega,peskin(xRes(1:2),omega,delta),'r');

plot(omegaRes1,velData1,'dr');

nom=['Mt2, ' name]
figure('Name',nom)
plot(omega,peskin(xObj(3:4),omega,delta),'--k');hold on;
plot(omega,peskin(xRes(3:4),omega,delta),'k')

plot(omegaRes2,velData2,'dr');