% DataGeneration fro a validation
clear;
delta=8e-3;

omega1=5*rand(150,1);
omega2=5*rand(150,1);

cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(150,1));
cAngle2=(omega1.*cAngle1)./omega2;
beta1=140;
gamma1=2415;


beta2=199;
gamma2=1794;

v0=4.2;

p1=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;
p2=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;



xObj = [140 2415 p1 199 1794 p2];

velData1=-brShrinkDirect(xObj(1:3),omega1,delta);
velData2=-brShrinkDirect(xObj(4:6),omega2,delta);

vUnload=v0;

save('DataShrink','velData1','velData2','cAngle1','cAngle2','omega1','omega2','vUnload');
