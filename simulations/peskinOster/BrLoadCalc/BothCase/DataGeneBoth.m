%Data generation

clear;

delta=8e-3;

omega1=5*rand(200,1);
omega2=5*rand(200,1);

cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(200,1));
cAngle2=(omega1.*cAngle1)./omega2;
omega2=omega1;
cAngle1=ones(200,1);
cAngle2=ones(200,1);

v0=4.2;

nbGS=100;

beta1=47; 
gamma1=2165;


beta2= 69;
gamma2=1934;
pGS=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;

pSG=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;
% GS


xObjGS=[160 0.4 69 1934 pGS];     % alphaG betaG betaS gammaS pS
velDataGS1=peskin(xObjGS(1:2),omega1(1:nbGS),delta);
velDataGS2=-brShrinkDirect(xObjGS(3:5),omega2(1:nbGS),delta);



% SG
xObjSG=[47 2165 pSG 201 0.05];
velDataSG1=-brShrinkDirect(xObjSG(1:3),omega1(nbGS+1:end),delta);
velDataSG2=peskin(xObjSG(4:5),omega2(nbGS+1:end),delta);


vUnload=v0;
save('DataBoth','velDataGS1','velDataGS2','velDataSG1','velDataSG2','cAngle1','cAngle2','omega1','omega2','vUnload');
