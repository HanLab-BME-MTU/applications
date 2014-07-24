function [velData1,velData2,velDataNoise51,velDataNoise52,velDataNoise451,velDataNoise452,cAngle1,cAngle2,omega1,omega2,v0,xObj]=DataGenerationShrink(dataSize)% DataGeneration fro a validation

delta=8e-3;


test=0;

while test==0
	omega1=5*rand(10*dataSize,1);
	cAngle1=cos(0.01+(pi/2-0.01-0.01)*rand(10*dataSize,1));
	cAngle2=cos(0.01+(pi/2-0.01-0.01)*rand(10*dataSize,1));
	
	omega2=(omega1.*cAngle1)./cAngle2;
    testNeg=find(omega2>0);
    testTooBig=find(omega2<5);
    if length(testNeg)>dataSize & length(testTooBig)>dataSize;
        test=1;
    end
end

omega1=omega1(intersect(testNeg,testTooBig));
omega2=omega2(intersect(testNeg,testTooBig));
cAngle1=cAngle1(intersect(testNeg,testTooBig));
cAngle2=cAngle2(intersect(testNeg,testTooBig));

omega1=omega1(1:dataSize);
omega2=omega2(1:dataSize);
cAngle1=cAngle1(1:dataSize);
cAngle2=cAngle2(1:dataSize);

beta1=1.1;
gamma1=122.2;


beta2=0.12;
gamma2=100;

v0=0.08;

p1=(beta1-(gamma1*delta*beta1)/v0)/gamma1+1;
p2=(beta2-(gamma2*delta*beta2)/v0)/gamma2+1;



xObj = [beta1 gamma1 p1 beta2 gamma2 p2];
xObj2 = [beta1 gamma1 beta2 gamma2];

velData1=-brShrinkDirect(xObj(1:3),omega1,delta);
velData2=-brShrinkDirect(xObj(4:6),omega2,delta);

velDataNoise451=velData1.*(0.55+(1.45-0.55)*rand(size(velData1)));
velDataNoise452=velData2.*(0.55+(1.45-0.55)*rand(size(velData2)));

velDataNoise51=velData1.*(0.95+(1.05-0.95)*rand(size(velData1)));
velDataNoise52=velData2.*(0.95+(1.05-0.95)*rand(size(velData2)));

vUnload=v0;
close all;
figure
plot(omega1,-velData1,'d');hold on;
plot(omega2,-velData2,'d');

save('DataShrink.mat','velData1','velData2','velDataNoise451','velDataNoise452','velDataNoise51','velDataNoise52','cAngle1','cAngle2','omega1','omega2','vUnload','xObj','xObj2');
