function [velData1,velData2,velDataNoise51,velDataNoise52,velDataNoise451,velDataNoise452,cAngle1,cAngle2,omega1,omega2,v01,v02,xObj]=DataGenerationGrow(dataSize);

delta=8e-3/13;

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

xObj=[25 0.15 45 0.3  ];
velData1=peskin(xObj(1:2),omega1,delta);
velData2=peskin(xObj(3:4),omega2,delta);
v01=peskin(xObj(1:2),0,delta);
v02=peskin(xObj(3:4),0,delta);
velDataNoise451=velData1.*(0.55+(1.45-0.55)*rand(size(velData1)));
velDataNoise452=velData2.*(0.55+(1.45-0.55)*rand(size(velData2)));
velDataNoise51=velData1.*(0.95+(1.05-0.95)*rand(size(velData1)));
velDataNoise52=velData2.*(0.95+(1.05-0.95)*rand(size(velData2)));
save('DataGrow','velData1','velData2','velDataNoise51','velDataNoise52','velDataNoise451','velDataNoise452','cAngle1','cAngle2','omega1','omega2','v01','v02','xObj');
% plot(omega1,velData1,'d');
% hold on;
% plot(omega2,velData2,'dr');
% plot([0:0.01:5],0,'r');
% 
% plot(omega1,velData1,'o');
% 
% plot(omega2,velData2,'or');
