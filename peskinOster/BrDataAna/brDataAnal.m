% -----------------brDataAnal

clear;close all;

load('velocityG1.mat');




indexGG0=find(mt1Vel>0 & mt2Vel>0);
indexSS0=find(mt1Vel<0 & mt2Vel<0);
indexGS0=find(mt1Vel>0 & mt2Vel<0);
indexSG0=find(mt1Vel<0 & mt2Vel>0);

mt0Index=[1:length(mt1Vel)];


%---------GENRAL DSITRIBUATION

xG=[0:5e-5:0.0025];
figure('Name','Cumulative distribution - Growing');
% nG=histc([mt1Vel(indexGG0) mt2Vel(indexGG0) mt1Vel(indexGS0) mt2Vel(indexSG0)],xG);
% bar(nG);

maxG=max([mt1Vel(indexGG0) mt2Vel(indexGG0) mt1Vel(indexGS0) mt2Vel(indexSG0)]);
xG2=[0:5e-5:maxG+5e-4];

for i=1:length(xG2)
    a(i)=length(find([mt1Vel(indexGG0) mt2Vel(indexGG0) mt1Vel(indexGS0) mt2Vel(indexSG0)]<xG2(i)));
end
plot(xG2,a);

% lower case
xG3=[0:1e-6:5e-4];

for i=1:length(xG3)
    a2(i)=length(find([mt1Vel(indexGG0) mt2Vel(indexGG0) mt1Vel(indexGS0) mt2Vel(indexSG0)]<xG3(i)));
end
plot(xG3,a2);




xS=[0:5e-5:0.036];
figure('Name','Cumulative distribution -Shrinking');
% nS=histc(abs([mt1Vel(indexSS0) mt1Vel(indexSG0) mt2Vel(indexSS0) mt2Vel(indexGS0) ]),xS);
% bar(nS);
maxS=max(abs([mt1Vel(indexSS0) mt1Vel(indexSG0) mt2Vel(indexSS0) mt2Vel(indexGS0)]));
xS2=[0:5e-5:maxS+5e-5];

for i=1:length(xS2)
    b(i)=length(find(abs([mt1Vel(indexSS0) mt1Vel(indexSG0) mt2Vel(indexSS0) mt2Vel(indexGS0)])<xS2(i)));
end
plot(xS2,b);

% lower case
xS3=[0:1e-6:5e-4];

for i=1:length(xS3)
    b2(i)=length(find(abs([mt1Vel(indexSS0) mt1Vel(indexSG0) mt2Vel(indexSS0) mt2Vel(indexGS0)])<xS3(i)));
end
plot(xS3,b2);


%---First discarding step  The upper limit G=0.02 S=0.02

indexGG=find(mt1Vel>0 & mt2Vel>0 & mt1Vel<0.02 & mt2Vel<0.02);
indexSS=find(mt1Vel<0 & mt2Vel<0 & mt1Vel<0.02 & mt2Vel<0.02);
indexGS=find(mt1Vel>0 & mt2Vel<0 & mt1Vel<0.02 & mt2Vel<0.02);
indexSG=find(mt1Vel<0 & mt2Vel>0 & mt1Vel<0.02 & mt2Vel<0.02);

velData1=mt1Vel([indexGG indexSS indexGS indexSG]);
velData2=mt2Vel([indexGG indexSS indexGS indexSG]);

discardRatio1=1-length(velData1)/length(mt1Vel);
% d
% % Discard the velocity lower than 1e-6 fro s cases
% indexGG=find(mt1Vel>0 & mt2Vel>0 & mt1Vel<0.0025 & mt2Vel<0.0025 );
% indexSS=find(mt1Vel<0 & mt2Vel<0 & abs(mt1Vel)<0.002 & abs(mt2Vel)<0.002 & abs(mt1Vel)>1e-6 & abs(mt2Vel)>1e-6);
% indexGS=find(mt1Vel>0 & mt2Vel<0 & mt1Vel<0.0025 & abs(mt2Vel)<0.002 & abs(mt2Vel)>1e-6);
% indexSG=find(mt1Vel<0 & mt2Vel>0 & mt1Vel<0.002 & abs(mt2Vel)<0.0025 & abs(mt1Vel)>1e-6 );
% mt1Index=[indexGG indexSS indexGS indexSG];
% mt2Index=[indexGG indexSS indexGS indexSG];
% velData1=mt1Vel(mt1Index);
% velData2=mt2Vel(mt2Index);
% 
% discardRatio2=1-length(velData1)/length(mt1Vel);
% 
% 
% 
% close all;
% % 
% % 
% delta=8e-3;
% opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',300,'TolFun',1e-5,'TolX',1e-5,'TolCon',1e-5);
% 
% x01=[0.1 0.01 0.2 0.50 0.2 0.15  0.01 0.20 1 0.4];
% x02=[1  0.1 0.1  5 0.1 0.02 0.1  1 11 0.1] ;
% x03=[20 0.01 0.2 28 0.3 2 0.01 1 25 0.2];
% x04=[0 0 0 0 0 0 0 0 0 0];
% [xRes1,objVal1,a1,b1]=fmincon(@brGeneralOptim,x01,[],[],[],[],[1e-6 1e-6 1e-6 1e-6 0.1 1e-6 1e-6 1e-6 1e-6 0.1],[1000 10 100 3000 0.9 1000 10 100 3000 0.9  ],@brGeneralOptimConst,opts,velData1,velData2,delta);
% save('resultat1','xRes1','objVal1');
% [xRes2,objVal2,a2,b2]=fmincon(@brGeneralOptim,x02,[],[],[],[],[1e-6 1e-6 1e-6 1e-6 0.1 1e-6 1e-6 1e-6 1e-6 0.1],[1000 10 100 3000 0.9 1000 10 100 3000 0.9  ],@brGeneralOptimConst,opts,velData1,velData2,delta);
% save('resultat2','xRes2','objVal2');
% [xRes3,objVal3,a3,b3]=fmincon(@brGeneralOptim,x03,[],[],[],[],[1e-6 1e-6 1e-6 1e-6 0.1 1e-6 1e-6 1e-6 1e-6 0.1],[1000 10 100 3000 0.9 1000 10 100 3000 0.9  ],@brGeneralOptimConst,opts,velData1,velData2,delta);
% save('resultat3','xRes3','objVal3');
% [xRes4,objVal4,a4,b4]=fmincon(@brGeneralOptim,x04,[],[],[],[],[1e-6 1e-6 1e-6 1e-6 0.1 1e-6 1e-6 1e-6 1e-6 0.1],[1000 10 100 3000 0.9 1000 10 100 3000 0.9  ],@brGeneralOptimConst,opts,velData1,velData2,delta);
% save('resultat4','xRes4','objVal4');
