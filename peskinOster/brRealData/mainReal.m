% main analyze real data
clear;

path='h:\MatlabSave\';

load([path 'data5SecFull2']);

delta=8e-3;

indexGG=find(velMt1>0 & velMt2>0);

velData1=velMt1(indexGG);
velData2=velMt2(indexGG);

cAngle1=cAngle1T(indexGG) ;
cAngle2=cAngle2T(indexGG) ;

velData1=velData1([1:8 10:end]);
velData2=velData2([1:8 10:end]);

cAngle1=cAngle1([1:8 10:end]) ;
cAngle2=cAngle2([1:8 10:end]) ;


opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);

x0=[0 0 0 0 ];

[xResGGT,objValT,exitFlag,outPut]=fminimax(@brBothGrowingOptimFct,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13);


[omega11,omega21,omega1Eff,omega2Eff,omegaDiff]=brOmegaGene(xResGGT,velData1,velData2,cAngle1,cAngle2,[],delta/13,'GG')
opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-5,'TolX',1e-5,'TolCon',1e-5);

omega2Eff=(omega1Eff+omega2Eff)*0.5;
[xRes3,objVal]=fmincon(@brDirectOptim,[xResGGT omega2Eff'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1))],[10000 10000  10000 10000 5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13,[],'GG','single','free');

[xRes4,objVal]=fminimax(@brDirectOptim,[xResGGT omega2Eff'],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1))],[10000 10000  10000 10000 5*ones(1,length(velData1))],@brDirectConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13,[],'GG','multi','free');

% 
% 
% cAngle1=cAngle1TdT(indexGG) ;
% cAngle2=cAngle2TdT(indexGG) ;
% 
% velData1=velMt1(indexGG);
% velData2=velMt2(indexGG);
% 
% velData1=velData1([1:8 10:end]);
% velData2=velData2([1:8 10:end]);
% 
% cAngle1=cAngle1([1:8 10:end]) ;
% cAngle2=cAngle2([1:8 10:end]) ;
% 
% 
% [xResGGTdT,objValTdT,exitFlag,outPut]=fminimax(@brBothGrowingOptimFct,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 100 1000 100],@brBothGrowingOptimFctConst,opts,velData1,velData2,cAngle1,cAngle2,delta/13);
% 
% 
% %%%%-----------------plot section
omega=[0:0.01:5];

figure('Name','Growth Growth case for t projection');


plot(omega,peskin(xResGGT(1:2),omega,delta/13),'r');hold on;
plot(omega,peskin(xResGGT(3:4),omega,delta/13),'k');
omega1=xRes3(5:end)./cAngle1';
omega2=xRes3(5:end)./cAngle2';
plot(omega1,velData1,'dr');
plot(omega1,peskin(xRes3(1:2),omega1,delta/13),'or');

plot(omega2,velData2,'dk');
plot(omega2,peskin(xRes3(3:4),omega2,delta/13),'ok');


% 
% figure('Name','Growth Growth case for t+dt projection');
% 
% plot(omega,peskin(xResGGTdT(1:2),omega,delta/13),'r');hold on;
% plot(omega,peskin(xResGGTdT(3:4),omega,delta/13),'k');
% if ~isempty(find(peskinInv(velData1,xResGGTdT(1:2),delta/13)<0))|~isempty(find(peskinInv(velData2,xResGGTdT(3:4),delta/13)<0))
%     warning('there is some negative omega for t+dt');
% end
% 
% plot(abs(peskinInv(velData1,xResGGTdT(1:2),delta/13)),velData1,'dr');
% plot(peskinInv(velData2,xResGGTdT(3:4),delta/13),velData2,'dk');
% 
% indexSS=find(velMt1<0 &velMt2<0);
% 
% vUnload=0.0260;
% 
% velData1=velMt1(indexSS);
% velData2=velMt2(indexSS);
% 
% cAngle1=cAngle1T(indexSS) ;
% cAngle2=cAngle2T(indexSS) ;
% 
% opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
% 
% x0=[1 10 1 10];
% [xResSS,objVal]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[100 1000 100 1000 ],@brBothShrinkOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
% x0=[xResGGT(1:2) xResSS(1:2) xResGGT(3:4) xResSS(3:4)];
% 
% [xRes,objVal,exitFlag,outPut]=fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velMt1,velMt2,cAngle1T,cAngle2T,delta,vUnload,'fixed');
% sGuess1=[0     0.01        0     0.01      ];
% sGuess2=[6      12    6      12    ];
% sGuess3=[12      24    12      24    ];
% sGuess4=[18    30       18    30    ];
% sGuess5=[24    60       24    60    ];
% 
% 
% sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];
% 
% a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1)];
% b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2)];
% c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3)];
% d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4)];
% 
% for i=1:5
%     guessTemp(1,1)=a(i);
%     for j=1:5
%         guessTemp(1,2)=b(j);
%         for k=1:5
%             guessTemp(1,3)=c(k);
%             for l=1:5
%                 guessTemp(1,4)=d(l);
%                 guess(l+(k-1)*2+4*(j-1)+8*(i-1),:)=guessTemp;
% 
%             end
%         end
%     end
% end
% 
% 
% 
% fid=fopen('H:\MatlabSave\guessSeekShrink1.txt','a');
% fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
% fclose(fid);
% [1:3 5:6 8:11 13:14 16 ]
% for i=[1:61]
%     x0=guess(i,:);
%     opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
%     i
%     [xRes(i,:),objVal(i),exitFlag(i),outPut(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[100 1000 100 1000 ],@brBothShrinkOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
%     fid=fopen('H:\MatlabSave\guessSeekShrink1.txt','a');
%     fprintf(fid,'\n x0: 12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',x0);
% 	fprintf(fid,'\n xRes: 12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',xRes(i,:));
%     fprintf(fid,'\n objVal: 12.8f ',objVal(i));
%     fprintf(fid,'\n Output: ');
%     fprintf(fid,'\n iteration : d ',outPut(i).iterations);
%     fprintf(fid,'\n funcCount : d ',outPut(i).funcCount);
%     fprintf(fid,'\n stepsize : d ',outPut(i).stepsize); 
% 	fprintf(fid,'\n--------------------------------------------------end of d ----------------------------------------------------',i);
% 	fclose(fid);
% end
% 
% fid=fopen('H:\MatlabSave\guessSeekShrink1.txt','a');
% fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
% fclose(fid);
% 
% save('resultats2.mat','xRes','objVal');