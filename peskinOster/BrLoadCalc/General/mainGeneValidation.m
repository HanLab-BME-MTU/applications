% Starting guess evalution

% main Test
clear;
close all;

delta=8e-3;

load('velDataGeneral.mat');
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

p2=(beta2-(gamma2*delta*beta2)/vUnload)/gamma2+1;
p1=(beta1-(gamma1*delta*beta1)/vUnload)/gamma1+1;

xObj=[alpha1 betaG1 beta1 gamma1  alpha2 betaG2 beta2 gamma2 ];
fct=brGeneralOptim(xObj,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
% x0=[10 1 1 20 10 1 1 20 ];
% 
% [c,ceq]=brGeneralOptimConst(xObj,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% fct=brGeneralOptim(xObj,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
% opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',300,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);
% [xRes,objVal,exitFlag,outPut]=fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
% 
% xRes=[xRes(1:4) p1 xRes(5:8) p2];
% x0=[x0(1:4) p1 x0(5:8) p2];
% 
% %Plot Section
% 
% xObj=[alpha1 betaG1 beta1 gamma1 p1 alpha2 betaG2 beta2 gamma2 p2];
% 
% omega=[0:0.01:5];
% figure('Name','Growing')
% 
% 
% plot(omega,peskin(xObj(1:2),omega,delta),'--r');hold on;
% plot(omega,peskin(xObj(6:7),omega,delta),'--k');
% 
% plot(omega,peskin(xRes(1:2),omega,delta),'r');
% plot(omega,peskin(xRes(6:7),omega,delta),'k');
% 
% 
% plot(omega,peskin(x0(1:2),omega,delta),':r');
% plot(omega,peskin(x0(6:7),omega,delta),':k');
% 
% plot(omega1(indexGG),velData1(indexGG),'d');
% plot(omega2(indexGG),velData2(indexGG),'d');
% plot(omega1(indexGS),velData1(indexGS),'d');
% plot(omega2(indexSG),velData2(indexSG),'d');
% 
% 
% figure('Name','Shrinking')
% 
% plot(omega,brShrinkDirect(xRes(3:5),omega,delta),'r');hold on
% plot(omega,brShrinkDirect(xRes(8:10),omega,delta),'k');
% 
% plot(omega,brShrinkDirect(xObj(3:5),omega,delta),'--r');
% plot(omega,brShrinkDirect(xObj(8:10),omega,delta),'--k');
% 
% plot(omega,brShrinkDirect(x0(3:5),omega,delta),':r');
% plot(omega,brShrinkDirect(x0(8:10),omega,delta),':k');
% 
% plot(omega1(indexSS),-velData1(indexSS),'d');
% plot(omega2(indexSS),-velData2(indexSS),'d');
% plot(omega2(indexGS),-velData2(indexGS),'d');
% plot(omega1(indexSG),-velData1(indexSG),'d');
% 



sGuess1=[9.0287    0.148    0.0363    6.0387    7.0165    0.090    1.2829    5.4725];
sGuess2=[ 7.0165    0.90    4.2829    5.4725 9.0287    0.0148    0.363    6.0387 ];
sGuess3=1.2*sGuess1;
sGuess4=0.8*sGuess1;
sGuess5=2*sGuess1;
sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];

a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1)];
b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2)];
c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3)];
d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4)];
e=[sGuess(1,5) sGuess(2,5) sGuess(3,5) sGuess(4,5) sGuess(5,5)];
f=[sGuess(1,6) sGuess(2,6) sGuess(3,6) sGuess(4,6) sGuess(5,6)];
g=[sGuess(1,7) sGuess(2,7) sGuess(3,7) sGuess(4,7) sGuess(5,7)];
h=[sGuess(1,8) sGuess(2,8) sGuess(3,8) sGuess(4,8) sGuess(5,8)];
for i=1:5
    guessTemp(1,1)=a(i);
    for j=1:5
        guessTemp(1,2)=b(j);
        for k=1:5
            guessTemp(1,3)=c(k);
            for l=1:5
                guessTemp(1,4)=d(l);
                for m=1:5
                    guessTemp(1,5)=e(m);
                    for n=1:5
                        guessTemp(1,6)=f(n); 
                          for o=1:5
                              guessTemp(1,7)=g(o);
                              for p=1:5
                                  guessTemp(1,8)=h(p);
                                  guess(p+(o-1)*2+(n-1)*4+(m-1)*8+(l-1)*16+(k-1)*32+64*(j-1)+128*(i-1),:)=guessTemp;
                              end
                          end
                      end
                  end
              end
          end
      end
  end
end


load('guess.mat');

fid=fopen('H:\matlab\brModelAnal\brLoadCalc\General\guessSeekShrink1.txt','a');
fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
fclose(fid);
for i=[1:1]
    x0=guess(i,:);
    x0=[10 0.1 0.1 10 10 0.1 0.1 10];
    opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
    [xRes(i,:),objVal(i),exitFlag(i),outPut(i)]=fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
    fid=fopen('H:\matlab\brModelAnal\brLoadCalc\General\guessSeekShrink1.txt','a');
    fprintf(fid,'\n x0: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',x0);
	fprintf(fid,'\n xRes: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',xRes(i,:));
    fprintf(fid,'\n objVal: %12.8f ',objVal(i));
    fprintf(fid,'\n Output: ');
    fprintf(fid,'\n iteration : %d ',outPut(i).iterations);
    fprintf(fid,'\n funcCount : %d ',outPut(i).funcCount);
    fprintf(fid,'\n stepsize : %d ',outPut(i).stepsize); 
	fprintf(fid,'\n--------------------------------------------------end of %d ----------------------------------------------------',i);
	fclose(fid);
end

fid=fopen('H:\matlab\brModelAnal\brLoadCalc\General\guessSeekShrink1.txt','a');
fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
fclose(fid);

save('resultats2.mat','xRes','objVal');