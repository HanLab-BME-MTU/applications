%  main function

clear;
delta=8e-3;
v0=0.03;

load('H:\MatlabSave\velDataGeneralGeneral');



opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',1500,'TolFun',1e-5,'TolX',1e-5,'TolCon',1e-5);

v0=0.9;

alpha1=420;
betaG1=0.2;
beta1=38; 
gamma1=281;

alpha2=382;
betaG2=0.9;
beta2= 84;
gamma2=528;




x0=[0 0 0 0 ];


xObj=[alpha1 betaG1 beta1 gamma1 alpha2 betaG2 beta2 gamma2];

[omega11,omega21,omega1Eff,omega2Eff,omegaDiff]=brOmegaGene([alpha1 betaG1 alpha2 betaG2],velData1(indexGG),velData2(indexGG),cAngle1(indexGG),cAngle2(indexGG),v0,delta,'GG');
%[xResGG,objVal]=fmincon(@brDirectOptim,[x0 ones(1,30)],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1(indexGG)))],[10000 10000   10000 10000  5*ones(1,length(velData1(indexGG)))],@brDirectConst,opts,velData1(indexGG),velData2(indexGG),cAngle1(indexGG),cAngle2(indexGG),delta,v0,'GG','single','fixed');

% [xResSS,objVal]=fmincon(@brDirectOptim,[60    600        60    600 ones(1,30)],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1(indexSS)))],[10000 10000   10000 10000  5*ones(1,length(velData1(indexSS)))],@brDirectConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,v0,'SS','single','fixed');

% [xResSS,objVal]=fminimax(@brDirectOptim,[60    600        60    600 ones(1,30)],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1(indexSS)))],[10000 10000   10000 10000  5*ones(1,length(velData1(indexSS)))],@brDirectConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,v0,'SS','multi','fixed');

sGuess1=[1     1        1     1     ];
sGuess2=[20   100      20     100     ];
sGuess3=[40    300      40    300   ];
sGuess4=[60    600        60    600     ];
sGuess5=[80    900       80    900    ];


sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];

% a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1)];
% b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2)];
% c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3)];
% d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4)];
% e=[sGuess(1,5) sGuess(2,5) sGuess(3,5) sGuess(4,5) sGuess(5,5)];
% f=[sGuess(1,6) sGuess(2,6) sGuess(3,6) sGuess(4,6) sGuess(5,6)];
% g=[sGuess(1,7) sGuess(2,7) sGuess(3,7) sGuess(4,7) sGuess(5,7)];
% h=[sGuess(1,8) sGuess(2,8) sGuess(3,8) sGuess(4,8) sGuess(5,8)];
a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1)];
b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2)];
c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3)];
d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4)];

% for i=1:5
%     guessTemp(1,1)=a(i);
%     for j=1:5
%         guessTemp(1,2)=b(j);
%         for k=1:5
%             guessTemp(1,3)=c(k);
%             for l=1:5
%                 guessTemp(1,4)=d(l);
%                 for m=1:5
%                     guessTemp(1,5)=e(m);
%                     for n=1:5
%                         guessTemp(1,6)=f(n); 
%                           for o=1:5
%                               guessTemp(1,7)=g(o);
%                               for p=1:5
%                                   guessTemp(1,8)=h(p);
%                                   guess(p+(o-1)*2+(n-1)*4+(m-1)*8+(l-1)*16+(k-1)*32+64*(j-1)+128*(i-1),:)=guessTemp;
%                               end
%                           end
%                       end
%                   end
%               end
%           end
%       end
%   end
% end

for i=1:5
    guessTemp(1,1)=a(i);
    for j=1:5
        guessTemp(1,2)=b(j);
        for k=1:5
            guessTemp(1,3)=c(k);
            for l=1:5
                guessTemp(1,4)=d(l);
                guess(l+(k-1)*2+4*(j-1)+8*(i-1),:)=guessTemp;

            end
        end
    end
end



fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
fclose(fid);

for i=[1:11 13:14 16:60]
    x0=guess(i,:);
    opts=optimset('Display','iter','MaxFunEvals',200000,'MaxIter',10000,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
    [xRes(i,:),objVal(i,:),exitFlag(i,:),outPut(i)]=fminimax(@brDirectOptim,[x0 ones(1,30)],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1(indexSS)))],[10000 10000   10000 10000  5*ones(1,length(velData1(indexSS)))],@brDirectConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,v0,'SS','multi','fixed');
    %fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
    fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
    fprintf(fid,'\n x0: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',x0);
	fprintf(fid,'\n xRes: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',xRes(i,:));
    fprintf(fid,'\n objVal: %12.8f ',objVal(i));
%     fprintf(fid,'\n Output: ');
%     fprintf(fid,'\n iteration : %d ',outPut(i).iterations);
%     fprintf(fid,'\n funcCount : %d ',outPut(i).funcCount);
%     fprintf(fid,'\n stepsize : %d ',outPut(i).stepsize); 
	fprintf(fid,'\n--------------------------------------------------end of %d ----------------------------------------------------',i);
	fclose(fid);
end

fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
fclose(fid);

fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
fclose(fid);
% 
% for i=[11:61]
%     x0=guess(i,:);
%     opts=optimset('Display','iter','MaxFunEvals',20000,'MaxIter',1000,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
%     [xRes(i,:),objVal(i),exitFlag(i,:),outPut(i)]=fmincon(@brDirectOptim,[x0 ones(1,30)],[],[],[],[],[1e-4 1e-4  1e-4 1e-4  1e-4*ones(1,length(velData1(indexSS)))],[10000 10000   10000 10000  5*ones(1,length(velData1(indexSS)))],@brDirectConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,v0,'SS','single','fixed');
%     %fmincon(@brGeneralOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 ],[1000 100 100 3000 1000 100 100 3000],@brGeneralOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');
%     fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
%     fprintf(fid,'\n x0: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',x0);
% 	fprintf(fid,'\n xRes: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',xRes(i,:));
%     fprintf(fid,'\n objVal: %12.8f ',objVal(i));
%     fprintf(fid,'\n Output: ');
%     fprintf(fid,'\n iteration : %d ',outPut(i).iterations);
%     fprintf(fid,'\n funcCount : %d ',outPut(i).funcCount);
%     fprintf(fid,'\n stepsize : %d ',outPut(i).stepsize); 
% 	fprintf(fid,'\n--------------------------------------------------end of %d ----------------------------------------------------',i);
% 	fclose(fid);
% end
% 
% fid=fopen('H:\MatlabSave\guessSeekShrinkMulti.txt','a');
% fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
% fclose(fid);
% 
% save('resultats2.mat','xRes','objVal');