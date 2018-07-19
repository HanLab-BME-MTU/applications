% Main: growing optimization validation
clear;
delta= 8e-3;
path= [];%'H:\matlabSave\'

%[velData1,velData2,velDataNoise51,velDataNoise52,velDataNoise451,velDataNoise452,cAngle1,cAngle2,omega1,omega2,v0,xObj]=DataGenerationShrink(30);% DataGeneration fro a validation
load('DataShrink')

vUnload=4.2;
pFixe='fixed';

v0=0.08;

x0=[100 2000 200 1000];
opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',1500,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
%[xRes,objVal,exitFlag,outPut]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[10 1000  10 100 ],[ 500 3000  500 3000 ],@brBothShrinkOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,'fixed');

% starting guess

sGuess1=[0.1 50  0.1 50];
sGuess2=[0.5 75 0.5 75];
sGuess3=[1  100 1  100];
sGuess4=[1.5 125 1.5 125];
sGuess5=[2   150 2   150];


sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];

a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1) ];
b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2) ];
c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3) ];
d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4) ];


for i=1:5
    guessTemp(1,1)=a(i);
    for j=1:5
        guessTemp(1,2)=b(j);
        for k=1:5
            guessTemp(1,3)=c(k);
            for l=1:5
                guessTemp(1,4)=d(l);   
                guess(l+(k-1)*2+(j-1)*4+(i-1)*8,:)=guessTemp;
              
            end
        end
    end
end



fid=fopen('guessSeekShrink2.txt','a');
fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
fclose(fid);
%1:30
for i=[1:30 32:61]
    x0=guess(i,:);
    
    opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
    [xRes(i,:),objVal(i),exitFlag(i),outPut(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-5 1e-5 1e-5 1e-5],[ 500 3000 500 3000 ],@brBothShrinkOptimConst,opts,velDataNoise451,velDataNoise452,cAngle1,cAngle2,delta,v0,pFixe);
    fid=fopen('guessSeekShrink2.txt','a');
    fprintf(fid,'\n x0: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',x0);
	fprintf(fid,'\n xRes: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f ',xRes(i,:));
    fprintf(fid,'\n objVal: %12.8f ',objVal(i));
    fprintf(fid,'\n Output: ');
    fprintf(fid,'\n iteration : %d ',outPut(i).iterations);
    fprintf(fid,'\n funcCount : %d ',outPut(i).funcCount);
    fprintf(fid,'\n stepsize : %d ',outPut(i).stepsize); 
	fprintf(fid,'\n--------------------------------------------------end of %d ----------------------------------------------------',i);
	fclose(fid);
end

fid=fopen('guessSeekShrink2.txt','a');
fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
fclose(fid);

save('resultats2Noise45_2.mat','xRes','objVal','exitFlag','omega1','omega2','velData1','velData2','v0','velDataNoise451','velDataNoise452','velDataNoise51','velDataNoise51','xObj');