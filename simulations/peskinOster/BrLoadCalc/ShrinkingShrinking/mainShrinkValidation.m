% Main: growing optimization validation
clear;
delta= 8e-3;

load('DataShrink');

vUnload=4.2;
pFixe='fixed';

% x0=[100 1500 0.1 100 1500 0.1];
% opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',1500,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
% [xRes,objVal,exitFlag,outPut]=myfmincon(@brBothShrinkOptim,x0,[],[],[],[],[10 1000 0 10 100 0],[ 500 3000 1 500 3000 1],@brBothShrinkOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);

%starting guess
% 
% sGuess1=[10 1000 0.01 10 1000 0.01];
% sGuess2=[20 1200 0.1 20 1200 0.1];
% sGuess3=[30 1600 0.3 30 1600 0.3];
% sGuess4=[60 2100 0.5 60 2100 0.5];
% sGuess5=[80 2500 0.7 80 2500 0.7];
% 
% 
% 
% sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];
% 
% 
% a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1) ];
% b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2) ];
% c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3) ];
% d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4) ];
% e=[sGuess(1,5) sGuess(2,5) sGuess(3,5) sGuess(4,5) sGuess(5,5) ];
% f=[sGuess(1,6) sGuess(2,6) sGuess(3,6) sGuess(4,6) sGuess(5,6) ];
% 
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
%                         guess((i-1)*32+(j-1)*16+(k-1)*8+(l-1)*4+(m-1)*2+n,:)=guessTemp;
%                     end
%                 end
%             end
%         end
%     end
% end


sGuess1=[10 1000  10 1000];
sGuess2=[20 1200  20 1200 ];
sGuess3=[30 1600  30 1600 ];
sGuess4=[60 2100  60 2100 ];
sGuess5=[80 2500  80 2500 ];
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



fid=fopen('H:\matlab\brModelAnal\brLoadCalc\Shrink\guessSeekShrink2.txt','a');
fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
fclose(fid);
for i=[1:200]
    x0=guess(i,:);
    opts=optimset('Display','iter','MaxFunEvals',2000,'MaxIter',250,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
    [xRes(i,:),objVal(i),exitFlag(i),outPut(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[10 1000  10 1000 ],[ 500 3000 500 3000 ],@brBothShrinkOptimConst,opts,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);
    fid=fopen('H:\matlab\brModelAnal\brLoadCalc\Shrink\guessSeekShrink2.txt','a');
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

fid=fopen('H:\matlab\brModelAnal\brLoadCalc\Shrink\guessSeekShrink2.txt','a');
fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
fclose(fid);

save('resultats2.mat','xRes','objVal');