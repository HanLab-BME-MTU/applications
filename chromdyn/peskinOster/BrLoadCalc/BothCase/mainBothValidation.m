%main

delta=8e-3;
load('DataBoth');


pFixe='fixed';
guess=[10 1000 100 0.1];

% fid=fopen('H:\matlab\brModelAnal\brLoadCalc\BothCase\guessSeek.txt','a');
% fprintf(fid,'\n------------------------------------------ ----START-------------------------------------------------');
% fclose(fid);

for i=1:1
    x0=guess(i,:);
    opts=optimset('Display','iter','MaxFunEvals',2000000,'MaxIter',100,'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10);
    [xRes,objVal,exitFlag,outPut]=fmincon(@brBothDiffOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4 ],[10000 10000 10000 10000],@brBothDiffOptimConst,opts,velDataSG1,velDataSG2,cAngle1,cAngle2,delta,vUnload,pFixe)
%     fprintf(fid,'\n x0: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',x0);
% 	fprintf(fid,'\n xRes: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',xRes(i,:));
%     fprintf(fid,'\n objVal: %12.8f ',objVal(i));
%     fprintf(fid,'\n Output: ');
%     fprintf(fid,'\n iteration : %d ',outPut(i).iterations);
%     fprintf(fid,'\n funcCount : %d ',outPut(i).funcCount);
%     fprintf(fid,'\n stepsize : %d ',outPut(i).stepsize); 
% 	fprintf(fid,'\n--------------------------------------------------end of %d----------------------------------------------------',i);
% 	fclose(fid);
end

% fid=fopen('H:\matlab\brModelAnal\brLoadCalc\BothCase\guessSeek.txt','a');
% fprintf(fid,'\n------------------------------------------ ----END-------------------------------------------------');
% fclose(fid);
