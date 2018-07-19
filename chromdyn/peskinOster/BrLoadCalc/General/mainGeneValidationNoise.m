% Starting guess evalution

% main Test
clear;
close all;

delta=8e-3;

load('h:\MatlabSave\velDataGeneralGeneral');
load('h:\MatlabSave\velDataGeneralGeneral2');
%v0=vUnload;
velData1o=velData1;
velData2o=velData2;
velData21o=velData21;
velData22o=velData22;

indexGG=find(velData1>0 &velData2>0);
indexSS=find(velData1<0 &velData2<0);

index2GG=find(velData21>0 &velData22>0);
index2SS=find(velData21<0 &velData22<0);

sGuess1=[1     1        1     1     ];
sGuess2=[20   100      20     100     ];
sGuess3=[40    300      40    300   ];
sGuess4=[60    600        60    600     ];
sGuess5=[80    900       80    900    ];


sGuess=[sGuess1; sGuess2;sGuess3; sGuess4;sGuess5];


a=[sGuess(1,1) sGuess(2,1) sGuess(3,1) sGuess(4,1) sGuess(5,1)];
b=[sGuess(1,2) sGuess(2,2) sGuess(3,2) sGuess(4,2) sGuess(5,2)];
c=[sGuess(1,3) sGuess(2,3) sGuess(3,3) sGuess(4,3) sGuess(5,3)];
d=[sGuess(1,4) sGuess(2,4) sGuess(3,4) sGuess(4,4) sGuess(5,4)];

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


[fct]=brBothGrowingOptimFct([420 0.2 382 0.9],velData1(indexGG),velData2(indexGG),cAngle1(indexGG),cAngle2(indexGG),delta/13);
opts=optimset('Display','iter','MaxFunEvals',3000,'MaxIter',2500,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);%,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);



%1->3 5-.11
for i=13:61
    
    x0=guess(i,:);
    i
    opts=optimset('Display','iter','MaxFunEvals',3000,'MaxIter',2500,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);
    load('H:\MatlabSave\velDataGeneralGeneralGOOD5.mat');
    [xRes5SS1(i,:),objVal5SS1(i),exit5Flag1(i),out5Put(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
    i
    [xRes5SS2(i,:),objVal5SS2(i),exit5Flag2(i),out5Put(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velDataNoise1(indexSS),velDataNoise2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
    i

    load('H:\MatlabSave\velDataGeneralGeneraGOOD45.mat')
    [xRes45SS1(i,:),objVal45SS1(i),exit45Flag1(i),out45Put(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
    i
    [xRes45SS2(i,:),objVal45SS2(i),exit45Flag2(i),out45Put(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velDataNoise1(indexSS),velDataNoise2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');

end


for i=1:61
    figure('Name',[num2str(i) 'Good, 5%']);hold on;
    a=['Obj: ' num2str(objVal5SS2(i))];
    b=['Obj: ' num2str(objVal5SS1(i))];
    c=(xRes5SS2(i,:)-[xObj(3:4) xObj(8:9)])./[xObj(3:4) xObj(8:9)];
    text(2,0.6,['noise: ' a]);
    text(2,0.3,['no noise: ' b]);
    text(1,0.1,['relative error: ' '[' num2str(c(1)) ' ; ' num2str(c(2)) ' ; ' num2str(c(3)) ' ; ' num2str(c(4)) ']']);
	omega=[0:0.01:5];
    pR1=(xRes5SS2(i,1)-(xRes5SS1(i,2)*delta*xRes5SS2(i,1))/v0)/xRes5SS2(i,2)+1;
    pR2=(xRes5SS2(i,3)-(xRes5SS2(i,4)*delta*xRes5SS2(i,3))/v0)/xRes5SS2(i,4)+1;
    pO1=xObj(5);
    pO2=xObj(10);
    
	plot(omega,brShrinkDirect([xRes5SS2(i,3:4) pR2],omega,delta),'k');
	plot(omega,brShrinkDirect([xObj(8:10)],omega,delta),'-.k');
	plot(omega2(indexSS),-velDataNoise2(indexSS),'dk');
	plot(omega,brShrinkDirect([xRes5SS2(i,1:2) pR1],omega,delta),'r');
	plot(omega,brShrinkDirect([xObj(3:5) ],omega,delta),'-.r');
	plot(omega1(indexSS),-velDataNoise1(indexSS),'dr');hold off;
        figure('Name',[num2str(i) 'Good, 45%']);hold on;
    a=['Obj: ' num2str(objVal45SS2(i))];
    b=['Obj: ' num2str(objVal45SS1(i))];
    c=(xRes45SS2(i,:)-[xObj(3:4) xObj(8:9)])./[xObj(3:4) xObj(8:9)];
    text(2,0.6,['noise: ' a]);
    text(2,0.3,['no noise: ' b]);
    text(1,0.1,['relative error: ' '[' num2str(c(1)) ' ; ' num2str(c(2)) ' ; ' num2str(c(3)) ' ; ' num2str(c(4)) ']']);
	omega=[0:0.01:5];
    pR1=(xRes45SS2(i,1)-(xRes45SS1(i,2)*delta*xRes45SS2(i,1))/v0)/xRes45SS2(i,2)+1;
    pR2=(xRes45SS2(i,3)-(xRes45SS2(i,4)*delta*xRes45SS2(i,3))/v0)/xRes45SS2(i,4)+1;
    pO1=xObj(5);
    pO2=xObj(10);
    
	plot(omega,brShrinkDirect([xRes45SS2(i,3:4) pR2],omega,delta),'k');
	plot(omega,brShrinkDirect([xObj(8:10)],omega,delta),'-.k');
	plot(omega2(indexSS),-velDataNoise2(indexSS),'dk');
	plot(omega,brShrinkDirect([xRes45SS2(i,1:2) pR1],omega,delta),'r');
	plot(omega,brShrinkDirect([xObj(3:5) ],omega,delta),'-.r');
	plot(omega1(indexSS),-velDataNoise1(indexSS),'dr');hold off;
end
    

% for i=1:25
%     DataGeneGeneral
%     x0=[20 200 10 200];
%     opts=optimset('Display','iter','MaxFunEvals',3000,'MaxIter',2500,'TolFun',1e-7,'TolX',1e-7,'TolCon',1e-5);,'DiffMaxChange',1e-2,'DiffMinChange',1e-10);
% 
%     [xResSS1(i,:),objValSS1(i),exitFlag1(i),outPut(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
%     [xResSS2(i,:),objValSS2(i),exitFlag2(i),outPut(i)]=fmincon(@brBothShrinkOptim,x0,[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[1000 1000 1000 1000 ],@brBothShrinkOptimConst,opts,velDataNoise1(indexSS),velDataNoise2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
%     
% 	figure('Name',num2str(i));hold on;
%     a=['Obj: ' num2str(objValSS2(i))];
%     b=['Obj: ' num2str(objValSS1(i))];
%     c=(xResSS2(i,:)-[xObj(3:4) xObj(8:9)])./[xObj(3:4) xObj(8:9)];
%     text(2,0.6,['noise: ' a]);
%     text(2,0.3,['no noise: ' b]);
%     text(1,0.1,['relative error: ' '[' num2str(c(1)) ' ; ' num2str(c(2)) ' ; ' num2str(c(3)) ' ; ' num2str(c(4)) ']']);
% 	omega=[0:0.01:5];
%     pR1=(xResSS2(i,1)-(xResSS1(i,2)*delta*xResSS2(i,1))/v0)/xResSS2(i,2)+1;
%     pR2=(xResSS2(i,3)-(xResSS2(i,4)*delta*xResSS2(i,3))/v0)/xResSS2(i,4)+1;
%     pO1=xObj(5);
%     pO2=xObj(10);
%     
% 	plot(omega,brShrinkDirect([xResSS2(i,3:4) pR2],omega,delta),'k');
% 	plot(omega,brShrinkDirect([xObj(8:10)],omega,delta),'-.k');
% 	plot(omega2(indexSS),-velDataNoise2(indexSS),'dk');
% 	plot(omega,brShrinkDirect([xResSS2(i,1:2) pR1],omega,delta),'r');
% 	plot(omega,brShrinkDirect([xObj(3:5) ],omega,delta),'-.r');
% 	plot(omega1(indexSS),-velDataNoise1(indexSS),'dr');hold off;
% end
% 


% [xResSS1,objValSS1,exitFlag1,outPut]=fmincon(@brBothShrinkOptim,[xObj(3:4) xObj(8:9)],[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[100 1000 100 1000 ],@brBothShrinkOptimConst,opts,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');
% [xResSS2,objValSS2,exitFlag2,outPut]=fmincon(@brBothShrinkOptim,[xObj(3:4) xObj(8:9)],[],[],[],[],[1e-4 1e-4 1e-4 1e-4  ],[100 1000 100 1000 ],@brBothShrinkOptimConst,opts,velDataNoise1(indexSS),velDataNoise2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,'fixed');

