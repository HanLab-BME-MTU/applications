function [c,ceq]=brGeneralOptimConst(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);



if isempty(find(x<0))==0
    x=[0.01 0.001 0.001 0.001 0.01 0.1 0.001 0.011 0.01 0.001];
end

% state separation

indexGG=find(velData1>0 & velData2>0);
indexSS=find(velData1<0 & velData2<0);
indexGS=find(velData1>0 & velData2<0);
indexSG=find(velData1<0 & velData2>0);

switch pFixe
    case 'fixed'
		paramG1=x(1:2);
		paramS1=x(3:4);
		paramG2=x(5:6);
		paramS2=x(7:8);
    case 'free'
		paramG1=x(1:2);
		paramS1=x(3:5);
		paramG2=x(6:7);
		paramS2=x(8:10);
end
%------------------------------- GG case

cGG=brBothGrowingOptimFctConst([paramG1 paramG2],velData1(indexGG),velData2(indexGG),cAngle1(indexGG),cAngle2(indexGG),delta);

%------------------------------- SS case

cSS=brBothShrinkOptimConst([paramS1 paramS2] ,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,pFixe);

%------------------------------- GS case

cGS=brBothDiffOptimConst([paramG1 paramS2],velData1(indexGS),velData2(indexGS),cAngle1(indexGS),cAngle2(indexGS),delta,vUnload,pFixe);
    


%------------------------------- SG case

cSG=brBothDiffOptimConst([paramS1 paramG2],velData1(indexSG),velData2(indexSG),cAngle1(indexSG),cAngle2(indexSG),delta,vUnload,pFixe);


% omegaDiscret=[0:0.001:5];
% 
% velDiscretS1Max=max(brShrinkDirect(paramS1,omegaDiscret,delta));
% velDiscretS2Max=max(brShrinkDirect(paramS2,omegaDiscret,delta));
% 
% velDiscretG1Max=delta*(x(1)-x(2));
% velDiscretG2Max=delta*(x(6)-x(7));
% 
% 
% velGMax1=max([velData1(indexGG) velData1(indexGS)]);
% velGMax2=max([velData2(indexGG) velData2(indexSG)]);
% 
% velSMax1=max([velData1(indexSS) velData1(indexSG)]);
% velSMax2=max([velData2(indexSS) velData2(indexGS)]);



c=[cGG' cSS cGS cSG];

v0=vUnload;
c(end+1)=(x(3)-(delta*x(4)*x(3))/v0)/x(4);
c(end+1)=(x(3)-(delta*x(4)*x(3))/v0)/x(4);
c(end+1)=-(x(3)-(delta*x(4)*x(3))/v0)/x(4)-1;
c(end+1)=-(x(3)-(delta*x(4)*x(3))/v0)/x(4)-1;
% c(end+1)=0.95*velDiscretS1Max-velSMax1;
% c(end+1)=velGMax1-1.3*velDiscretG1Max;
% c(end+1)=0.7*velDiscretG1Max-velGMax1;
% c(end+1)=velGMax2-1.3*velDiscretG2Max;
% c(end+1)=0.7*velDiscretG2Max-velGMax2;

% c(end+1:end+5)=(x(1:5)-x(6:10))./(0.5*x(1:5)+x(6:10))-10;


c=c';
ceq=[];
