function fct=brGeneralOptim(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);
% brBothShrinkOptim resuturn a obj fct for genral case
% INPUT
%     x         :  paramters of depol [pol depol pol depol]
%     velData1  :  velocity of the MT1 with sign
%     velData2  :  velocity of MT2 with sign
%     cAngle1   :  cosines of angle between MT1 and inter-centromere axis
%     cAngle2   :  cosines of angle between MT2 and inter-centromere axis
%     delta     :  1/13 of subunit length
%     vUnload   :  unload velocity of the shrink case
%     pFixe     :  fixed or free regarding if you optimize with p or
%     without
% OUPUT 
%     fct       :  objective function
% COMMENT : the free case is not always update

%DEBUG
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

fctGG=brBothGrowingOptimFct([paramG1 paramG2],velData1(indexGG),velData2(indexGG),cAngle1(indexGG),cAngle2(indexGG),delta/13);

%------------------------------- SS case

fctSS=brBothShrinkOptim([paramS1 paramS2] ,velData1(indexSS),velData2(indexSS),cAngle1(indexSS),cAngle2(indexSS),delta,vUnload,pFixe);

%------------------------------- GS case

fctGS=brBothDiffOptim([paramG1 paramS2],velData1(indexGS),velData2(indexGS),cAngle1(indexGS),cAngle2(indexGS),delta,vUnload,pFixe);
    
    
%------------------------------- SG case

fctSG=brBothDiffOptim([paramS1 paramG2],velData1(indexSG),velData2(indexSG),cAngle1(indexSG),cAngle2(indexSG),delta,vUnload,pFixe);
   
fct=sum(fctGG)+fctSS+fctGS+fctSG;

