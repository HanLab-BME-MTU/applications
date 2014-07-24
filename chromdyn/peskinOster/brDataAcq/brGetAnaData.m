function [mt1State,mt2State,mt1Velocity,mt2Velocity,cAngle1,cAngle2,trajMatMt1,trajMatMt2,dataSpb]=brGetAnaData(idlisttrack,dataProperties,lastResult)
%BRGETANADAT put the data in the needed format for the optimzation for 4
%tags

%SYNOPSIS
%[mt1State,mt2State,mt1Velocity,mt2Velocity,cAngle1,cAngle2,trajMatMt1,trajMatMt2,dataSpb]=brGetAnaData(idlisttrack,dataProperties,lastResult)
%    
%
%INPUT 
%       idlisttrack     : list from jonas tracking program
%       dataProperties  : properties from jonas tracking program
%       lastResult      : last result from jonas tracking program
%
%OUTPUT
%       mt1Velocity     : vector with the velocity for first MT
%       mt2Velocity     : vector with the velocity for second MT
%       cAngle1         : vector with the cosinus of the angle between the
%                         first MT and the inter-centromere axis
%       cAngle2         : vector with the cosinus of the angle between the
%                         second MT and the inter-centromere axis
%       tarjMatMt1      : matrix returning all the info about the
%                         of the first MT see trajectoryAnalysis.m for more information.
%       tarjMatMt2      : matrix returning all the info about the trajectory
%                         of the second MT see trajectoryAnalysis.m for more information.
%       dataSpb         : matrix returning all the info about the trajectory
%                         of the spindle pole body axis see trajectoryAnalysis.m for more information.
%COMMENTS: 




%-------- tag-name link---------------

spb1Tag=[];
spb2Tag=[];
cen1Tag=[];
cen2Tag=[];


spb1Tag=find(strcmp('spb1',idlisttrack(1).stats.labelcolor)==1);
spb2Tag=find(strcmp('spb2',idlisttrack(1).stats.labelcolor)==1);
cen1Tag=find(strcmp('cen1',idlisttrack(1).stats.labelcolor)==1);
cen2Tag=find(strcmp('cen2',idlisttrack(1).stats.labelcolor)==1);


dataMt1 = calculateTrajectoryFromIdlist(idlisttrack,dataProperties,cen1Tag,spb1Tag,[]);
ioOpt.verbose = 1;
ioOpt.saveTxt=0;
ioOpt.saveMat=0;
ioOpt.clusterData=0;
trajectoryDescriptionMt1 = trajectoryAnalysis(dataMt1,ioOpt,[]);

dataMt2 = calculateTrajectoryFromIdlist(idlisttrack,dataProperties,cen2Tag,spb2Tag,[]);

dataSpb = calculateTrajectoryFromIdlist(idlisttrack,dataProperties,spb1Tag,spb2Tag,[]);


trajectoryDescriptionMt2 = trajectoryAnalysis(dataMt2,ioOpt,[]);
close all;

% put the trajectory at the same dimension: some vel are calculated over
% more than two pts -> take an average val.
frameMatMt1=cat(1,trajectoryDescriptionMt1.individualStatistics.dataListGroup);
frameMatMt1=frameMatMt1(:,[1:5 7]);

frameMatMt2=cat(1,trajectoryDescriptionMt2.individualStatistics.dataListGroup);
frameMatMt2=frameMatMt2(:,[1:5 7]);




[trajMatMt1]=brTrajectoryModif(trajectoryDescriptionMt1);
[trajMatMt2]=brTrajectoryModif(trajectoryDescriptionMt2);

%  removing frames no present in both mts
%  (the non growing & non shrinking have already been
%  discarded in  brTrajectoryModif
indiceOKMt1=[];
indiceOKMt2=[];


if size(trajMatMt1)<size(trajMatMt2)
	iLength=size(trajMatMt1);
    for i=1:iLength
        if ~isempty(find(trajMatMt1(i,1)==trajMatMt2(:,1)))
		    indiceOKMt1=[indiceOKMt1 i];
            indiceOKMt2=[indiceOKMt2 find(trajMatMt1(i,1)==trajMatMt2(:,1))];
        end
	end
else
    iLength=size(trajMatMt2);
        for i=1:iLength
        if ~isempty(find(trajMatMt2(i,1)==trajMatMt1(:,1)))
		    indiceOKMt2=[indiceOKMt2 i];
            indiceOKMt1=[indiceOKMt1 find(trajMatMt2(i,1)==trajMatMt1(:,1))];
        end
	end
end




trajMatMt1=trajMatMt1(indiceOKMt1,:);
trajMatMt2=trajMatMt2(indiceOKMt2,:);

if (isempty(spb1Tag)|isempty(spb2Tag)|isempty(cen1Tag)|isempty(cen2Tag))
    disp('All the four tags have not be indentified');
end

%--------  anaDat acq ------------------------------
anaDat = adgui_calc_anaDat(idlisttrack,dataProperties,lastResult);

%------ Calculation of the change of length /dt-----
for i=1:length(anaDat)
    mt1Length(i)=anaDat(i).distanceMatrix(spb1Tag,cen1Tag);
    mt2Length(i)=anaDat(i).distanceMatrix(spb2Tag,cen2Tag);
end
% elimination of the zero
% mt1Length=mt1Length(find(mt1Length));
% mt2Length=mt2Length(find(mt2Length));



  
  
mt1Velocity=diff(mt1Length')./diff(cat(1,anaDat.time));
mt2Velocity=diff(mt2Length')./diff(cat(1,anaDat.time));

mt1State(find(mt1Velocity>0))=1;
mt1State(find(mt1Velocity<0))=-1;
mt1State(find(mt1Velocity==0))=0;

mt2State(find(mt2Velocity>0))=1;
mt2State(find(mt2Velocity<0))=-1;
mt2State(find(mt2Velocity==0))=0;

%------- Calculation of the mt direction------------
%-------  cent direction ---------------------------

for i=1:length(anaDat)
    mt1Coord(i,:)=anaDat(i).coord(cen1Tag,:)-anaDat(i).coord(spb1Tag,:);
    mt2Coord(i,:)=anaDat(i).coord(cen2Tag,:)-anaDat(i).coord(spb2Tag,:);
    cenCoord(i,:)=anaDat(i).coord(cen2Tag,:)-anaDat(i).coord(cen1Tag,:);
end





%------Calculation vel vect (rate of change on unity vect MT)---
% the change of length is scalar and the vect is directed to cen

 % the change of length is scalar and the vect is directed to cen
 
% abs of the scalar product. abs -> get the angle between 0 and pi/2
dotVect1=(sum((mt1Coord.*cenCoord)')');
normCen=sum((cenCoord.*cenCoord)').^0.5;
normMt1=sum((mt1Coord.*mt1Coord)').^0.5;
normVects1=(normCen.*normMt1)';

normVects1(find(normVects1==0))=NaN;

cAngle1=dotVect1./normVects1;


dotVect2=sum((mt2Coord.*(-cenCoord))')';
normCen=sum((cenCoord.*cenCoord)').^0.5;
normMt2=sum((mt2Coord.*mt2Coord)').^0.5;
normVects2=(normCen.*normMt2)';
normVects2(find(normVects2==0))=NaN;
cAngle2=dotVect2./normVects2;






