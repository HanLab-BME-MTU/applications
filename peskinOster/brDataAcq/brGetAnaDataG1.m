function [trajectoryDescription,data,mt1Velocity,frameNb,mt1VelCons,frameCons]=brGetAnaDataG1(idlisttrack,dataProperties,lastResult)
%BRGETANADAT put the data in the needed format for the optimzation for 2
%tag

%SYNOPSIS
%[trajectoryDescription,data,mt1Velocity,frameNb,mt1VelCons,frameCons]=brGetAnaDataG1(idlisttrack,dataProperties,lastResult)
%    
%
%INPUT 
%       idlisttrack     : list from jonas tracking program
%       dataProperties  : properties from jonas tracking program
%       lastResult      : last result from jonas tracking program
%
%OUTPUT
%       trajectoryDescription      : matrix returning all the info about the
%                                    of the first MT see trajectoryAnalysis.m for more information.
%       data                       : matrix returning all the info about
%                                    the trajectory from
%                                    calculateTrajectoryFromIdlist.m
%       mt1Velocity                : vector with the velocities
%       frameNb                    : number of frame analyzed
%       mt1VelCons                 : vector of consequtive velocities
%       frame Cons                 : vector of consequtive frames

%-------- tag-name link---------------

spb1Tag=[];
cen1Tag=[];



spb1Tag=find(strcmp('spb1',idlisttrack(1).stats.labelcolor)==1);
cen1Tag=find(strcmp('cen1',idlisttrack(1).stats.labelcolor)==1);

data = calculateTrajectoryFromIdlist(idlisttrack,dataProperties,cen1Tag,spb1Tag,[]);
ioOpt.verbose = 1;
ioOpt.saveTxt=0;
ioOpt.saveMat=0;
ioOpt.clusterData=0;
trajectoryDescription = trajectoryAnalysis(data,ioOpt,[]);

if (isempty(spb1Tag)|isempty(cen1Tag))
    disp('All two tags have not be indentified');
end

%--------  anaDat acq ------------------------------
anaDat = adgui_calc_anaDat(idlisttrack,dataProperties,lastResult);

%------ Calculation of the change of length /dt-----
for i=1:length(anaDat)
    mt1Length(i)=anaDat(i).distanceMatrix(spb1Tag,cen1Tag);
end
  
mt1Velocity=diff(mt1Length')./diff(cat(1,anaDat.time));


mt1State(find(mt1Velocity>0))=1;
mt1State(find(mt1Velocity<0))=-1;
mt1State(find(mt1Velocity==0))=0;

%----- elimiatin of the vel calculated over more than 2 frames
frameNb=cat(1,anaDat.timePoint);
frameConseqIndex=find(diff(frameNb)==1);
frameCons=frameNb(frameConseqIndex);
for i=1:length(frameConseqIndex)
    mt1DeltaLength(i)=anaDat(frameConseqIndex(i)+1).distanceMatrix(spb1Tag,cen1Tag)-anaDat(frameConseqIndex(i)).distanceMatrix(spb1Tag,cen1Tag);
end
mt1VelCons=mt1DeltaLength'./(cat(1,anaDat(frameConseqIndex+1).time)-cat(1,(anaDat(frameConseqIndex).time)));

