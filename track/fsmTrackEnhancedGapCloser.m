function [M,closedGap]=fsmTrackEnhancedGapCloser(M,strg,threshold,userPath,firstIndex)
% fsmTrackMain uses the interpolated vector field to refine the gapcloser
%
% SYNOPSIS   [M,gapClosed]=fsmTrackEnhancedGapCloser(M,strg,threshold,userPath,firstIndex)
%
% INPUT      M          : M stack returned by the tracker
%            strg       : format string for the correct numbering of files
%            threshold  : radius of the region searched by the tracker for matching speckles
%            userPath   : work directory (where the interpolated vector field is saved)
%            firstIndex : index of the first corresponding image
%
% OUTPUT     M          : modified M stack
%            closedGap  : total number of closed gap in the M stack  
%
% DEPENDENCES   fsmTrackEnhancedGapCloser uses { fsmTrackPropSpecklePos ; fsmTrackGapCloser }
%               fsmTrackEnhancedGapCloser is used by { fsmTrackMain }
%
% Aaron Ponti, January 14th, 2003

if nargin~=5
    error('Wrong number of input arguments');
end

for c1=1:size(M,3)-1
    
    % Current index
    currentIndex=c1+firstIndex-1;

    % Initialize counter
    closedGap=0;
    
    % Load saved (interpolated) vector field
    indxStr=sprintf(strg,currentIndex); % For propagation of position at t1 (F(t1->t2))
    eval(['load ',userPath,filesep,'vectors',filesep,'vectors',indxStr,'.mat;']); % gapList
    vectors1=vectors;
    indxStr=sprintf(strg,currentIndex+1); % For backpropagation of position at t3 (-F(t2->t3))
    eval(['load ',userPath,filesep,'vectors',filesep,'vectors',indxStr,'.mat;']); % gapList
    vectors2=vectors; % Pass 'vectors' as is, it will be rotated 180 degrees when propagating
    % backwards in fsmTrackPropSpecklePos
    
    % Extract source and target positions to be (back)propagated
    source=M(:,1:2,c1);
    target=M(:,3:4,c1+1);
    
    % Keep track of the row in source and target which contain the speckles
    fS=find(source(:,1)~=0);
    fT=find(target(:,1)~=0);
    
    % Forward propagate
    pSource=fsmTrackPropSpecklePos(source(fS,:),vectors1,'forward');
    % Backward propagate
    pTarget=fsmTrackPropSpecklePos(target(fT,:),vectors2,'backward');
    
    % Create the copy of M (only three time-points) to pass to the gap closer
    Mcopy=M(:,:,c1:c1+1);
    Mcopy(fS,1:2,1)=pSource;
    Mcopy(fT,3:4,2)=pTarget;
    
    % Pass Mcopy to the gap closer
    [pMcopy,gapClosed]=fsmTrackGapCloser(Mcopy,threshold,strg,userPath,firstIndex,currentIndex+1);
    closedGap=closedGap+gapClosed;
    
    % Replace propagated positions with original positions and write back
    % to M (with gap closed)
    M(:,:,c1:c1+1)=pMcopy;
    M(:,1:2,c1)=source;
    M(:,3:4,c1+1)=target;
    
end

% Save also gapList for the first and last timepoints
gapList=[0 0];
indxStr=sprintf(strg,firstIndex);
gapListFirstFileName=[userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat'];
if ~exist(gapListFirstFileName); % Do not overwrite if it exists from a previous experiment
    eval(['save ',gapListFirstFileName,' gapList;']); 
end
    
indxStr=sprintf(strg,currentIndex+2);
gapListLastFileName=[userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat'];
if ~exist(gapListLastFileName); % Do not overwrite if it exists from a previous experiment
    eval(['save ',gapListLastFileName,' gapList;']); 
end
