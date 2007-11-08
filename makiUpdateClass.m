function dataStruct = makiUpdateClass(dataStruct)
%MAKIUPDATECLASS updates frame and kinetochore classification using tracking and sister information
%
% SYNOPSIS: dataStruct = makiUpdateClass(dataStruct)
%
% INPUT dataStruct (opt): maki data structure. If empty, it will be loaded
%                         via GUI
%
% OUTPUT dataStruct.updatedClass: structure of length nTimepoints with the
%                                 fields:
%               .inlierIdx       index into initCoord of all spots that are
%                                considered to be inliers
%               .unalignedIdx    index into initCoord of all spots that are
%                                considered to belong to lagging chromosomes 
%                                (occurs late prometaphase and anaphase)
%               .laggingIdx      index into initCoord of all spots that are
%                                considered to belong to lagging chromosomes
%                                (occurs only in anaphase frames)
%               .phase           'e' = prophase/early prometaphase (no plane)
%                                'p' = late prometaphase (plane, unaligned Ks)
%                                'm' = metaphase (plane, no unaligned Ks)
%                                'a' = anaphase (EITHER plane, 2 groups of Ks,
%                                some can be unaligned or lagging, OR no
%                                plane but after 'm' in end of movie)
%
% REMARKS The code assumes that the tracks do not contain merges and splits
%
% created by: kjaqaman
% DATE: 01-Nov-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input

%load dataStruct if not input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

%get number of frames in movie
numFrames = dataStruct.dataProperties.movieSize(end);

%get number of spots in each frame
numSpots = vertcat(dataStruct.initCoord.nSpots);

%read information from dataStruct
planeFit = dataStruct.planeFit;
tracks = dataStruct.tracks;
sisterList = dataStruct.sisterList;

%put all tracks in one big matrix
[dummy,tracksIndxMatrix] = convStruct2MatNoMS(tracks);

%% output

updatedClass = repmat(struct('inlierIdx',[],'unalignedIdx',[],'laggingIdx',[],...
    'phase',[],'tracksUnaligned',[],'tracksLagging',[],'sistersUnaligned',...
    [],'sistersLagging',[]),numFrames,1);
for iFrame = 1 : numFrames
    updatedClass(iFrame).unalignedIdx = planeFit(iFrame).unalignedIdx;
    updatedClass(iFrame).laggingIdx = planeFit(iFrame).laggingIdx;
end

%% outlier expansion based on sister pairs

if ~isempty(sisterList) && ~isempty(sisterList(1).trackPairs)

    %get number of sisters
    numSisters = length(sisterList);

    %go over all sisters ...
    for iSister = 1 : numSisters
       
        %get tracks making this sister pair
        trackIndx12 = sisterList(1).trackPairs(iSister,1:2);
        
        %get feature indices in these tracks
        featIndx1 = tracksIndxMatrix(trackIndx12(1),:)';
        featIndx2 = tracksIndxMatrix(trackIndx12(2),:)';
        
        %find frames where sister pair exists
        %place zeros where sister pair does not exist
        framesExist = find(~isnan(sisterList(iSister).distances(:,1)))';
        
        %go over frames and look for unaligned and lagging kinetochores
        %if one of the sisters is unaligned or lagging, add its sister to
        %the list
        for iFrame = framesExist
            
            %get unaligned and lagging kinetochores in this frame
            unalignedIdx = updatedClass(iFrame).unalignedIdx;
            laggingIdx = updatedClass(iFrame).laggingIdx;
        
            %check for unaligned kinetochores
            if any(unalignedIdx==featIndx1(iFrame)) || ...
                    any(unalignedIdx==featIndx2(iFrame))
                unalignedIdx = unique([unalignedIdx featIndx1(iFrame) ...
                    featIndx2(iFrame)]);
                updatedClass(iFrame).unalignedIdx = unalignedIdx;
            end
            
            %check for lagging kinetochores
            if any(laggingIdx==featIndx1(iFrame)) || ...
                    any(laggingIdx==featIndx2(iFrame))
                laggingIdx = unique([laggingIdx featIndx1(iFrame) ...
                    featIndx2(iFrame)]);
                updatedClass(iFrame).laggingIdx = laggingIdx;
            end
            
        end %(for iFrame = framesExist)
            
    end %(for iSister = 1 : numSisters)
    
end %(if ~isempty(sisterList))

%% outlier examination based on tracking information

tracksMatrixUnaligned = tracksIndxMatrix;
tracksMatrixLagging = tracksIndxMatrix;

%go over all frames and give a minus sign to outlier features
for iFrame = 1 : numFrames

    %get the unaligned and lagging kinetochores
    unalignedIdx = updatedClass(iFrame).unalignedIdx;
    laggingIdx = updatedClass(iFrame).laggingIdx;

    %find their locations in the matrix of tracks
    [dummy,unalignedLoc] = intersect(tracksIndxMatrix(:,iFrame),unalignedIdx);
    [dummy,laggingLoc] = intersect(tracksIndxMatrix(:,iFrame),laggingIdx);

    %put a negative sign in their location
    tracksMatrixUnaligned(unalignedLoc,iFrame) = ...
        -tracksMatrixUnaligned(unalignedLoc,iFrame);
    tracksMatrixLagging(laggingLoc,iFrame) = ...
        -tracksMatrixLagging(laggingLoc,iFrame);
    
end

%find tracks with lifetime > 4 frames
trackLifetime = getTrackSEL(tracks);
trackLifetime = trackLifetime(:,3);
goodTracks = find(trackLifetime > 4)';

%go over these tracks and see whether outlier features have other outlier 
%features in the same track
for iTrack = goodTracks
    
    %get number of outlier features in this track
    numOutliers = length(find(tracksMatrixUnaligned(iTrack,:)<0)) + ...
        length(find(tracksMatrixLagging(iTrack,:)<0));
    
    %change outliers back into inliers if number of outliers < 3
    %change is made by convering indices back into positive integers
    if numOutliers < 3
        tracksMatrixUnaligned(iTrack,:) = abs(tracksMatrixUnaligned(iTrack,:));
        tracksMatrixLagging(iTrack,:) = abs(tracksMatrixLagging(iTrack,:));
    end
    
end

%store unaligned and lagging kinetochore indices per frame
%store inlier indices and frame phase
for iFrame = 1 : numFrames
    unalignedIdx = abs(tracksMatrixUnaligned(tracksMatrixUnaligned(:,iFrame)<0,iFrame)');
    laggingIdx = abs(tracksMatrixLagging(tracksMatrixLagging(:,iFrame)<0,iFrame)');
    inlierIdx = setdiff((1:numSpots(iFrame)),[unalignedIdx laggingIdx]);
    framePhase = planeFit(iFrame).phase;
    if ~strcmp(framePhase,'a') && ~strcmp(framePhase,'e')
        if isempty(unalignedIdx)
            framePhase = 'm';
        else
            framePhase = 'p';
        end
    end
    updatedClass(iFrame).inlierIdx = inlierIdx;
    updatedClass(iFrame).unalignedIdx = unalignedIdx;
    updatedClass(iFrame).laggingIdx = laggingIdx;
    updatedClass(iFrame).phase = framePhase;
end

%% identification of tracks and sister pairs with outliers

%find tracks with unaligned and lagging kinetochores
tracksUnaligned = find(any(tracksMatrixUnaligned<0,2));
tracksLagging = find(any(tracksMatrixLagging<0,2));
updatedClass(1).tracksUnaligned = tracksUnaligned';
updatedClass(1).tracksLagging = tracksLagging';

if ~isempty(sisterList) && ~isempty(sisterList(1).trackPairs)

    %list the track indices making up sister pairs
    %look only at first track of each pair since either both tracks have
    %outliers or neither one has an outlier
    sisterTracks = sisterList(1).trackPairs(:,1);

    %find sisters with unaligned and laggin kinetochores
    [dummy,sistersUnaligned] = intersect(sisterTracks,tracksUnaligned);
    [dummy,sistersLagging] = intersect(sisterTracks,tracksLagging);
    updatedClass(1).sistersUnaligned = sistersUnaligned';
    updatedClass(1).sistersLagging = sistersLagging';
    
end

%% frame phase re-assignment

%get frame phases
framePhase = vertcat(updatedClass.phase);

%convert letters to numbers
framePhaseNum = zeros(size(framePhase));
framePhaseNum(framePhase=='p') = 1;
framePhaseNum(framePhase=='m') = 2;
framePhaseNum(framePhase=='a') = 3;

%find last frame that has a nonzero phase to start reassignment from 
lastFrame = find(framePhaseNum~=0,1,'last');

%go backwards over all frames and assign correct phase
for t = lastFrame-1 : -1 : 1
    framePhaseNum(t) = min(framePhaseNum(t:lastFrame));
end

%convert back to letters
framePhase(framePhaseNum==0) = 'e';
framePhase(framePhaseNum==1) = 'p';
framePhase(framePhaseNum==2) = 'm';
framePhase(framePhaseNum==3) = 'a';

%put phases back in structure
for t = 1 : numFrames
    updatedClass(t).phase = framePhase(t);
end

%% final output to dataStruct

dataStruct.updatedClass = updatedClass;


%% ~~~ the end ~~~