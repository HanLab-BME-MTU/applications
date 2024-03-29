function sisterList = makiConstructAlignedSisters(dataStruct,removeNetDisp,...
    randomize,correctStd,chooseLeftRight,verbose)
%MAKICONSTRUCTALIGNEDSISTERS extracts aligned coordinates of sister pairs
%
%SYNOPSIS sisterList = makiConstructAlignedSisters(dataStruct,removeNetDisp,...
%    randomize)
%
%INPUT  dataStruct: dataStruct as generated by makiMovieAnalysis. Must have
%                   all tasks (1-9, except 4) performed.
%       removeNetDisp: 1 - Make center of each sister pair at zero
%                      throughout the whole movie, 0 otherwise.
%                      Optional. Default: 0.
%       randomize: 1 - Randomize sister pairing, 0 otherwise.
%                  Optional. Default: 0.
%       correctStd: 1 to correct stds (because they are underestimated in
%                   initCoord), 0 otherwise. Optional. Default: 0.
%       chooseLeftRight: 1 to make sister 1 always sister on the left, 0
%                   to keep things not ordered. Optional. Default: 0.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%       
%OUTPUT sisterList: Same as input (dataStruct.sisterList) but with
%                   additional fields (.coords1Aligned, .coords2Aligned 
%                   and .distanceAligned) that contain the aligned
%                   coordinates and distance, with whatever additional 
%                   modifications requested in the input.
%
%Khuloud Jaqaman, February 2008

%% input
if nargin < 2 || isempty(removeNetDisp)
    removeNetDisp = 0;
end

if nargin < 3 || isempty(randomize)
    randomize = 0;
end

if nargin < 4 || isempty(correctStd)
    correctStd = 0;
end

if nargin < 5 || isempty(chooseLeftRight)
    chooseLeftRight = 0;
end

if nargin < 6 || isempty(verbose)
    verbose = 0;
end

%copy fields out of dataStruct
sisterList = dataStruct.sisterList;
tracks = dataStruct.tracks;
frameAlignment = dataStruct.frameAlignment;
% initCoord = dataStruct.initCoord;

%get number of sisters
numSisters = length(sisterList);

%% get aligned coordinates

%go over all sisters
for iSister = 1 : numSisters

    %find track indices
    tracksIndx = sisterList(1).trackPairs(iSister,1:2);

    %determine frame where each track starts
    trackStart = [tracks(tracksIndx(1)).seqOfEvents(1,1) ...
        tracks(tracksIndx(2)).seqOfEvents(1,1)];

    %find number of frames and frames where pair "exists"
    goodFrames = ~isnan(sisterList(iSister).distances(:,1));
    numFrames = length(goodFrames);
    goodFrames = find(goodFrames);

    %find feature indices making up sisters
    sisterIndx1 = NaN(numFrames,1);
    sisterIndx2 = NaN(numFrames,1);
    sisterIndx1(goodFrames) = tracks(tracksIndx(1))...
        .tracksFeatIndxCG(goodFrames-trackStart(1)+1);
    sisterIndx2(goodFrames) = tracks(tracksIndx(2))...
        .tracksFeatIndxCG(goodFrames-trackStart(2)+1);

    %get aligned sister coordinates
    coords1 = NaN(numFrames,6);
    coords2 = NaN(numFrames,6);
    for iFrame = goodFrames'
        coords1(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx1(iFrame),:);
        coords2(iFrame,:) = frameAlignment(iFrame).alignedCoord(sisterIndx2(iFrame),:);
    end
    
    %correct position standard deviation if necessary
    if correctStd
        coords1(:,4:6) = coords1(:,4:6) * 1.5;
        coords2(:,4:6) = coords2(:,4:6) * 1.5;
    end
    
    %sort sisters to left and right
    if chooseLeftRight
        
        %calculate the average coordinate of each sister along the normal
        %to the metaphase plate
        meanCoord1Normal = nanmean(coords1(:,1));
        meanCoord2Normal = nanmean(coords2(:,1));
        
        %put the sister with the smaller average coordinate on the "left"
        %negative is smaller than positive, no matter the magnitude
        if meanCoord2Normal < meanCoord1Normal
            tmp = coords2;
            coords2 = coords1;
            coords1 = tmp;
        end
    
    end
    
    %store aligned coordinates in sisterList
    sisterList(iSister).coords1Aligned = coords1;
    sisterList(iSister).coords2Aligned = coords2;
    sisterList(iSister).distanceAligned = sqrt(sum((coords1(:,1:3) - coords2(:,1:3)).^2,2));
    
end

%% remove net displacement

if removeNetDisp

    %go over all sisters
    for iSister = 1 : numSisters

        %get coordinates
        coords1T = sisterList(iSister).coords1Aligned;
        coords2T = sisterList(iSister).coords2Aligned;

        %calculate center of mass
        centerMass = (coords1T(:,1:3) + coords2T(:,1:3)) / 2;

        %shift coordinates such that center of mass is always at 0
        sisterList(iSister).coords1Aligned(:,1:3) = coords1T(:,1:3) - centerMass;
        sisterList(iSister).coords2Aligned(:,1:3) = coords2T(:,1:3) - centerMass;

    end

end

%% randomize kinetochore pairing

if randomize

    %first, list pairs such that sister 1 is the "left" kinetochore and
    %sister 2 is the "right" kinetochore
    for iSister = 1 : numSisters
        
        %calculate the average coordinate of each sister along the normal
        %to the metaphase plate
        meanCoord1Normal = nanmean(sisterList(iSister).coords1Aligned(:,1));
        meanCoord2Normal = nanmean(sisterList(iSister).coords2Aligned(:,1));
        
        %put the sister with the smaller average coordinate on the "left"
        %negative is smaller than positive, no matter the magnitude
        if meanCoord2Normal < meanCoord1Normal
            tmp = sisterList(1).trackPairs(iSister,1);
            sisterList(1).trackPairs(iSister,1) = sisterList(1).trackPairs(iSister,2);
            sisterList(1).trackPairs(iSister,2) = tmp;
            tmp = sisterList(iSister).coords1Aligned;
            sisterList(iSister).coords1Aligned = sisterList(iSister).coords2Aligned;
            sisterList(iSister).coords2Aligned = tmp;
        end
        
    end

    %     %second, shuffle the "right" sisters to make random pairs
    %     indxRand = randsample(numSisters,numSisters);
    %
    %     %store original sister list in a new variable
    %     sisterListOriginal = sisterList;
    %
    %     %go over all sisters
    %     for iSister = 1 : numSisters
    %
    %         %assign new "right" sister
    %         sisterList(1).trackPairs(iSister,2) = sisterListOriginal(1).trackPairs(indxRand(iSister),2);
    %
    %         %move corresponding coordinates
    %         sisterList(iSister).coords2Aligned = sisterListOriginal(indxRand(...
    %             iSister)).coords2Aligned;
    %
    %         %calculate new distance vector and store
    %         sisterList(iSister).distanceAligned = sqrt(sum((sisterList(...
    %             iSister).coords1Aligned(:,1:3) - sisterList(iSister...
    %             ).coords2Aligned(:,1:3)).^2,2));
    %
    %     end

    %second, shuffle sisters such that "left" kinetochores are paired and
    %"right" kinetochores are paired
    indx1Rand = randsample(numSisters,numSisters);
    indx2Rand = randsample(numSisters,numSisters);

    %if there is an odd number of sisters, remove the last
    if mod(numSisters,2)
        indx1Rand = indx1Rand(1:end-1);
        indx2Rand = indx2Rand(1:end-1);
        numSisters = numSisters - 1;
    end

    %store original sister list in a new variable
    sisterListOriginal = sisterList;

    %go over all sisters
    for iSister = 1 : numSisters/2

        %assign new sister pair
        sisterList(1).trackPairs(iSister,1) = sisterListOriginal(...
            1).trackPairs(indx1Rand(2*iSister-1),1);
        sisterList(1).trackPairs(iSister,2) = sisterListOriginal(...
            1).trackPairs(indx1Rand(2*iSister),1);

        %move corresponding coordinates
        sisterList(iSister).coords1Aligned = sisterListOriginal(indx1Rand(...
            2*iSister-1)).coords1Aligned;
        sisterList(iSister).coords2Aligned = sisterListOriginal(indx1Rand(...
            2*iSister)).coords1Aligned;

        %calculate new distance vector and store
        sisterList(iSister).distanceAligned = sqrt(sum((sisterList(...
            iSister).coords1Aligned(:,1:3) - sisterList(iSister...
            ).coords2Aligned(:,1:3)).^2,2));

    end
    for iSister = numSisters/2+1 : numSisters

        %assign new sister pair
        sisterList(1).trackPairs(iSister,1) = sisterListOriginal(...
            1).trackPairs(indx2Rand(2*(iSister-numSisters/2)-1),2);
        sisterList(1).trackPairs(iSister,2) = sisterListOriginal(...
            1).trackPairs(indx2Rand(2*(iSister-numSisters/2)),2);

        %move corresponding coordinates
        sisterList(iSister).coords1Aligned = sisterListOriginal(indx2Rand(...
            2*(iSister-numSisters/2)-1)).coords2Aligned;
        sisterList(iSister).coords2Aligned = sisterListOriginal(indx2Rand(...
            2*(iSister-numSisters/2))).coords2Aligned;

        %calculate new distance vector and store
        sisterList(iSister).distanceAligned = sqrt(sum((sisterList(...
            iSister).coords1Aligned(:,1:3) - sisterList(iSister...
            ).coords2Aligned(:,1:3)).^2,2));

    end

end %(if randomize)

%% plots

if verbose
    
    %get project name
    fileName = dataStruct.projectName;
    
    %get number of frames and time lapse
    numFrames = dataStruct.dataProperties.movieSize(end);
    timeLapse = dataStruct.dataProperties.timeLapse;
    
    for iSister = 1 : numSisters
        
        %open figure and write title
        figFileName = [fileName ' - Sister Pair ' num2str(iSister)];
        figHandle = figure('Name',figFileName,'NumberTitle','off');
        
        %get sister coordinates along normal to metaphase plate
        coord1 = sisterList(iSister).coords1Aligned(:,[1 4]);
        coord2 = sisterList(iSister).coords2Aligned(:,[1 4]);
        
        %calculate center position along normal
        coordMean = mean([coord1(:,1) coord2(:,1)],2);
        coordMean(:,2) = 0.5 * sqrt( coord1(:,2).^2 + coord2(:,2).^2 );
        
        %calculate sister separation as projected on normal
        coordDiff = coord2(:,1) - coord1(:,1);
        coordDiff(:,2) = sqrt( coord1(:,2).^2 + coord2(:,2).^2 );
        
        %calculate full sister separation
        sisterVec = sisterList(iSister).coords1Aligned(:,1:3) - ...
            sisterList(iSister).coords2Aligned(:,1:3);
        sisterVecVar = sisterList(iSister).coords1Aligned(:,4:6).^2 + ...
            sisterList(iSister).coords2Aligned(:,4:6).^2;
        sisterDist = sqrt(sum(sisterVec.^2,2));
        sisterDistStd = sqrt( sum( sisterVecVar .* sisterVec.^2 ,2) ) ...
            ./ sisterDist;
        
        %plot positions
        subplot(2,1,1)
        hold on
        plot((0:numFrames-1)*timeLapse,coord1(:,1),'r')
        myErrorbar((0:numFrames-1)*timeLapse,coord1(:,1),coord1(:,2))
        plot((0:numFrames-1)*timeLapse,coord2(:,1),'g')
        myErrorbar((0:numFrames-1)*timeLapse,coord2(:,1),coord2(:,2))
        plot((0:numFrames-1)*timeLapse,coordMean(:,1),'k')
        myErrorbar((0:numFrames-1)*timeLapse,coordMean(:,1),coordMean(:,2))
        legend({'Left sister','Right sister','Center'})
        xlabel('Time (s)')
        ylabel('Position along normal (um)')
        
        %plot sister separation
        subplot(2,1,2)
        hold on
        plot((0:numFrames-1)*timeLapse,sisterDist,'b')
        myErrorbar((0:numFrames-1)*timeLapse,sisterDist,sisterDistStd)
        plot((0:numFrames-1)*timeLapse,coordDiff(:,1),'m')
        myErrorbar((0:numFrames-1)*timeLapse,coordDiff(:,1),coordDiff(:,2))
        legend({'Full','Projection'})
        xlabel('Time (s)')
        ylabel('Sister separation (um)')
        
    end
    
end

%% ~~~ the end ~~~