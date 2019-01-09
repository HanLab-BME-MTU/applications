function [tracks,idxFinalTracks,idOtherOverlappingTracks] = filterOverlappingTracks(tracks)
%function curStartingAmpG1 = filterOverlappingTracks(curStartingAmpG1)
%filters out overlapping tracks.
% input: tracks      tracks
% output: tracks     tracks filtered
% Sangyoon Han Jan 2017

% Find out overlapping tracks
% Definition of overlapping: whenever there are five frames 
% err = 1e-2;

% Per track, try to find other tracks that overlap with the track
numOverFrames=5;
% pp=0;
% notCompletelyChecked=true;
idxlTracksLeft=true(numel(tracks),1);
idxFinalTracks=false(numel(tracks),1);
idOtherOverlappingTracks=zeros(numel(tracks),1);
allX = arrayfun(@(x) x.xCoord,tracks,'UniformOutput',false);
allY = arrayfun(@(x) x.yCoord,tracks,'UniformOutput',false);
while any(idxlTracksLeft)
%     pp=pp+1;
    curInd = find(idxlTracksLeft,1);
%     curInd = leftIndices(pp);
    curTrack = tracks(curInd);
    curX = curTrack.xCoord;
    curY = curTrack.yCoord;
%     curAmp = curTrack.ampTotal;
    
    % matching
    matchXY = cellfun(@(x,y) round(x,3)==round(curX,3) & round(y,3)==round(curY,3),allX,allY,'UniformOutput',false);
    numMatchXY = cellfun(@sum,matchXY);
    curIdxOver = numMatchXY>numOverFrames;
    % Merging these overlapping tracks based on earliest and latest 
    curOverX = (allX(curIdxOver));
%     curOverY = (allY(curIdxOver));
    
%     figure, plot(curOverX(1,:),curOverY(1,:))
%     hold on, plot(curOverX(2,:),curOverY(2,:))
%     plot(curOverX(3,:),curOverY(3,:))
    earlyX = cellfun(@(x) find(~isnan(x),1),curOverX);
    lateX = cellfun(@(x) find(~isnan(x),1,'last'),curOverX);
    % assuming the overlapping tracks end always at the same frame...
    indexCurOver=find(curIdxOver);
    if length(unique(earlyX))>1
        [~,indEarliest]=min(earlyX);
        if length(unique(lateX))>1
            disp([num2str(indexCurOver') 'th tracks have ' num2str(length(unique(lateX))) ' endings.'])
        end
        idxFinalTracks(indexCurOver(indEarliest))=true;
        idOtherOverlappingTracks(curIdxOver)=(indexCurOver(indEarliest));
    else
        [~,indLatest]=min(lateX);
        if length(unique(earlyX))>1
            disp([num2str(indexCurOver') 'th tracks have ' num2str(length(unique(earlyX))) ' beginnings.'])
        end
        idxFinalTracks(indexCurOver(indLatest))=true;
        idOtherOverlappingTracks(curIdxOver)=(indexCurOver(indLatest));
    end
    idxlTracksLeft(curIdxOver)=false;
end
disp(['Originally ' num2str(numel(tracks)) ' tracks were found to have ' num2str(sum(~idxFinalTracks)) ' overlaps.'])
tracks=tracks(idxFinalTracks);
disp(['Returning ' num2str(numel(tracks)) ' tracks as non-overlapping tracks.'])
end
