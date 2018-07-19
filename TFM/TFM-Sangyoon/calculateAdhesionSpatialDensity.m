function [tracksNA] = calculateAdhesionSpatialDensity(tracksNA, r)
% [tracksNA] = calculateAdhesionSpatialDensity(tracksNA, r)
% calculates adhesion density within vicinity radius r for entire tracks
% with respect to the mean position of each track.
% input
%           tracksNA:           tracks that have member of .maturing
%           (created by separateMatureAdhesionTracks after colocalizationAdhesionsWithTFM)
%           r:                       neighbor radius
%           
% output
%           tracksNA with new members such as:
%                   .numNeighbors:            the number of adhesions
%                   .rNeigh:                       vicinity radius used for neighbor quantification  
%                   .adhDensity                 adhesion density = numNeighbors/(rNeigh^2*pi)
%           
% Sangyoon Han April 2015

%% Neighbor sample statistics - failing adhesions
% for failing adhesions
meanXall = arrayfun(@(x) nanmean(x.xCoord),tracksNA);
meanYall = arrayfun(@(x) nanmean(x.yCoord),tracksNA);
for ii = 1:numel(tracksNA)
    % find the neighbors in r
    [idx, ~] = KDTreeBallQuery([meanXall, meanYall],...
        [nanmean(tracksNA(ii).xCoord),nanmean(tracksNA(ii).yCoord)],r);
%         [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(NAFrames),tracksNA(presentID)),...
%             arrayfun(@(x) x.yCoord(NAFrames),tracksNA(presentID))],...
%             [tracksNA(ii).xCoord(NAFrames),tracksNA(ii).yCoord(NAFrames)],r);
    if ~isempty(idx{1})
        % inspect their fates
%         nNei = length(idx{1});
        neiIDs = idx{1};
        tracksNA(ii).numNeighbors = length(neiIDs);
        tracksNA(ii).rNeigh = r;
        tracksNA(ii).adhDensity =  tracksNA(ii).numNeighbors/(r^2*pi); % in #/pixel
    else
        tracksNA(ii).numNeighbors = 0;
        tracksNA(ii).rNeigh = r;
        tracksNA(ii).adhDensity =  0; % in #/pixel
    end
end

end