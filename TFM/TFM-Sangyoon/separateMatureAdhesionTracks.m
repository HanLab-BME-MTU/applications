function [tracksNA,tracksNAfailing,tracksNAmaturing,lifeTimeNAfailing,lifeTimeNAmaturing,maturingRatio] = separateMatureAdhesionTracks(tracksNA, outputPath)
% [tracksNA,lifeTimeNA] = separateMatureAdhesionTracks
% separates failing and maturing NA tracks from existing tracksNA, obtain life time of each NA tracks

% Sangyoon Han April 2014

% Set up the output file path
dataPath = [outputPath filesep 'data'];
if ~exist(dataPath,'dir') 
    mkdir(dataPath);
end
%% Lifetime analysis
minLifetime = 5;
maxLifetime = 61;
p=0;
idx = false(numel(tracksNA),1);
for k=1:numel(tracksNA)
    % look for tracks that had a state of 'BA' and become 'NA'
    firstNAidx = find(strcmp(tracksNA(k).state,'NA'),1,'first');
    firstFCidx = find(strcmp(tracksNA(k).state,'FC'),1,'first');
    firstFAidx = find(strcmp(tracksNA(k).state,'FA'),1,'first');
    % see if the state is 'BA' before 'NA' state
%     if (~isempty(firstNAidx) && firstNAidx>1 && strcmp(tracksNA(k).state(firstNAidx-1),'BA')) %%|| (~isempty(firstNAidx) &&firstNAidx==1)
    if ~isempty(firstNAidx) && (isempty(firstFCidx) || firstNAidx<firstFCidx) && (isempty(firstFAidx) || firstNAidx<firstFAidx)
        p=p+1;
        idx(k) = true;
        tracksNA(k).emerging = true;
        tracksNA(k).emergingFrame = firstNAidx;
    else
        tracksNA(k).emerging = false;
    end        
end
%% Analysis of those whose force was under noise level: how long does it take
% Analysis shows that force is already developed somewhat compared to
% background. 
% Filter out any tracks that has 'Out_of_ROI' in their status (especially
% after NA ...)
trNAonly = tracksNA(idx);
indMature = false(numel(tracksNA),1);
indFail = false(numel(tracksNA),1);
p=0; q=0;

for k=1:numel(tracksNA)
    if tracksNA(k).emerging 
        % maturing NAs
        if (any(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'FC')) || ...
                any(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'FA'))) && ...
                sum(tracksNA(k).presence)>8
            
            tracksNA(k).maturing = 1;
            indMature(k) = true;
            p=p+1;
            % lifetime until FC
            lifeTimeNAmaturing(p) = sum(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'NA'));
            % it might be beneficial to store amplitude time series. But
            % this can be done later from trackNAmature
        elseif sum(tracksNA(k).presence)>minLifetime && sum(tracksNA(k).presence)<maxLifetime 
        % failing NAs
            tracksNA(k).maturing = 0;
            indFail(k) = true;
            q=q+1;
            % lifetime until FC
            lifeTimeNAfailing(q) = sum(strcmp(tracksNA(k).state(tracksNA(k).emergingFrame:end),'NA'));
        else
            tracksNA(k).maturing = 2; % it didn't mature for a long time as a NA
        end
    else % this means they are already FC or FA, 
        tracksNA(k).maturing = 3; % already matured at the starting of the movie
    end
end
maturingRatio = p/(p+q);
tracksNAmaturing = tracksNA(indMature);
tracksNAfailing = tracksNA(indFail);
save([dataPath filesep 'failingMaturingTracks.mat'], 'trNAonly', 'tracksNA', 'tracksNAfailing','tracksNAmaturing','maturingRatio','lifeTimeNAfailing','lifeTimeNAmaturing')

end