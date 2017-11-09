function [tracksNA,idx,idxFC] = filterOutNonEmergingNA(tracksNA)
%  [tracksNA] = filterTracksNA(tracksNA)
% filter out tracks whose state is 'NA' without 'BA' in the previous time
% point(s).

% Sangyoon Han September 2015
%% Pre-processing
% re-status the NA based on startingFrameExtra and endingFrameExtra
disp('Pre-processing')
tic
for k=1:numel(tracksNA)
    for ii=tracksNA(k).startingFrameExtra:tracksNA(k).endingFrameExtra
        if strcmp(tracksNA(k).state(ii),'BA') || strcmp(tracksNA(k).state{ii},'ANA')
            tracksNA(k).state{ii} = 'NA';
        end
    end
end
toc
%% filtering
idx = false(numel(tracksNA),1);
idxFC = false(numel(tracksNA),1); % When adhesion already starts as FA or FC

for k=1:numel(tracksNA)
    % look for tracks that had a state of 'BA' and become 'NA'
    firstNAidx = find(strcmp(tracksNA(k).state,'NA'),1,'first');
    firstFCidx = find(strcmp(tracksNA(k).state,'FC'),1,'first');
    firstFAidx = find(strcmp(tracksNA(k).state,'FA'),1,'first');
    firstAdhIdx = min([firstNAidx firstFCidx firstFAidx]);
    % see if the state is no 'BA' before any adhesion state
    if (~isempty(firstNAidx) && firstAdhIdx>1 && strcmp(tracksNA(k).state(firstAdhIdx-1),'BA')) %%|| (~isempty(firstNAidx) &&firstNAidx==1)
%     if (~isempty(firstNAidx) && firstNAidx==1) || (~isempty(firstNAidx) && firstNAidx>1 && ~strcmp(tracksNA(k).state(firstNAidx-1),'BA')) %%|| (~isempty(firstNAidx) &&firstNAidx==1)
%     if ~isempty(firstNAidx) && (isempty(firstFCidx) || firstNAidx<firstFCidx) && (isempty(firstFAidx) || firstNAidx<firstFAidx)
        idx(k) = true;
    end    
    if (~isempty(firstFCidx) && firstFCidx==1) || (~isempty(firstFAidx) && firstFAidx==1)
        idxFC(k)=true;
    end
end
%% Analysis of those whose force was under noise level: how long does it take
% Analysis shows that force is already developed somewhat compared to
% background. 
% Filter out any tracks that has 'Out_of_ROI' in their status (especially
% after NA ...)
tracksNA = tracksNA(idx);