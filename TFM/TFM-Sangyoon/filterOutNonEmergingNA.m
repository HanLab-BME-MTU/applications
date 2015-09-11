function [tracksNA] = filterOutNonEmergingNA(tracksNA)
%  [tracksNA] = filterTracksNA(tracksNA)
% filter out tracks whose state is 'NA' without 'BA' in the previous time
% point(s).

% Sangyoon Han September 2015
%% filtering
idx = true(numel(tracksNA),1);
for k=1:numel(tracksNA)
    % look for tracks that had a state of 'BA' and become 'NA'
    firstNAidx = find(strcmp(tracksNA(k).state,'NA'),1,'first');
    % see if the state is no 'BA' before 'NA' state
    if (~isempty(firstNAidx) && firstNAidx==1) || (~isempty(firstNAidx) && firstNAidx>1 && ~strcmp(tracksNA(k).state(firstNAidx-1),'BA')) %%|| (~isempty(firstNAidx) &&firstNAidx==1)
%     if ~isempty(firstNAidx) && (isempty(firstFCidx) || firstNAidx<firstFCidx) && (isempty(firstFAidx) || firstNAidx<firstFAidx)
        idx(k) = false;
    end        
end
%% Analysis of those whose force was under noise level: how long does it take
% Analysis shows that force is already developed somewhat compared to
% background. 
% Filter out any tracks that has 'Out_of_ROI' in their status (especially
% after NA ...)
tracksNA = tracksNA(idx);