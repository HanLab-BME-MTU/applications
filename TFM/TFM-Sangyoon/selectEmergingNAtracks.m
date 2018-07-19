function [tracksNA, idx] = selectEmergingNAtracks(tracksNA)
% selectEmergingNAtracks filters out already matured tracks and give you
% only emerging nascent adhesion tracks. 
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
tracksNA = tracksNA(idx);
end
