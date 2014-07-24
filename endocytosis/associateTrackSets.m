% Francois Aguet, May 2011

function [finalIdx noAssocIdx] = associateTrackSets(trackSet1, trackSet2, varargin)

mCh = 1;

% mean positions of tracks
mean_set1 = arrayfun(@(t) [mean(t.x(mCh,:)) mean(t.y(mCh,:))], trackSet1, 'UniformOutput', false);
mean_set1 = vertcat(mean_set1{:});

mean_set2 = arrayfun(@(t) [mean(t.x(mCh,:)) mean(t.y(mCh,:))], trackSet2, 'UniformOutput', false);
mean_set2 = vertcat(mean_set2{:});

% 5 pixel search radius; indexes returned from set1
[idx, d] = KDTreeBallQuery(mean_set1, mean_set2, 8*ones(size(mean_set2)));

idxM = [arrayfun(@(x) x, 1:numel(idx), 'UniformOutput', false)' idx];

% remove tracks with no association
noAssocIdx = find(cellfun(@(x) isempty(x), idx));
idxM(noAssocIdx,:) = [];

% tracks in set2 associated to set1
% assocIdx_set2 = find(~cellfun(@(x) isempty(x), idx));

% tracks in set2 not associated to set1
% noAssocIdx_set2 = cellfun(@(x) isempty(x), idx);
% lengths of these tracks
% l2 = arrayfun(@(t) length(t.t), trackSet2(noAssocIdx_set2));
% figure; hist(l2)

for ti = 1:size(idxM,1)
    % matches
    matches = idxM{ti,2};
    nm = length(matches);
    validOverlap = zeros(1,nm);
    for m = 1:nm
        validOverlap(m) = max([trackSet2(idxM{ti,1}).start trackSet1(matches(m)).start]) < ...
                          min([trackSet2(idxM{ti,1}).end trackSet1(matches(m)).end]);
    end
    idxM{ti,2}(validOverlap==0) = [];
end
rmIdx = cellfun(@(x) isempty(x), idxM(:,2));

idxM(rmIdx,:) = [];


finalIdx = idxM;



% 
% % retain only matches with overlap
% % finalIdx = cell(length(assocIdx_set2),2);
% for ti = 1:length(assocIdx_set2)
%     matches = idx{assocIdx_set2(ti)}; % indexes in set1
%     nm = length(matches);
%     validOverlap = zeros(1,nm);
%     for m = 1:nm
%         validOverlap(m) = max(trackSet2(ti).start, trackSet1(matches(m)).start) <...
%                           min(trackSet2(ti).end, trackSet1(matches(m)).end);
%     end
%     finalIdx{ti,1} = assocIdx_set2(ti);
%     finalIdx{ti,2} = matches(validOverlap==1);
% end
% finalIdx(cellfun(@(x) isempty(x), finalIdx(:,2)),:) = [];

% outputs:
% unmatched