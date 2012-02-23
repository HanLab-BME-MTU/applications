
function [tracksAcc tracksRej] = runTrackPostProcessing(data, tracks)

% retain only tracks with lifetime >= 4 frames
lft = round([tracks.lifetime_s] / data.framerate);
tracks = tracks(lft>=4);


validGaps = arrayfun(@(t) max([t.gapStatus 4]), tracks)==4;
singleIdx = [tracks.nSeg]==1;
vis = [tracks.visibility];

idx_Ia = singleIdx & validGaps & vis==1;

tracksX = tracks(idx_Ia);
% kLevel = norminv(1-0.05/2, 0, 1);

N = numel(tracksX);
% Eliminate all tracks with significant signal in buffers
invIdx = zeros(1,N);
for k = 1:N
    %bgcorr = nanmean(tracks(k).c);
    
    %bg = [tracksX(k).startBuffer.sigma_r tracksX(k).endBuffer.sigma_r]*kLevel;
    %A = [tracksX(k).startBuffer.A tracksX(k).endBuffer.A];
    %A_pstd = [tracksX(k).startBuffer.A_pstd tracksX(k).endBuffer.A_pstd];
    
    pval_Ar = [tracksX(k).startBuffer.pval_Ar tracksX(k).endBuffer.pval_Ar];
    hval_Ar = pval_Ar < 0.05;
    % if hval_Ar == 1, reject A<=r, i.e., significant signal
    if any(hval_Ar)
        invIdx(k) = 1;
    end
end
tracksAcc = tracksX(invIdx==0);
tracksRej = tracksX(invIdx==1);


% Reset .isPSF field based on mask size
load([data.source 'Detection' filesep 'detection_v2.mat']);

[f_ecdf, x_ecdf] = ecdf([frameInfo.mask_Ar]);
T = interp1(f_ecdf(2:end), x_ecdf(2:end), 0.99);

for k = 1:numel(tracksAcc)
    tmp = tracksAcc(k).maskN < T;
    tmp(isnan(tmp)) = false;
    tmp2 = tracksAcc(k).isPSF;
    tmp2(isnan(tmp2)) = false;
    tracksAcc(k).isPSF = tmp2 | tmp;
end







% disp('a');




% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('data', @isstruct);
% ip.addParamValue('Buffer', 5, @isscalar);
% ip.addParamValue('FileName', 'trackAnalysis.mat', @ischar);
% ip.addParamValue('Tracks', []);
% ip.parse(data, varargin{:});
% buffer = ip.Results.Buffer;
% tracks = ip.Results.Tracks;
% 
% mCh = find(strcmp(data.source, data.channels));
% nc = numel(data.channels);
% 
% % Thresholds
% T1 = norminv(1-0.05/2.0, 0, 1); % ~2 std above background
% T2 = norminv(1-0.1/2.0, 0, 1);
% 
% if isempty(tracks)
%     load([data.source 'Tracking' filesep ip.Results.FileName]);
% end
% 
% 
% gapStatusAll = arrayfun(@(t) max([t.gapStatus{:} 4]), tracks);
% 
% segFields = {'t', 'x', 'y', 'A', 'c', 'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd',...
%     'sigma_r', 'SE_sigma_r', 'pval_Ar', 'pval_KS', 'isPSF',...
%     'gapStatus', 'gapStarts', 'gapEnds', 'gapLengths', 'gapVect', 'segmentStarts', 'segmentEnds',...
%     'MSD', 'MSDstd', 'totalDisplacement'};
% 
% bufferFields = {'t', 'x', 'y', 'A', 'c', 'A_pstd', 'sigma_r'};
% 
% % I. Remove short segments (1-2 frames) from compound tracks
% segLengths = arrayfun(@(t) t.segmentEnds - t.segmentStarts + 1, tracks, 'UniformOutput', false);
% singletonIdx = cellfun(@(t) find(t==1), segLengths, 'UniformOutput', false);
% nSingleton = cellfun(@(t) numel(t), singletonIdx);
% idx = find(nSingleton>0 & [tracks.nSeg]>1);
% 
% for k = idx
%     for f = 1:length(segFields)
%         if ~isempty(tracks(k).(segFields{f}))
%             tracks(k).(segFields{f})(singletonIdx{k}) = [];
%         end
%     end
%     
%     for f = 1:length(bufferFields)
%         if ~isempty(tracks(k).startBuffer.(bufferFields{f}))
%             tracks(k).startBuffer.(bufferFields{f})(singletonIdx{k}) = [];
%             tracks(k).endBuffer.(bufferFields{f})(singletonIdx{k}) = [];
%         end
%     end
%     
%     tracks(k).tracksFeatIndxCG(singletonIdx{k},:) = [];
%     sev = tracks(k).seqOfEvents;
%     for i = singletonIdx{k}
%         sev(sev(:,3) == i,:) = [];
%     end
%     tracks(k).seqOfEvents = sev;
%     tracks(k).nSeg = tracks(k).nSeg - nSingleton(k);
% end
% 
% ns = [tracks.nSeg];
% 
% 
% % II. Split compound tracks with 2-3 segments, if spatially isolated
% idx = find((ns==2 | ns==3) & [tracks.status]==1 & [tracks.start]>=1+buffer & [tracks.end]<=data.movieLength-buffer & gapStatusAll==4);
% for k = 1:length(idx)
%     ti = idx(k);
%     
%     % segment mask
%     %segmask = ~isnan(catTrackFields(tracks(ti), data.movieLength, 'x', 1));
%     
%     % for all segments, test: start buffer,end buffer
%     A1 = vertcat(tracks(ti).startBuffer.A{:});
%     sr1 = vertcat(tracks(ti).startBuffer.sigma_r{:});
%     A2 = vertcat(tracks(ti).endBuffer.A{:});
%     sr2 = vertcat(tracks(ti).endBuffer.sigma_r{:});
%     validSeg = sum(A1 < T2*sr1,2) >= 3 & sum(A2 < T2*sr2,2) >= 3;
%     validSeg = validSeg(mCh:nc:end);
%     %nv = numel(validSeg);
%     
%     starts = tracks(ti).segmentStarts;
%     ends = tracks(ti).segmentEnds;
%     % for each valid segment, check whether it is spatially isolated in overlapping regions
%     for si = find(validSeg)
%         %remIdx = setdiff(find(validSeg), si);
%         remIdx = setdiff(1:ns(ti), si);
%         
%         overlapStarts = max(starts(si), starts(remIdx));
%         overlapEnds = min(ends(si), ends(remIdx));
%         %overlap = overlapEnds - overlapStarts + 1;
%         
%         
%         
%     end
%     
%     
% end


% overlap+5 buffer frames

% lifetimes_s = [tracks.lifetime_s];



% II. Cut tracks corresponding to multiple sequential events into singletons


% III. Remove tracks that are not diffraction-limited

