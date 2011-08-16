%M = catTrackFields(data, tracks, fieldname, ch) concatenates a selected field from all tracks into a matrix

% Francois Aguet, 06/24/11 (Last modified: 08/12/11)

function [M segStarts segEnds startIndex endIndex] = catTrackFields(tracks, movieLength, fieldname, ch)

if nargin<4
    ch = 1;
end

nt = length(tracks);

% tracks must be sorted as regular followed by merging/splitting
idx = find(arrayfun(@(t) ~iscell(t.x), tracks)==1, 1, 'last');
rIdx = 1:idx; % regular tracks
cIdx = (idx+1):nt; % compound tracks

M1 = NaN(idx, movieLength);
starts1 = [tracks(rIdx).start]';
ends1 = [tracks(rIdx).end]';
for t = 1:idx
    M1(t, starts1(t):ends1(t)) = tracks(t).(fieldname)(ch,:);
end

% number of segments in each compound track
ns = arrayfun(@(t) length(t.x), tracks(cIdx));
M2 = NaN(sum(ns), movieLength);

i = 1;
starts2 = NaN(sum(ns),1);
ends2 = NaN(sum(ns),1);
for t = 1:length(cIdx)
    for s = 1:ns(t)
        ti = tracks(cIdx(t));
        starts2(i) = ti.segmentStarts{s};
        ends2(i) = ti.segmentEnds{s};
        M2(i, starts2(i):ends2(i)) = ti.(fieldname){s}(ch,:);
        i = i+1;
    end
end

M = [M1; M2];
segStarts = [starts1; starts2];
segEnds = [ends1; ends2];
if ~isempty(ns)
    si2 = cumsum([1 ns(1:end-1)]);
    ei2 = si2+ns-1;
else
    si2 = [];
    ei2 = [];
end
if ~isempty(rIdx)
    startIndex = [rIdx rIdx(end)+si2]';
    endIndex = [rIdx rIdx(end)+ei2]';
else
    startIndex = si2';
    endIndex = ei2';
end
