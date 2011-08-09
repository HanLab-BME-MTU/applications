%M = catTrackFields(data, tracks, fieldname, ch)

% Francois Aguet, 06/24/11

function [M starts ends] = catTrackFields(tracks, movieLength, fieldname, ch)

if nargin<4
    ch = 1;
end

nt = length(tracks);

if ~any(arrayfun(@(t) iscell(t.x), tracks))
    
    M = NaN(nt, movieLength);
    
    starts = [tracks.start];
    ends = [tracks.end];
    
    for t = 1:nt
        M(t, starts(t):ends(t)) = tracks(t).(fieldname)(ch,:);
    end
elseif all(arrayfun(@(t) iscell(t.x), tracks))
    ns = arrayfun(@(t) length(t.x), tracks);
    M = NaN(sum(ns), movieLength);
    i = 1;
    starts = NaN(sum(ns),1);
    ends = NaN(sum(ns),1);
    for t = 1:nt
        for s = 1:ns(t)
            starts(i) = tracks(t).segmentStarts{s};
            ends(i) = tracks(t).segmentEnds{s};
            M(i, starts(i):ends(i)) = tracks(t).(fieldname){s}(ch,:);
            i = i+1;
        end
    end
end