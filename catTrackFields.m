%M = catTrackFields(data, tracks, fieldname, ch)

% Francois Aguet, 06/24/11

function M = catTrackFields(data, tracks, fieldname, ch)

if nargin<4
    ch = 1;
end

nt = length(tracks);

M = NaN(nt, data.movieLength);

starts = [tracks.start];
ends = [tracks.end];

for t = 1:nt
    M(t, starts(t):ends(t)) = tracks(t).(fieldname)(ch,:);
end