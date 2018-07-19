function z = nanZscore(vec)
% jungsik noh
% 2016/07/11
    m = nanmean(vec);
    s = nanstd(vec);
    if s == 0
        z = nan(size(vec));
    else
        z = (vec - m*ones(size(vec))) ./s;
    end
end

