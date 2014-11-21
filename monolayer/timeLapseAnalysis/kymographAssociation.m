%% Calculates temporal order between two physical properties, based on spatiotemporal kymographs
%   1. gather all corresponding kymogrpahs bins according to different time lags
%   2. calculate and output correlations + max correlation + time lag

function [correlations] = kymographAssociation(kymographs1,kymogrpahs2,timeShifts,nTimeFrame)

nKymographs = length(kymographs1);

% Sanity check
assert(nKymographs == length(kymogrpahs2));

nTimeShifts = length(timeShifts);



correlations = nan(1,nTimeShifts);
pvals = zeros(1,nTimeShifts);

for d = 1 : nTimeShifts
    shift = timeShifts(d);
    kemos1 = [];
    kemos2 = [];
    for k = 1 : nKymographs
        kymograph1 = kymographs1{k};
        kymograph2 = kymogrpahs2{k};
        kymograph1 = kymograph1(1:12,1:nTimeFrame);
        kymograph2 = kymograph2(1:12,1:nTimeFrame);
        kymograph1 = imresize(kymograph1,[size(kymograph1,1)/3,size(kymograph1,2)]);
        kymograph2 = imresize(kymograph2,[size(kymograph2,1)/3,size(kymograph2,2)]);
        if shift < 0
            kymo1 = kymograph1(:,abs(shift)+1:end);
            kymo2 = kymograph2(:,1:end-abs(shift));
        else
            kymo1 = kymograph1(:,1:end-shift);
            kymo2 = kymograph2(:,abs(shift)+1:end);
        end
        naninds = isnan(kymo1) | isnan(kymo2) | isinf(kymo1) | isinf(kymo2);
        kymo1 = kymo1(~naninds);
        kymo2 = kymo2(~naninds);
        kymos1 = [kemos1; kymo1];
        kymos2 = [kemos2; kymo2];
    end
    [rho, pval] = corr(kymos1,kymos2);
    pvals(d) = pval;
    if pval < 0.05 %&& rho > 0.3
        correlations(d) = rho;
    end
end

inds = find(~isnan(correlations));
if ~enoughInARow(inds)
    correlations = nan(1,nTimeShifts);
end

end

function [yes] = enoughInARow(inds)
yes = false;

if isempty(inds)
    return;
end

nMax = 1;
nCur = 1;

iCur = 1;
while iCur < length(inds)
    while iCur+1 <= length(inds) && inds(iCur+1) == inds(iCur) + 1
        nCur = nCur + 1;
        iCur = iCur + 1;
    end
    if nCur > nMax
        nMax = nCur;
    end
    iCur = iCur + 1;
end

yes = nCur > 2; 

end