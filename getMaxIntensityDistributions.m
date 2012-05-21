function [maxIntDist maxIntDistCat] = getMaxIntensityDistributions(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('lb', [3:10 11 16 21 41 61 81 101 141]);
ip.addOptional('ub', [3:10 15 20 40 60 80 100 140 200]);
ip.addParamValue('Tlft_f', 4); % lifetime threshold, in frames
ip.addParamValue('XTicks', 0:40:360);
ip.parse(varargin{:});
lb = ip.Results.lb;
ub = ip.Results.ub;
xa = ip.Results.XTicks;

da = xa(2)-xa(1);

dxi = da/8;
xi = 0:dxi:xa(end)+da;

nc = numel(lb);
nd = numel(data);
mCh = strcmp(data(1).source, data(1).channels);

firstN = 3:20;

maxIntDist = struct([]);
for i = 1:nd
    load([data(i).source 'Tracking' filesep 'ProcessedTracks.mat']);
    tracks = tracks([tracks.catIdx]==1);
    
    lft = [tracks.lifetime_s];
    
    % concatenate amplitudes of master channel into matrix
    trackLengths = [tracks.end]-[tracks.start]+1;
    nt = numel(tracks);
    M = NaN(nt, data(i).movieLength);
    for k = 1:nt
        M(k,1:trackLengths(k)) = tracks(k).A(mCh,:);
    end
    
    for k = 1:nc
        % indexes within cohorts
        cidx = lb(k)<=lft & lft<=ub(k);
        
        maxIntDist(i).maxA{k} = nanmax(M(cidx,:),[],2);
        
        for n = firstN
            maxIntDist(i).(['maxA_f' num2str(n)]){k} = nanmax(M(cidx,1:n),[],2);
        end
        
        % lifetimes for given cohort
        maxIntDist(i).lft{k} = lft(cidx);
        
    end
    maxIntDist(i).lft_all = lft;
    maxIntDist(i).maxA_all = nanmax(M,[],2);
    maxIntDist(i).firstN = firstN;
end


% Rescale EDFs (correction for FP-fusion expression level)
a = rescaleEDFs({maxIntDist.maxA_all}, 'Display', false);

% apply scaling
for i = 1:nd
    maxIntDist(i).maxA = cellfun(@(x) a(i)*x, maxIntDist(i).maxA, 'UniformOutput', false);
    for n = firstN
        fname = ['maxA_f' num2str(n)];
        maxIntDist(i).(fname) = cellfun(@(x) a(i)*x, maxIntDist(i).(fname), 'UniformOutput', false);
    end
end



% Concatenate
maxIntDistCat = struct([]);
for k = 1:numel(lb)
    tmp = vertcat(maxIntDist.maxA);
    maxIntDistCat(k).maxA = vertcat(tmp{:,k});
    
    ni = hist(maxIntDistCat(k).maxA, xi);
    maxIntDistCat(k).ni = ni/sum(ni)/dxi;
    
    for n = firstN
        fname = ['maxA_f' num2str(n)];
        tmp = vertcat(maxIntDist.(fname));
        maxIntDistCat(k).(fname) = vertcat(tmp{:,k});
        
        ni = hist(maxIntDistCat(k).(fname), xi);
        maxIntDistCat(k).(['ni_f' num2str(n)]) = ni/sum(ni)/dxi;
    end
    
    tmp = vertcat(maxIntDist.lft);
    maxIntDistCat(k).lft = horzcat(tmp{:,k});
    
    %mu = mean(pres(k).maxA);
    %mu3 = mean((maxA-mu).^3);
    %mu4 = mean((maxA-mu).^4);
    %sigma = std(pres(k).maxA);
    %pres(k).skew = mu3/sigma^3;
    %pres(k).kurt = mu4/sigma^4-3;
    
    if lb(k)==ub(k)
        maxIntDistCat(k).cohortLabel = num2str(lb(k));
    else
        maxIntDistCat(k).cohortLabel = [num2str(lb(k)) ' - ' num2str(ub(k))];
    end
end




