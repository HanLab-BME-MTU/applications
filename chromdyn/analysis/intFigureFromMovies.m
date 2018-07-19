function intFigureFromMovies
%INTFIGUREFROMMOVIES produces an intFigure from movies similar to the one generated for spindle images with measureSpindles

% this is a one-time hack

% read idlists, dataProperties. Only if 3+ tags, no '?'.
idlistList = loadIdlistList(cdBiodata(4),...
    'length(idlist(1).stats.labelcolor) > 2 && isempty(strmatch(''?'',idlist(1).stats.labelcolor)) ');

nData = length(idlistList);

% loop through the data. For every 9 frames, read 
% [nSpots, ds1s2, mean(cenInt)/mean(spbInt)]

data1 = NaN(30*nData,3); % use measureSpindles-nomenclature
ct = 1;

for iData = 1:nData
    % find order of tags
    spbIdx = strmatch('spb',idlistList(iData).idlist(1).stats.labelcolor);
    cenIdx = strmatch('cen',idlistList(iData).idlist(1).stats.labelcolor);
    % loop through linklists 1:9:end. Check if linklist exists, check
    % nTags, fill data1
    for t=1:9:length(idlistList(iData).idlist)
        if ~isempty(idlistList(iData).idlist(t).linklist)
            ll = idlistList(iData).idlist(t).linklist;
            % mask all tracked/estimated spots. We want this to be an exact
            % comparison!
            ll(ll(:,3)>0|ll(:,5)>0,8)=NaN;
            spbAmp = ll(spbIdx,8);
            cenAmp = ll(cenIdx,8);
            if any(isfinite(cenAmp)) && all(isfinite(spbAmp))
                data1(ct,:) = ...
                    [sum(isfinite(cenAmp))+2,...
                    norm(diff(ll(spbIdx,9:11))),...
                    nanmean(cenAmp)/nanmean(spbAmp)];
                ct = ct+1;
            end
        end
    end
end

% remove superfluous entries in data1
data1(ct:end,:) = [];

% plot figure
d31 = data1(:,1)==3 & data1(:,1) > 1;
d41 = data1(:,1) == 4;

figure('Name','relative centromere intensity (3=b/4=r)')
plot(data1(d31,2),data1(d31,3),'.b',data1(d41,2),data1(d41,3),'.r')
xlabel('Spindle Length (\mum)')
ylabel('$$\frac{\Sigma(cenInt)}{\Sigma(spbInt)}$$','Interpreter','latex')