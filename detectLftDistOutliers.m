% Francois Aguet, 01/28/2012

function lftDataOut = detectLftDistOutliers(lftData, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('lftData');
ip.addParamValue('EndIdx', find(lftData.meanHist~=0, 1, 'last'));
ip.addParamValue('Display', true, @islogical);
ip.parse(lftData, varargin{:});
endIdx = ip.Results.EndIdx;

N = numel(lftData.lftHist);

ecdfMat = cumsum(vertcat(lftData.lftHist{:}), 2);
ecdfMat = ecdfMat(:,1:endIdx);
t = lftData.t(1:endIdx);

A = ecdfMat(:,1);
for i = 1:numel(A)
    ecdfMat(i,:) = (ecdfMat(i,:)-A(i))/(1-A(i));
end



%============================================
% Plot raw distributions
%============================================
figure; hold on;
fset = loadFigureSettings();
% set(gca, 'ColorOrder', C);
plot(t, ecdfMat', 'k', 'LineWidth', 1);
axis([0 100 0 1.01]);
set(gca, fset.sfont{:}, 'LineWidth', 2, 'Layer', 'top',...
    'YTick', 0:0.1:1, 'XTick', 0:10:lftData.t(end));
xlabel('Time (s)');
ylabel('F(t)');


cmb = pcombs(1:N);
nc = size(cmb,1);
D = zeros(1,nc);
for i = 1:nc
    D(i) = max(abs(ecdfMat(cmb(i,1),:)-ecdfMat(cmb(i,2),:)));
end

%============================================
% Identification of lifetime subpopulations (clusters) 
%============================================

% check whether a two-component mixture fits the data
% obj = gmdistribution.fit(D(:),2, 'replicates', 20);
% mu = obj.mu;
% sigma = sqrt(squeeze(obj.Sigma));
% if abs(diff(mu)) > 2*sum(sigma)
%     k = 2;
% else
%     k = 1;
% end

ns = sqrt(mean(lftData.nSamples)/2);
k = 0;
isnormal = 0;
while k<=2 && ~all(isnormal)
    k = k+1;
    
    [idx] = kmeans(D, k, 'replicates', 20);
    % test each cluster for normality
    isnormal = NaN(1,k);
    for c = 1:k
        dstr = D(idx==c);
        
        [f_ecdf, t_ecdf] = ecdf(ns*dstr);
        
        % Kolmogorov CDF
        ni = 50;
        i = 1:ni;
        nx = numel(t_ecdf);
        xM = repmat(t_ecdf', [ni 1]);
        iM = repmat(i(:), [1 nx]);        
        M = exp(-(2*iM-1).^2*pi^2./(8*xM.^2))*sqrt(2*pi)./xM;
        F = sum(M,1);
        
        D_ks = max(abs(F-f_ecdf')); % test statistic
        
        pval = qks((sqrt(ns) + 0.12 + 0.11/sqrt(ns))*D_ks);
        hval = pval < 0.05;
        
        % test for normality (usually when numel(dstr) is small (<10))
        if hval==1
            dstr = (dstr-mean(dstr))/std(dstr);
            [hval] = kstest(dstr,[], 0.1);
        end
        %if numel(dstr)>5
        %    hval = adtest(dstr);
        %else
        %    hval = 0; % assume normality
        %end
        isnormal(c) = ~hval;
    end
end
if ~all(isnormal)
    error('More than two characteristic distributions detected.');
end



if k == 2
    [idx ctr] = kmeans(D, k, 'replicates', 20);
    
    minIdx = find(ctr==min(ctr));
    
    % idx for the cluster with smaller distances contains two sets
    cand = cmb(idx==minIdx,:);
    sets{1} = cand(1,:);
    cand(1,:) = [];
    
    newset = 0;
    while ~isempty(newset)
        idx2 = ismember(cand(:,1), sets{1}) | ismember(cand(:,2), sets{1});
        newset = cand(idx2==1,:);
        sets{1} = unique([sets{1} newset(:)']);
        cand(idx2==1,:) = [];
    end
    tmp = unique(cand(:))';
    if ~isempty(tmp)
        sets{2} = tmp;
    end
else
    sets{1} = 1:N;
end

%============================================
% Outliers in each set: median ± 3*sigma
%============================================
nset = numel(sets);

medianECDF = cell(1,nset);
sigma = cell(1,nset);
ub = cell(1,nset);
lb = cell(1,nset);
outlierIdx = cell(1,nset);
inlierIdx = cell(1,nset);
for k = 1:nset
    N = numel(sets{k});
    
    medianECDF{k} = median(ecdfMat(sets{k},:), 1);
    sigma{k} = 1/norminv(0.75) * mad(ecdfMat(sets{k},:), 1, 1);
    %sigma{k} = std(ecdfMat(sets{k},:), [], 1);
    
    ub{k} = medianECDF{k}+3*sigma{k};
    lb{k} = medianECDF{k}-3*sigma{k};
    
    x = 1:numel(t);
    chub = convhull(x, ub{k});
    chub = interp1(x(chub(2:end)), ub{k}(chub(2:end)), x, 'spline');
    %chlb = convhull(x, lb{k});
    %chlb = interp1(x(chlb(2:end)), lb{k}(chlb(2:end)), x, 'spline');
    
    du = chub - medianECDF{k};
    %dl = medianECDF{k} - chlb;
    %dm = (du+dl)/2;
    
%     figure; hold on;
%     plot(x, medianECDF{k}, 'k');
%     plot(x, medianECDF{k}+dm, 'c');
%     plot(x, medianECDF{k}-dm, 'c');
%     plot(x, ub{k}); plot(x, chub, 'r');
%     plot(x, lb{k}); plot(x, chlb, 'r');
    
    ub{k} = medianECDF{k}+du;
    lb{k} = medianECDF{k}-du;
    
    outlierIdx{k} = zeros(1,N);
    for i = 1:N
        tmp = ecdfMat(sets{k}(i),:) > ub{k} | ecdfMat(sets{k}(i),:) < lb{k};
        outlierIdx{k}(i) =  sum(tmp)/numel(tmp) > 0.05;
    end
    outlierIdx{k} = find(outlierIdx{k});
    inlierIdx{k} = setdiff(1:N, outlierIdx{k});
end

if nset>1
    fprintf('%d lifetime distribution clusters identified.\n', nset);
    for i = 1:nset
        if ~isempty(outlierIdx{i})
            fprintf('Outliers in cluster %d: %s\n', i, num2str(outlierIdx{i}));
        end
    end
else
    if ~isempty(outlierIdx{1})
        fprintf('Outliers: %s\n', num2str(outlierIdx{1}));
    end
end


%============================================
% Plot results
%============================================
hp = NaN(1,4);
if ip.Results.Display
    fset = loadFigureSettings();
    figure; hold on;
    for k = 1:numel(sets)
        
        hp(4) = fill([t t(end:-1:1)], [lb{k} ub{k}(end:-1:1)], 0.6*[1 1 1],...
            'EdgeColor', 'none');
    
        if ~isempty(outlierIdx{k})
            hx = plot(t, ecdfMat(sets{k}(outlierIdx{k}),:), 'Color', [1 0 0], 'LineWidth', 2);
            hp(2) = hx(1);
        end
        hx = plot(t, ecdfMat(sets{k}(inlierIdx{k}),:), 'Color', [0 1 0], 'LineWidth', 2);
        hp(1) = hx(1);
        hp(3) = plot(t, medianECDF{k}, 'k--', 'LineWidth', 1.5);
    end
    axis([0 100 0 1.01]);
    set(gca, fset.sfont{:}, 'LineWidth', 2, 'Layer', 'top',...
        'YTick', 0:0.1:1, 'XTick', 0:10:lftData.t(end));
    xlabel('Time (s)');
    ylabel('F(t)');
    if isnan(hp(2))
        hp(2) = [];
        hl = legend(hp, 'Inlier distribution', 'Median distribution', 'Critical area', 'Location', 'SouthEast');
    else
        hl = legend(hp, 'Inlier distribution', 'Outlier distribution', 'Median distribution', 'Critical area', 'Location', 'SouthEast');
    end
    set(hl, fset.tfont{:}, 'Box', 'off');
    pos = get(hl, 'Position');
    pos(2) = 0.2;
    set(hl, 'Position', pos);
end


%============================================
% Create output structures
%============================================
lftDataOut = struct([]);
%lftDataOut(1:nset) = struct('source', [], 't', [], 'lftHist', [], 'meanHist', []);
fnames = fieldnames(lftData);
idx = find(structfun(@(i) iscell(i), lftData));
for k = 1:nset
    lftDataOut(k).t = lftData.t;
    for fi = 1:numel(idx) % loop through cell fields
        f = idx(fi);
        lftDataOut(k).(fnames{f}) = lftData.(fnames{f})(inlierIdx{k});
        if ~isempty(outlierIdx{k})
            lftDataOut(k).outliers.(fnames{f}) = lftData.(fnames{f})(outlierIdx{k});
        end
    end
    lftDataOut(k).meanHist = mean(vertcat(lftDataOut(k).lftHist{:}),1);
end



function q = qks(z)
% if z < 0.0
%     error('z value for KS dist. must be positive.');
% end
if z == 0.0
    q = 1.0;
elseif z < 1.8
    y = exp(-pi^2/(8*z^2));
    q = 1.0 - sqrt(2*pi) / z * (y + y^9 + y^25 + y^49);
else
    x = exp(-2*z^2);
    q = 2.0*(x - x^4 - x^9);
end
