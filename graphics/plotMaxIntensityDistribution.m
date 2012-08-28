function plotMaxIntensityDistribution(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Mode', 'pdf', @(x) any(strcmpi(x, {'pdf', 'cdf'})));
ip.addParamValue('XTick', []);
ip.addParamValue('YTick', []);
ip.addParamValue('FirstNFrames', 5);
ip.addParamValue('CohortLB', [5  11 16 21 41 61 81]);
ip.addParamValue('CohortUB', [10 15 20 40 60 80 120]);
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'screen', 'print'})));
ip.addParamValue('ShowSignificance', false, @islogical);
ip.addParamValue('ShowGaussians', false, @islogical);
ip.addParamValue('ShowFirstFrame', false, @islogical);
ip.addParamValue('Cutoff_f', 5);
ip.parse(varargin{:});
mode = ip.Results.Mode;

[~, figPath] = getCellDir(data(1));
figPath = [figPath 'Figures' filesep];
[~,~] = mkdir(figPath);

mCh = strcmp(data(1).source, data(1).channels);

lb = ip.Results.CohortLB;
ub = ip.Results.CohortUB;
nc = numel(lb);

ny = nc;

lftData = getLifetimeData(data);
A = arrayfun(@(i) i.A(i.lifetime_s(i.catIdx==1)>=ip.Results.Cutoff_f,:,mCh), lftData, 'UniformOutput', false);
maxA_all = cellfun(@(i) nanmax(i,[],2)', A, 'UniformOutput', false);

% Rescale EDFs (correction for FP-fusion expression level)
a = rescaleEDFs(maxA_all, 'Display', false);

% apply scaling
nd = numel(data);
for i = 1:nd
    maxA_all{i} = a(i) * maxA_all{i};
    A{i} = a(i)*A{i};
end


% Concatenate maximum intensity and lifetime data
maxA = [maxA_all{:}];
A = vertcat(A{:});
f = ip.Results.FirstNFrames;

maxAFirstN = nanmax(A(:,1:f), [], 2);

lifetime_s = [lftData.lifetime_s];
lifetime_s = lifetime_s([lftData.trackLengths]>=ip.Results.Cutoff_f & [lftData.catIdx]==1);

% x-axis
xa = ip.Results.XTick;
if isempty(xa)
    db = prctile(maxA, [0.1 99.9]);
    mag = 10^floor(log10(db(2)));
    da = ceil(db(2)/mag)*mag/10;
    xa = (0:ceil(db(2)/da))*da;
else
    da = xa(2)-xa(1);
end
dxi = da/8;
xi = 0:dxi:xa(end)+da;

%%
% verify whether scaling of A_pstd is correct
% figure;
% hold on;
% 
% dx = 1;
% xi = 0:dx:max(res(medIdx).maxA_all_pstd);
% ni = hist(res(medIdx).maxA_all_pstd, xi);
% ni = ni/sum(ni)/dx;
% %[ni,xi] = ksdensity(samples{medIdx}, 'npoints', 1000);
% plot(xi, ni, 'r.-', 'LineWidth', 2);
% idx = setdiff(1:nd, medIdx);
% for i = 1:nd-1
%     ni = hist(res(idx(i)).maxA_all_pstd*a(idx(i)), xi);
%     ni = ni/sum(ni)/dx;
%     %[ni,xi] = ksdensity(samples{idx(i)}, 'npoints', 1000);
%     plot(xi, ni, 'k.-', 'LineWidth', 2);
% end


%%
%========================================
% Generate cohorts
%========================================
maxAcohort = cell(1,nc);
maxAcohortFirstN = cell(1,nc);
A0cohort = cell(1,nc);
ni = cell(1,nc);
niFirstN = cell(1,nc);
niFirst = cell(1,nc);
for k = 1:numel(lb)
    cidx = lb(k)<=lifetime_s & lifetime_s<ub(k);
    maxAcohort{k} = maxA(cidx);
    maxAcohortFirstN{k} = maxAFirstN(cidx);
    A0cohort{k} = A(cidx,1);

    ni0 = hist(maxAcohort{k}, xi);
    ni{k} = ni0/sum(ni0);
    ni0 = hist(maxAcohortFirstN{k}, xi);
    niFirstN{k} = ni0/sum(ni0);
    ni0 = hist(A0cohort{k}, xi);
    niFirst{k} = ni0/sum(ni0);
end

if isempty(ip.Results.YTick)
    di = 3;
    ymax = max(cellfun(@(i) max(i), niFirstN));
    mag = 10^floor(log10(ymax/di));
    dy = ceil(ymax/mag/di)*mag;
    ya = (0:di)*dy;
else
    ya = ip.Results.YTick;
end

%%
cf0 = [1 1 1]*0.6;
ce0 = [1 1 1]*0.3;
fset = loadFigureSettings('print');

pos = get(0, 'DefaultFigurePosition');
pos(3:4) = [550 700];

axPos = fset.axPos;
aw = axPos(3);
ah = 0.2*aw;
xo = 1.5;
yo = 1.5;
sh = ah/3;

cmap = jet(nc);

% figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');
figure('Units', 'centimeters', 'Position', [5 5 8 ny*ah+(ny-1)*sh+2.5], 'Color', [1 1 1], 'PaperPositionMode', 'auto');
hbg = axes('Units', 'centimeters', 'Position', [0 0 2*xo+aw+2 2*yo+ny*ah+(ny-1)*sh]);
hold(hbg, 'on');
axis(hbg, [0 2*xo+aw+2 0 2*yo+ny*ah+(ny-1)*sh]);
for k = 1:nc
    
    iPos = [xo yo+(ny-k)*(ah+sh) aw ah];
    % plot grid below data
    hi = axes(fset.axOpts{:}, 'Position', iPos,...
        'XLim', [xa(1) xa(end)], 'XTick', xa, 'XTickLabel', [], 'YLim', [ya(1) ya(end)], 'YTickLabel', [],...
        'TickLength', [0 0], 'Color', 'none');
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':', 'LineWidth', 1);
    
    % data axes
    hi = axes(fset.axOpts{:}, 'Position', iPos,...
        'Color', 'none');
    hold on;
    box off;
    
    if strcmpi(mode, 'pdf')
        
        if ip.Results.ShowFirstFrame
            bar(xi, niFirst{k}, 'BarWidth', 1, 'FaceColor', 0.9*[1 1 1], 'EdgeColor', 0.45*[1 1 1], 'LineWidth', 0.75);
        end
        bar(xi, niFirstN{k}, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 0.75);
        bar(xi, ni{k}, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 0.75);
        stairsXT(xi, niFirstN{k}, 'EdgeColor', ce0, 'LineWidth', 0.75);
        if ip.Results.ShowFirstFrame
            stairsXT(xi, niFirst{k}, 'EdgeColor', ce0, 'LineWidth', 0.75);
        end
        
        if ip.Results.ShowGaussians
            [mu_g,sigma_g,xg,g] = fitGaussianModeToHist(xi, niFirstN{k}/dxi);
            %[mu_g,sigma_g,xg,g] = fitGaussianModeToCDF(maxAcohortFirstN{k});
            plot(xg, g*dxi, 'Color', cmap(k,:), 'LineWidth', 1);
            plot(norminv(0.99, mu_g, sigma_g)*[1 1], [0 0.15], '--', 'Color', cmap(k,:), 'LineWidth', 1);
        end
        
        %plot(repmat(pct, [2 1]), repmat([0 3/5*ya(end)]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
        axis([xa(1) xa(end) ya(1) ya(end)]);
        
        % cohort label
        if lb(k)==ub(k)
            text(xa(end), ya(end), ['' num2str(lb(k)) ' s'], 'BackgroundColor', [1 1 1],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:});
        else
            text(xa(end), ya(end), ['[' num2str(lb(k)) '...' num2str(ub(k)) '] s'],...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:},...
                'BackgroundColor', [1 1 1]);
        end
        
    else % CDF plot
        [f_ecdf, t_ecdf] = ecdf(maxIntDistCat(k).maxA);
        %plot(t_ecdf, f_ecdf*aVec(k)+1-aVec(k), 'k', 'LineWidth', 1.5);
        plot(t_ecdf, f_ecdf, 'k', 'LineWidth', 1.5);
        ya = 0:0.2:1;
        axis([xa(1) xa(end) 0 1.05]);
        plot(xVec, FVec{k}, 'r');
        
        if lb(k)==ub(k)
            text(xa(end), ya(1), ['t_L = ' num2str(lb(k))],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', fset.sfont{:});
        else
            text(xa(end), ya(1), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', fset.sfont{:});
        end
        
    end
    if k==1
        text(xa(end), ya(end)*1.35, 'Lifetime cohort',...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:});
    end
    yal = [0 arrayfun(@(i) num2str(i, '%.2f'), ya(2:end), 'UniformOutput', false)];
    set(hi, 'YTick', ya, 'YTickLabel', yal, 'XTick', xa, 'XTickLabel', []);
    
    
    % plot vertical axis label (once)
    ncol = ceil(numel(lb)/ny);
    if mod(k,ny)==0
        set(hi, 'TickLength', [0.015 0], 'XTickLabel', xa);
        if floor(k/ny)==floor(ncol/2)+1
            hx = xlabel('Max. fluo. intensity (A.U.)', fset.lfont{:});
            if mod(ncol,2)==0
                xpos = get(hx, 'Position');
                xpos(1) = xa(1);
                xpos(2) = 1.2*xpos(2);
                set(hx, 'Position', xpos);
            end
        end
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==ceil(ny/2)
        hy = ylabel('Frequency', fset.lfont{:});
        ypos = get(hy, 'Position');
        if mod(ny,2)==0
            ypos(2) = 0;
        end
        %ypos(1) = 1.1*ypos(1);
        set(hy, 'Position', ypos);
    end
    
    if k>1 && ip.Results.ShowSignificance % indicate that the distributions are the same
        [hval pval] = kstest2(maxAcohortFirstN{k-1}, maxAcohortFirstN{k});
        if hval==0 
            % x: after box:
            x0 = xo + aw + 0.5;
            y0 = (ny-k)*(ah+sh) + yo + (ah+sh/2)/2;
            % line width in cm
            lw = 1/(72/2.54);
            plot(hbg, [x0 x0], y0+[0 ah+sh/2], 'Color', ce0, 'LineWidth', 1);
            plot(hbg, [x0-0.1 x0], [y0 y0]+ah+sh/2-lw/2, 'Color', ce0, 'LineWidth', 1);
            plot(hbg, [x0-0.1 x0], [y0 y0]+lw/2, 'Color', ce0, 'LineWidth', 1);
            text(x0+0.1, y0+(ah+sh/2)/2-0.1, '*', fset.lfont{:}, 'Color', ce0, 'Parent', hbg, 'VerticalAlignment', 'middle')            
        end
    end
end

% extra axes for outside legend
% if ~isprint
%     hi = axes('Units', 'pixels', 'Position', [xo+aw (ny-mod(k,ny)-1)*(ah+sh)+yo aw ah]);%'Color', 'none', 'XTick', [], 'YTick', []);
%     hold on;
%     bar(0, 1, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 1);
%     bar(0, 1, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 1);
%     axis([2 3 2 3]);
%     hl = legend('Full lifetime', ['First ' num2str(ip.Results.FirstNFrames*data(1).framerate) ' s of lifetime'], 'Location', 'NorthWest');
%     set(hl, 'Box', 'off');
%     set(gca, 'Visible', 'off');
% end


% bring background axes to front
% ch = get(gcf, 'Children');
% set(gcf, 'Children', [hbg; setdiff(ch, hbg)]);

axis(hbg, 'off');


if numel(data)==1
    filename = ['_' getCellDir(data)];
else
    filename = '_pooled';
end
% print('-depsc2', '-loose', [figPath 'maxIntensityDist' filename '.eps']);

%%
return
%========================================================
% Time until max. intensity is reached, per cohort
%========================================================
xa = 0:20:120;

timeUntilMaxInt = cell(1,nd);
for i = 1:nd
    nt = size(lftData(i).A,2);
    mask = lftData(i).A==repmat(max(lftData(i).A, [], 2), [1 nt]);
    tmp = repmat((0:nt-1)*data(1).framerate, [size(lftData(i).A,1) 1]);
    timeUntilMaxInt{i} = sum(tmp.*mask,2);
end
timeUntilMaxInt = vertcat(timeUntilMaxInt{:});


figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');

for k = 1:numel(lb)
    
    % plot grid below data
    hi = axes('Units', 'pixels', 'Position', [xo yo+(ny-mod(k-1,ny)-1)*(ah+sh) aw ah],...
        'XLim', [xa(1) xa(end)], 'XTick', xa, 'XTickLabel', [], 'YLim', [ya(1) ya(end)], 'YTickLabel', [],...
        'TickLength', [0 0], 'Color', 'none');
    set(hi, 'XGrid', 'on', 'GridLineStyle', ':');
    
    % data axes
    hi = axes('Units', 'pixels', 'Position', [xo yo+(ny-mod(k-1,ny)-1)*(ah+sh) aw ah],...
        'Color', 'none');
    hold on;
    box off;
    
    cidx = lb(k)<=lifetime_s & lifetime_s<ub(k);
    ctime = timeUntilMaxInt(cidx);
    ti = (0:200)*data(1).framerate;
    ni = hist(ctime, ti);
    ni = ni/sum(ni);
    bar(ti, ni, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 0.5);
    
    %bar(xi, niFirstN{k}, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 1);
    %bar(xi, ni{k}, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 1);
    %stairsXT(xi, niFirstN{k}, 'EdgeColor', ce0, 'LineWidth', 1);
    
    
    %plot(repmat(pct, [2 1]), repmat([0 3/5*ya(end)]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
    axis([0 120 ya(1) ya(end)]);
    
    if lb(k)==ub(k)
        text(xa(end), ya(end), ['t_L = ' num2str(lb(k)) ' s'], 'BackgroundColor', [1 1 1],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:});
    else
        text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) '] s'],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.sfont{:},...
            'BackgroundColor', [1 1 1]);
    end
    if k==1
        text(xa(end), ya(end)*1.35, 'Lifetime cohort',...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.lfont{:});
    end
    
    set(hi, 'YTick', ya, 'LineWidth', 1.5, fset.sfont{:}, 'Layer', 'top', 'TickDir', 'out',...
        'TickLength', [0.015 0], 'XTick', xa, 'XTickLabel', []);
    
    
    ncol = ceil(numel(lb)/ny);
    if mod(k,ny)==0
        set(hi, 'TickLength', [0.015 0], 'XTickLabel', xa);
        if floor(k/ny)==floor(ncol/2)+1
            hx = xlabel('Time (s)', fset.lfont{:});
            if mod(ncol,2)==0
                xpos = get(hx, 'Position');
                xpos(1) = xa(1);
                xpos(2) = 1.2*xpos(2);
                set(hx, 'Position', xpos);
            end
        end
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==ceil(ny/2)
        hy = ylabel('Frequency', fset.lfont{:});
        if mod(ny,2)==0
            ypos = get(hy, 'Position');
            ypos(2) = 0;
            ypos(1) = 1.2*ypos(1);
            set(hy, 'Position', ypos);
        end
    end
end    
    
% extra axes for outside legend
if ~isprint
    hi = axes('Units', 'pixels', 'Position', [50+floor((k)/ny)*(170+115) (ny-mod(k,ny)-1)*100+70 aw 80]);%'Color', 'none', 'XTick', [], 'YTick', []);
    hold on;
    bar(0, 1, 'BarWidth', 1, 'FaceColor', fset.cfB, 'EdgeColor', fset.ceB, 'LineWidth', 1);
    % bar(0, 1, 'BarWidth', 1, 'FaceColor', cf0, 'EdgeColor', ce0, 'LineWidth', 1);
    axis([2 3 2 3]);
    hl = legend('Time until max. intensity', 'Location', 'NorthWest');
    set(hl, 'Box', 'off');
    set(gca, 'Visible', 'off');
end


%%

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = arrayfun(@(i) i.A(:,1:f), lftData, 'UniformOutput', false);
M = vertcat(M{:});

% Growth in the first 4 frames
T = 10;
MA = M(T<=lifetime_s,:);
MB = M(lifetime_s<T,:);

% percentiles
pRef = prctile(MA, [2.5 50 95],1);

fi = 2:4;
ti = (fi-1)*data(1).framerate;

ya = 0:20:140;

figure;
hold on;
fill([ti ti(end:-1:1)], [pRef(1,fi) pRef(3,fi(end:-1:1))], 'r', 'EdgeColor', 'none');
plot(ti, pRef(2,fi), '.-');
xlabel('Time (s)');
ylabel('Fluo. intensity (A.U.)');
axis([0 4 ya(1) ya(end)]);

p = prctile(MB, [2.5 50 95],1);

figure;
hold on;
fill([ti ti(end:-1:1)], [p(1,fi) p(3,fi(end:-1:1))], 'r', 'EdgeColor', 'none');
plot(ti, p(2,fi), '.-');
xlabel('Time (s)');
ylabel('Fluo. intensity (A.U.)');
axis([0 4 ya(1) ya(end)]);

% index of tracks with higher intensity in first X frames
idx = sum(MB(:,fi)>repmat(pRef(3,fi),[size(MB,1) 1]),2);

%%

% first 10 frames, all tracks, lifetime cohorts in 10s intervals

xa = 1:f;
cb = 0:10:100;
nc = numel(cb)-1;

cmap = jet(nc);
cv = rgb2hsv(cmap);
cv(:,2) = 0.2;
cv = hsv2rgb(cv);

figure;
hold on;
for c = 1:nc
    Mx = M(cb(c)<=lifetime_s & lifetime_s<cb(c+1),:);
    
    stdVec = nanstd(Mx,[],1);
    plot(xa,nanmean(Mx,1)+stdVec, '--', 'Color', cv(c,:));
    plot(xa,nanmean(Mx,1)-stdVec, '--', 'Color', cv(c,:));
    plot(xa,nanmean(Mx,1), 'Color', cmap(c,:));
end

%%
% first 10 frames, tracks above max. intensity threshold, lifetime cohorts in 10s intervals
T = 110;

figure;
hold on;
for c = 1:nc

    idx = cb(c)<=lifetime_s & lifetime_s<cb(c+1) & maxA'>T;
    Mx = M(idx,:);
    
    stdVec = nanstd(Mx,[],1);
    plot(xa,nanmean(Mx,1)+stdVec, '--', 'Color', cv(c,:));
    plot(xa,nanmean(Mx,1)-stdVec, '--', 'Color', cv(c,:));
    
    plot(xa,nanmean(Mx,1), 'Color', cmap(c,:));
end

%%
% Growth rate during the first 3 frames (2:4)
% M: matrix of concatenated intensity traces
fi = 2:5;
nf = numel(fi);
X = repmat(fi,[size(M,1) 1]);
Y = M(:,fi);
Sxx = nansum((X-repmat(nanmean(X,2),[1 nf])).^2,2)/(nf-1);
Sxy = nansum((X-repmat(nanmean(X,2),[1 nf])).*(Y-repmat(mean(Y,2), [1 nf])),2)/(nf-1);
D = Sxy./Sxx;
b0 = mean(Y,2)-mean(X,2).*D;
% D(isnan(D)) = [];
% figure; hist(D, 100)

% k = 144;
% figure; plot(X(k,:), Y(k,:)); hold on; plot(X(k,:), D(k)*X(k,:)+b0(k,:), 'r')

% Reference distribution of initial growth rate
figure;
hold on;
di = 2;
xi = -100:di:100;
for c = 1:nc
    cidx = cb(c)<=lifetime_s & lifetime_s<cb(c+1) & maxA'>110;
    iD = D(cidx);
    iD = iD(~isnan(iD) & iD~=0);
    ni = hist(iD, xi);
    ni = ni/sum(ni)/di;
    plot(xi, ni, '-', 'Color', cmap(c,:));
end



%%
t = [0 10 20 30 40 50];
nt = numel(t);
cmap = jet(nt);


figure;
hold on;
for c = 1:nt

    idxA = t(c)<=lifetime_s & maxA'>T;
    idxB = t(c)>lifetime_s & maxA'>T;
    plot(xa,nanmean(M(idxA,:),1), 'Color', cmap(c,:));
    plot(xa,nanmean(M(idxB,:),1), '--', 'Color', cmap(c,:));
    
end


%%
figure('Position', pos, 'Color', [1 1 1], 'PaperPositionMode', 'auto');

xa = (0:ip.Results.FirstNFrames-1)*data(1).framerate;
ya = 40:5:80;
ya = -10:5:10;

for k = 1:numel(lb)
    aw = 6*40;
    hi = axes('Units', 'pixels', 'Position', [80+floor((k-1)/ny)*(170+115) (ny-mod(k-1,ny)-1)*100+70 aw 80]);
    hold on;
    box off;
    %set(hi, 'XGrid', 'on', 'GridLineStyle', ':');
    
    % first X frames in current cohort
    Mx = M(lb(k)<=lifetime_s & lifetime_s<=ub(k),:);
    
    %plot(xa,Mx(1,:));
    
    %plot(xa,nanmean(Mx,1));
    plot(xa([1 end]), [0 0], 'k--');
    plot(xa(1:end-1),nanmedian(diff(Mx,1,2),1));

    
    
    % for each lifetime cohort, get cumulative intensity of the first 4 frames
    %M = arrayfun(@(i) vertcat(i.int4{k}{:}), res(1), 'UniformOutput', false);
    %M = vertcat(M{:});   
    
    %Tindex = arrayfun(@(i) vertcat(i.maxA{k}), res(1), 'UniformOutput', false);
    %Tindex = horzcat(Tindex{:});
    
    
    %xa = 1:4;
    %ya = 0:200:1000;
    
    %plot(M', 'k');
    
    %plot(M(Tindex>=T_int,:)', 'g-');
    %plot(M(Tindex<T_int,:)', 'r--');

    
    %plot(cumsum(T_int*ones(1,4)), 'k', 'LineWidth', 2);
   
%     %plot(repmat(pct, [2 1]), repmat([0 0.03]', [1 numel(pct)]), 'r--', 'LineWidth', 1.5);
%     axis([0.5 4.5 ya(1) ya(end)]);
%     %set(gca, 'XLim', [0.5 4.5]);
%     
%     if lb(k)==ub(k)
%         text(5, ya(end), ['t_L = ' num2str(lb(k))],...
%             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', sfont{:});
%     else
%         text(5, ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
%             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', sfont{:});
%     end
%     
%     


%     axis([xa(1) xa(end) ya(1) ya(end)]);
        
    if lb(k)==ub(k)
        text(xa(end), ya(end), ['t_L = ' num2str(lb(k))],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.tfont{:});
    else
        text(xa(end), ya(end), ['t_L \in [' num2str(lb(k)) '...' num2str(ub(k)) ']'],...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', fset.tfont{:});
    end
    
    
    set(hi, 'LineWidth', 2, fset.ifont{:}, 'Layer', 'bottom', 'TickDir', 'out',...
        'TickLength', [0.02 0], 'XTick', xa, 'XLim', [xa(1) xa(end)], 'XTickLabel', [],...
        'YLim', [ya(1) ya(end)], 'YTick', ya);
    
    
    
    if mod(k,ny)==0 || k == numel(lb)
        xlabel('Time (s)', fset.sfont{:});
        set(hi, 'TickLength', [0.02 0], 'XTickLabel', xa);
    end
    if k>ny
        set(hi, 'YTickLabel', []);
    end
    
    if k==(ny+1)/2
       ylabel('Fluo. intensity', lfont{:}); 
    end
end

