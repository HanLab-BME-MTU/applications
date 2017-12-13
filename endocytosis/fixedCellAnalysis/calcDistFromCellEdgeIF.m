%[out, res] = calcDistFromCellEdgeIF(data, varargin) calculates and plots
% distance and amplitude distributions from the cell edge for spot signals in
% all channels selected in the input
%
% Input:
%    data : structure with fields   
%          .channels : cell array of paths to TIFF files
%          .results  : path to results file
%
% Options:
%   'Channels' : channel indexes to analyze. Default: [2 3]
%   'BinSize'  : histogram bin size, in pixels. Default: 5
%
% Output:
%     out : structure with fields
%          .dist:           distance to cell edge
%          .A:              amplitude
%          .dfeHists:       distance from edge histogram      
%          .ampHists:       amplitude histogram (as a function of distance)
%          .dfeHistMean_bc: average of bias-corrected dfeHists
%          .dfeHistSD_bc:   s.d. of bias-corrected dfeHists
%                           These are the values plotted
%          .binc:           histogram bin center coordinates
%          .ampHistMean_bc: average of bias-corrected ampHists  
%          .ampHistSD_bc:   s.d. of bias-corrected ampHists
%
%     Note: 'dist' and 'A' are cell arrays of values for each channel
%
%
%     res : structure saved by processFramesIF(), with additional field
%          .ex : coordinates (Nx2) of the cell edge

% Francois Aguet, 01/2014

function [out, res] = calcDistFromCellEdgeIF(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParameter('BinSize', 5, @isposint);
ip.addParameter('Channels', [2 3]);
ip.addParameter('Display', true, @islogical);
ip.addParameter('Axis', []);
ip.addParameter('Name', '');
ip.addParameter('Hues', []);
ip.addParameter('Names', []);
ip.addParameter('Normalized', true, @islogical);
ip.addParameter('DisplayMode', 'screen');
ip.addParameter('PixelSize', 0.065, @isscalar);
ip.parse(data, varargin{:});
chIdx = ip.Results.Channels;
nc = numel(chIdx);

nd = numel(data);

hues = ip.Results.Hues;
if isempty(hues)
    hues = [0.0 0.3];
end


% load data sets
for i = 1:nd
    res(i) = load(data(i).results);
end

[ny,nx] = arrayfun(@(i) size(i.mask), res);
xv = 0:ip.Results.BinSize:max([nx ny]/2);
binc = ip.Results.BinSize/2:ip.Results.BinSize:xv(end);

se = strel('disk', 1);
for i = 1:nd
    [ny, nx] = size(res(i).mask);
    
    % cell boundary (including image boundary)
    B = bwboundaries(res(i).mask);
    B = vertcat(B{:}); % [x y] coordinates
    bmask = false(ny,nx);
    bmask(sub2ind([ny nx], B(:,1), B(:,2))) = true;
    
    % use only detected cell edge (as opposed to image boundary)
    % for distance calculation
    % emask: cell edge
    % bmask: cell edge in image boundary
    borderIdx = [1:ny 2*ny:ny:(nx-1)*ny nx*ny:-1:(nx-1)*ny+1 (nx-2)*ny+1:-ny:ny+1];
    borderMask = false(ny,nx);
    borderMask(borderIdx) = true;
    emask = bmask;
    emask(borderIdx) = false;
    emask = bwmorph(emask, 'skel');
    emask = emask | (imdilate(emask,se) & borderMask); % add edge pixel(s) at border
    bmask = bmask & ~emask;
    
    [by,bx] = find(bmask);
    [ey,ex] = find(emask);
    res(i).ex = [ex ey];
    
    % distance transform of mask (for normalization)
    D = bwdist(~res(i).mask);
    D2 = bwdist(~padarray(res(i).mask, [1 1],0));
    D2 = D2(2:end-1,2:end-1);
    D(D~=D2) = 0;
    D = D-1;
    D(D<0) = 0;
    
    % relative probability of observing specific distances from edge
    tmp = D(D~=0);
    dfeNorm = histc(tmp(tmp<=xv(end)), xv);
    dfeNorm = dfeNorm(1:end-1);
    
    if ip.Results.Normalized
        dfeNorm = dfeNorm/sum(dfeNorm);
    else
        dfeNorm = dfeNorm/max(dfeNorm);
    end
    
    
    for c = 1:nc
        % query cell boundary (free edge + image boundary)
        [idx, dist] = KDTreeBallQuery([ex ey; bx by], [res(i).ps(chIdx(c)).x' res(i).ps(chIdx(c)).y'], max(nx,ny)/2);
        
        % index of points with match;
        % results are sorted by increasing distance -> retain first point only
        edgeMatchIdx = ~cellfun(@isempty, idx);
        idx = cellfun(@(i) i(1), idx(edgeMatchIdx));
        dist = cellfun(@(i) i(1), dist(edgeMatchIdx));
        
        % retain only points closer to cell edge than image boundary
        dist = dist(idx<=numel(ex));
        
        out.dist{chIdx(c),i} = dist; % distance to cell edge
        A = res(i).ps(chIdx(c)).A(idx<=numel(ex))';
        out.A{chIdx(c),i} = A; % amplitude of retained points
        
        [dfeHist, hidx] = histc(dist, xv);
        dfeHist = dfeHist(1:end-1)./ dfeNorm;
        
        dfeHist(isnan(dfeHist)) = 0;
        if ip.Results.Normalized
            dfeHist = dfeHist/sum(dfeHist);
        end
        out.dfeHists{chIdx(c)}(i,:) = dfeHist;
        
        % for each bin, calculate average amplitude
        ahist = zeros(size(dfeHist));
        tmp = arrayfun(@(i) mean(A(hidx==i)), unique(hidx));
        ahist(unique(hidx)) = tmp;
        out.ampHists{chIdx(c)}(i,:) = ahist;
    end
end
for c = chIdx
    
    % correct for observation bias (cell size-dependent), same for all hists/channel
    w = sum(out.dfeHists{c}~=0,1);
    
    M = out.dfeHists{c}./repmat(w,[nd 1])*nd;
    mu = mean(M,1);
    s = std(M,[],1)/sqrt(nd);
    nidx = isnan(mu);
    mu(nidx) = [];
    s(nidx) = [];
    
    out.dfeHistMean_bc{c} = mu;
    out.s{c} = s;
    out.binc{c} = binc(~nidx);
    
    % amplitude
    M = out.ampHists{c}./repmat(w,[nd 1])*nd;
    mu = mean(M,1);
    s = std(M,[],1)/sqrt(nd);
    nidx = isnan(mu);
    mu(nidx) = [];
    s(nidx) = [];
    
    out.ampHistMean_bc{c} = mu;
    out.ampHistSD_bc{c} = s;
end


if ip.Results.Display
    % #objects/cell as a function of distance from cell edge
    setupFigure('DisplayMode', ip.Results.DisplayMode, 'Name', ip.Results.Name);
    hp = zeros(1,nc);
    for c = 1:nc
        ch = chIdx(c);
        xi = out.binc{ch}*ip.Results.PixelSize;
        mu = out.dfeHistMean_bc{ch};
        s = out.s{ch};
        
        fill([xi xi(end:-1:1)], [mu+s mu(end:-1:1)-s(end:-1:1)], hsv2rgb([hues(c) 0.5 1]),...
            'EdgeColor', 'none');
        hp(c) = plot(xi, mu, 'Color', hsv2rgb([hues(c) 1 0.9]), 'LineWidth', 1.5);
    end
    
    xlabel(['Distance from cell edge (' char(181) 'm)']);
    if ip.Results.Normalized
        ylabel('Frequency (normalized)');
    else
        ylabel('# objects/cell');
    end
    
    if ~isempty(ip.Results.Axis)
        axis(ip.Results.Axis{1});
    end
    if ~isempty(ip.Results.Names)
        hl = legend(hp, ip.Results.Names{:}, 'Location', 'NorthEast');
        set(hl, 'Box', 'off');
    end
    
    % Average spot amplitudes as a function of distance from cell edge
    setupFigure('DisplayMode', ip.Results.DisplayMode, 'Name', ip.Results.Name);
    for c = 1:nc
        ch = chIdx(c);
        xi = out.binc{ch}*ip.Results.PixelSize;
        mu = out.ampHistMean_bc{ch};
        s = out.ampHistSD_bc{ch};
        fill([xi xi(end:-1:1)], [mu+s mu(end:-1:1)-s(end:-1:1)], hsv2rgb([hues(c) 0.5 1]),...
            'EdgeColor', 'none');
        hp(c) = plot(xi, mu, 'Color', hsv2rgb([hues(c) 1 0.9]), 'LineWidth', 1.5);
    end
    xlabel(['Distance from cell edge (' char(181) 'm)']);
    ylabel('Average intensity (A.U.)');
    if ~isempty(ip.Results.Axis)
        axis(ip.Results.Axis{2});
    end
    if ~isempty(ip.Results.Names)
        hl = legend(hp, ip.Results.Names{:}, 'Location', 'NorthEast');
        set(hl, 'Box', 'off');
    end
end
