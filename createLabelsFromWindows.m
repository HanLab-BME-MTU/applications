function labels = createLabelsFromWindows(winPoly, imSize, method)

% method is a string option which can be:
% 'band: 1 label per band, i.e. winPoly(iBand, :)
% 'sector': 1 label per sector, i.e winPoly(:, iSector)
% 'window': 1 label per window, i.e. winPoly(iBand,iSlice)

[nBands,nSectors] = size(winPoly);

labels = zeros(imSize);

%
% Band labeling
%

if strcmp(method,'band') || strcmp(method,'window')
    bandMasks = false(imSize(1),imSize(2), nBands);
    
    for iBand = 1:nBands
        e = cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(iBand,:).outerBorder}) & ...
            cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(iBand,:).innerBorder});
        
        firstSector = find(e, 1, 'first');
        lastSector = find(e, 1, 'last');
        
        % This is a little patch to avoid non-contigous windows in a
        % band. Should not happen anymore.
        if isempty(firstSector) || isempty(lastSector) || ~min(e(firstSector:lastSector))
            continue;
        end
        
        n = 2 * numel(firstSector:lastSector);
        
        X = cell(1, n);
        Y = cell(1, n);
        
        cnt = 1;
        
        for iSector = firstSector:lastSector
            X{cnt} = winPoly(iBand,iSector).outerBorder(1,:);
            Y{cnt} = winPoly(iBand,iSector).outerBorder(2,:);
            cnt = cnt + 1;
        end
        
        for iSector = lastSector:-1:firstSector
            X{cnt} = winPoly(iBand,iSector).innerBorder(1,end:-1:1);
            Y{cnt} = winPoly(iBand,iSector).innerBorder(2,end:-1:1);
            cnt = cnt + 1;
        end
        
        X = horzcat(X{:});
        Y = horzcat(Y{:});
        
        bandMasks(:,:,iBand) =  poly2mask(X, Y, imSize(1), imSize(2));
    end
end

%
% Sector labeling
%

if strcmp(method,'sector') || strcmp(method,'window')
    sectorMasks = false(imSize(1),imSize(2), nSectors);
    
    for iSector = 1:nSectors
        e = cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(:,iSector).outerBorder}) & ...
            cellfun(@(x) ~(isempty(x) || any(isnan(x(:)))), {winPoly(:,iSector).innerBorder});
        
        firstBand = find(e, 1, 'first');
        lastBand = find(e, 1, 'last');
        
        % This is a little patch to avoid non-contigous windows in a sector.
        % Should not happen anymore.
        if isempty(firstBand) || isempty(lastBand) || ~min(e(firstBand:lastBand))
            continue;
        end
        
        n = 2 * (numel(firstBand+1:lastBand-1) + 1);
        
        X = cell(1, n);
        Y = cell(1, n);
        
        cnt = 1;
        
        X{cnt} = winPoly(firstBand, iSector).outerBorder(1, :);
        Y{cnt} = winPoly(firstBand, iSector).outerBorder(2, :);
        
        for iBand = firstBand+1:lastBand-1
            cnt = cnt + 1;
            X{cnt} = winPoly(iBand, iSector).outerBorder(1, end);
            Y{cnt} = winPoly(iBand, iSector).outerBorder(2, end);
        end
        
        cnt = cnt + 1;
        X{cnt} = winPoly(lastBand, iSector).innerBorder(1, end:-1:1);
        Y{cnt} = winPoly(lastBand, iSector).innerBorder(2, end:-1:1);
        
        for iBand = lastBand-1:-1:firstBand+1
            cnt = cnt + 1;
            X{cnt} = winPoly(iBand, iSector).outerBorder(1, 1);
            Y{cnt} = winPoly(iBand, iSector).outerBorder(2, 1);
        end
        
        X = horzcat(X{:});
        Y = horzcat(Y{:});
        
        sectorMasks(:,:,iSector) =  poly2mask(X, Y, imSize(1), imSize(2));
    end
end

if exist('bandMasks','var') && exist('sectorMasks','var')
    
    % Window labeling: multiply both mask types
    [iBandFull, iSectorFull] = meshgrid(1:nBands,1:nSectors);
    ind = sub2ind([nBands,nSectors], iBandFull(:), iSectorFull(:));
    cnt = 1;
    for iBand = 1:nBands
        for iSector = 1:nSectors
            labels(bandMasks(:,:,iBand) & sectorMasks(:,:,iSector)) = ind(cnt);
            cnt = cnt + 1;
        end
    end
    
elseif exist('bandMasks','var')
    
    % Band labeling
    for iBand = 1:nBands
        labels(bandMasks(:,:,iBand)) = iBand;
    end
    
elseif exist('sectorMasks','var')
    
    % Sector labeling
    for iSector = 1:nSectors
        labels(sectorMasks(:,:,iSector)) = iSector;
    end
end