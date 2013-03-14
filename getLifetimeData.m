% Francois Aguet, 05/2012

function [lftData, rmIdx] = getLifetimeData(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('LifetimeData', 'lifetimeData.mat', @ischar);
ip.addParamValue('ReturnValidOnly', true, @islogical);
ip.addParamValue('Cutoff_f', [], @isscalar);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.addParamValue('Scale', false, @islogical);
ip.addParamValue('DisplayScaling', false, @islogical);
ip.addParamValue('RemoveOutliers', false, @islogical);
ip.parse(varargin{:});

nCh = numel(data(1).channels);
nd = numel(data);
rescale = ip.Results.Scale;
if numel(rescale)==1
    rescale = repmat(rescale, [nCh 1]);
end

fnames = {'lifetime_s', 'trackLengths', 'start', 'catIdx', 'A', 'A_pstd',...
    'sigma_r', 'SE_sigma_r', 'sbA', 'ebA', 'sbSigma_r', 'ebSigma_r', 'gapMat_Ia'};
lftData(1:nd) = cell2struct(cell(size(fnames)), fnames, 2);
vnames = fnames(1:4);
mnames = fnames(5:end);

maxA = cell(nCh,nd);
for i = 1:nd
    fpath = [data(i).source 'Analysis' filesep ip.Results.LifetimeData];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        
        tracks = loadTracks(data(i), 'Mask', true, 'Category', 'all', 'Cutoff_f', 0, 'FileName', ip.Results.ProcessedTracks);
        
        % concatenate amplitudes of master channel into matrix
        trackLengths = [tracks.end]-[tracks.start]+1;
        
        lftData(i).lifetime_s = [tracks.lifetime_s]';
        lftData(i).trackLengths = trackLengths';
        lftData(i).start = [tracks.start]';
        lftData(i).catIdx = [tracks.catIdx]';
        if isfield(tracks, 'significantSignal')
            lftData(i).significantSignal = [tracks.significantSignal]';
        end
        
        % store intensities of cat. Ia tracks
        idx_Ia = find([tracks.catIdx]==1);
        tracks = tracks(idx_Ia);
        
        nt = numel(tracks);
        b = numel(tracks(1).startBuffer.t);
        
        % store intensity matrices
        nf = data(i).movieLength;
        lftData(i).A = NaN(nt,nf,nCh);
        lftData(i).A_pstd = NaN(nt,nf);
        lftData(i).sigma_r = NaN(nt,nf,nCh);
        lftData(i).SE_sigma_r = NaN(nt,nf);

        lftData(i).sbA = NaN(nt,b,nCh);
        lftData(i).ebA = NaN(nt,b,nCh);
        lftData(i).sbSigma_r = NaN(nt,b,nCh);
        lftData(i).ebSigma_r = NaN(nt,b,nCh);
        %lftData(i).RSS = NaN(nt,nf,nCh);
        lftData(i).gapMat_Ia = false(nt,nf);
        for k = 1:nt
            range = 1:trackLengths(idx_Ia(k));
            lftData(i).A(k,range,:) = tracks(k).A';
            lftData(i).A_pstd(k,range) = tracks(k).A_pstd(1,:);
            lftData(i).sigma_r(k,range,:) = tracks(k).sigma_r';
            lftData(i).SE_sigma_r(k,range) = tracks(k).SE_sigma_r(1,:);
            lftData(i).sbA(k,:,:) = tracks(k).startBuffer.A';
            lftData(i).ebA(k,:,:) = tracks(k).endBuffer.A';
            lftData(i).sbSigma_r(k,:,:) = tracks(k).startBuffer.sigma_r';
            lftData(i).ebSigma_r(k,:,:) = tracks(k).endBuffer.sigma_r';
            %lftData(i).RSS(k,range,:) = tracks(k).RSS';
            lftData(i).gapMat_Ia(k,range) = tracks(k).gapVect';
        end
    else
        tmp = load(fpath);
        if isfield(tmp, 'significantSignal')
            lftData(i).significantSignal = [];
        end
        if isfield(tmp, 'RSS')
            tmp = rmfield(tmp, 'RSS');
        end
        lftData(i) = tmp;
    end
end
% save outside of parfor
for i = 1:nd
    [~,~] = mkdir([data(i).source 'Analysis']);
    fpath = [data(i).source 'Analysis' filesep ip.Results.LifetimeData];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        iData = lftData(i); %#ok<NASGU>
        save(fpath, '-struct', 'iData');
    end
end

if isfield(lftData(1), 'significantSignal')
    vnames = [vnames 'significantSignal'];
    fnames = [vnames mnames];
end
for i = 1:nd
    
    if rescale
       for c = 1:nCh
           maxA{c,i} = nanmax(lftData(i).A(:,:,c),[],2);
       end
    end
    
    % apply frame cutoff to all fields
    if ~isempty(ip.Results.Cutoff_f)
        idx = lftData(i).trackLengths(lftData(i).catIdx==1)>=ip.Results.Cutoff_f;
        for f = 1:numel(mnames)
            lftData(i).(mnames{f}) = lftData(i).(mnames{f})(idx,:,:);
        end
        
        idx = lftData(i).trackLengths>=ip.Results.Cutoff_f;
        for f = 1:numel(vnames)
            lftData(i).(vnames{f}) = lftData(i).(vnames{f})(idx,:);
        end
    end
    
    if ~ip.Results.ReturnValidOnly
        for f = 1:numel(vnames)
            lftData(i).([vnames{f} '_all']) = lftData(i).(vnames{f});
        end
    end
    
    % remaining fields: retain category==1
    idx = lftData(i).catIdx==1;
    for f = 1:numel(vnames)
        lftData(i).(vnames{f}) = lftData(i).(vnames{f})(idx,:);
    end
    
    % remove visitors
    if ip.Results.ExcludeVisitors
        vidx = getVisitorIndex(lftData(i));
        for f = 1:numel(fnames)
            lftData(i).visitors.(fnames{f}) = lftData(i).(fnames{f})(vidx{1},:,:);
            lftData(i).(fnames{f}) = lftData(i).(fnames{f})(~vidx{1},:,:);
        end
    end
end

av = zeros(nCh,nd);
rmIdx = [];
for c = 1:nCh
    if rescale(c)
        %maxA(c,:) = arrayfun(@(i) nanmax(i.A(:,:,c),[],2), lftData, 'UniformOutput', false);
        [a, offset, refIdx] = scaleEDFs(maxA(c,:), 'Display', ip.Results.DisplayScaling);
        av(c,:) = a;
        movieLength = min([data.movieLength]);
        for i = 1:nd
            lftData(i).A = lftData(i).A(:,1:movieLength,:);
            maxA{c,i} = a(i) * maxA{c,i};
            lftData(i).A(:,:,c) = a(i) * lftData(i).A(:,:,c);
            lftData(i).sbA(:,:,c) = a(i) * lftData(i).sbA(:,:,c);
            lftData(i).ebA(:,:,c) = a(i) * lftData(i).ebA(:,:,c);
            % Standard deviations are not scaled
            %lftData(i).A_pstd
            %lftData(i).sigma_r(:,:,mCh) = a(i) * lftData(i).sigma_r(:,:,mCh);
            %lftData(i).SE_sigma_r
            %lftData(i).sbSigma_r
            %lftData(i).ebSigma_r
        end
    end
    
    if ip.Results.RemoveOutliers && nd>=5
        outlierIdx = detectEDFOutliers(maxA(c,:), offset, refIdx);
        if ~isempty(outlierIdx)
            fprintf('Outlier data sets:\n');
            for i = 1:numel(outlierIdx)
                fprintf('Index %d: %s\n', outlierIdx(i), getShortPath(data(outlierIdx(i))));
            end
            rmv = input('Remove outliers? (y/n) ', 's');
            if strcmpi(rmv, 'y') || isempty(rmv)
                rmIdx = [rmIdx outlierIdx]; %#ok<AGROW>
            end
        end
    end
end
if ~isempty(rmIdx)
    lftData(rmIdx) = [];
    av(:,rmIdx) = [];
    maxA(:,rmIdx) = [];
end

if rescale(1)
    a = mat2cell(av,nCh,ones(1,numel(lftData)));
    [lftData.a] = deal(a{:});    
    %[lftData.maxA] = deal(maxA{:});
end
