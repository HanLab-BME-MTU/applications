% Francois Aguet, 05/2012

function lftData = getLifetimeData(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('InputName', 'ProcessedTracks.mat', @ischar);
ip.addParamValue('OutputName', 'lifetimeData.mat', @ischar);
ip.addParamValue('ReturnValidOnly', true, @islogical);
ip.addParamValue('Cutoff_f', [], @isscalar);
ip.addParamValue('ExcludeVisitors', true, @islogical);
ip.parse(varargin{:});

nd = numel(data);
lftData(1:nd) = struct('lifetime_s', [], 'trackLengths', [], 'start', [], 'catIdx', [],...
    'A', [], 'A_pstd', [], 'sigma_r', [], 'SE_sigma_r', [],...
    'sbA', [], 'ebA', [], 'sbSigma_r', [], 'ebSigma_r', [], 'gapMat_Ia', []);
fnames = fieldnames(lftData);
for i = 1:nd
    fpath = [data(i).source 'Analysis' filesep ip.Results.OutputName];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        
        tracks = loadTracks(data(i), 'Mask', true, 'Category', 'all', 'Cutoff_f', 0, 'FileName', ip.Results.InputName);
        nCh = size(tracks(1).A,1);
        
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
        A = NaN(nt,nf,nCh);
        A_pstd = NaN(nt,nf);
        sbA = NaN(nt,b,nCh);
        ebA = NaN(nt,b,nCh);
        sbSigma_r = NaN(nt,b,nCh);
        ebSigma_r = NaN(nt,b,nCh);
        sigma_r = NaN(nt,nf,nCh);
        SE_sigma_r = NaN(nt,nf);
        gapMat_Ia = false(nt,nf);
        for k = 1:nt
            range = 1:trackLengths(idx_Ia(k));
            A(k,range,:) = tracks(k).A';
            A_pstd(k,range) = tracks(k).A_pstd(1,:);
            sigma_r(k,range,:) = tracks(k).sigma_r';
            SE_sigma_r(k,range) = tracks(k).SE_sigma_r(1,:);
            sbA(k,:,:) = tracks(k).startBuffer.A';
            ebA(k,:,:) = tracks(k).endBuffer.A';
            sbSigma_r(k,:,:) = tracks(k).startBuffer.sigma_r';
            ebSigma_r(k,:,:) = tracks(k).endBuffer.sigma_r';
            gapMat_Ia(k,range) = tracks(k).gapVect';
        end
        
        lftData(i).A = A;
        lftData(i).A_pstd = A_pstd;
        lftData(i).sigma_r = sigma_r;
        lftData(i).SE_sigma_r = SE_sigma_r;

        lftData(i).sbA = sbA;
        lftData(i).ebA = ebA;
        lftData(i).sbSigma_r = sbSigma_r;
        lftData(i).ebSigma_r = ebSigma_r;
        lftData(i).gapMat_Ia = gapMat_Ia;

        [~,~] = mkdir([data(i).source 'Analysis']);
        iData = lftData(i);
        save(fpath, '-struct', 'iData');
    else
        tmp = load(fpath);
        if isfield(tmp, 'significantSignal')
            lftData(i).significantSignal = [];
        end
        lftData(i) = tmp;
    end
    if ip.Results.ReturnValidOnly
        lftData(i).lifetime_s = lftData(i).lifetime_s(lftData(i).catIdx==1);
        lftData(i).trackLengths = lftData(i).trackLengths(lftData(i).catIdx==1);
        lftData(i).start = lftData(i).start(lftData(i).catIdx==1);
        lftData(i).catIdx = lftData(i).catIdx(lftData(i).catIdx==1);
        if isfield(lftData(i), 'significantSignal')
            lftData(i).significantSignal = lftData(i).significantSignal(lftData(i).catIdx==1);
        end
        
        idx = true(numel(lftData(i).lifetime_s),1);
        if ip.Results.ExcludeVisitors
            vidx = getVisitorIndex(lftData(i));
            idx = idx & ~vidx{1};
        end
        if ~isempty(ip.Results.Cutoff_f)
            idx = idx & lftData(i).trackLengths>=ip.Results.Cutoff_f;
        end
        for f = 1:numel(fnames)
            lftData(i).(fnames{f}) = lftData(i).(fnames{f})(idx,:,:);
        end
        if isfield(lftData(i), 'significantSignal')
            lftData(i).significantSignal = lftData(i).significantSignal(idx);
        end
    end
end