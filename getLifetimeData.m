% Francois Aguet, 05/2012

function lftData = getLifetimeData(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Overwrite', false, @islogical);
ip.parse(varargin{:});

nd = numel(data);
% lftData(1:nd) = struct('lifetime_s', [], 'trackLengths', [], 'start', [], 'catIdx', [],...
%     'intMat_Ia', [], 'startBuffer_Ia', [], 'endBuffer_Ia', []);
for i = 1:nd
    fpath = [data(i).source 'Analysis' filesep 'lifetimeData.mat'];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        
        tracks = loadTracks(data(i), 'Mask', true, 'Category', 'all', 'Cutoff_f', 0);
        nCh = size(tracks(1).A,1);
        
        % concatenate amplitudes of master channel into matrix
        trackLengths = [tracks.end]-[tracks.start]+1;
        
        lftData(i).lifetime_s = [tracks.lifetime_s];
        lftData(i).trackLengths = trackLengths;
        lftData(i).start = [tracks.start];
        lftData(i).catIdx = [tracks.catIdx];
        if isfield(tracks, 'significantSignal')
            lftData(i).significantSignal = [tracks.significantSignal];
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
        %xMat_Ia = NaN(nt,data(i).movieLength);
        %yMat_Ia = NaN(nt,data(i).movieLength);
        sbA = NaN(nt,b,nCh);
        ebA = NaN(nt,b,nCh);
        sbSigma_r = NaN(nt,b,nCh);
        ebSigma_r = NaN(nt,b,nCh);
        sigma_r = NaN(nt,nf,nCh);
        SE_sigma_r = NaN(nt,nf);
        gapMat_Ia = false(nt,nf);
        for k = 1:nt
            A(k,1:trackLengths(idx_Ia(k)),:) = tracks(k).A';
            A_pstd(k,1:trackLengths(idx_Ia(k))) = tracks(k).A_pstd(1,:);
            %xMat_Ia(k,1:trackLengths(idx_Ia(k))) = tracks(k).x(1,:)';
            %yMat_Ia(k,1:trackLengths(idx_Ia(k))) = tracks(k).y(1,:)';
            %sigma_r(k,1:trackLengths(idx_Ia(k))+2*b,:) = [tracks(k).startBuffer.sigma_r'; tracks(k).sigma_r'; tracks(k).endBuffer.sigma_r'];
            sigma_r(k,1:trackLengths(idx_Ia(k)),:) = tracks(k).sigma_r';
            SE_sigma_r(k,1:trackLengths(idx_Ia(k))) = tracks(k).SE_sigma_r(1,:);
            sbA(k,:,:) = tracks(k).startBuffer.A';
            ebA(k,:,:) = tracks(k).endBuffer.A';
            sbSigma_r(k,:,:) = tracks(k).startBuffer.sigma_r';
            ebSigma_r(k,:,:) = tracks(k).endBuffer.sigma_r';
            gapMat_Ia(k,1:trackLengths(idx_Ia(k))) = tracks(k).gapVect';
        end
        
        lftData(i).A = A;
        %lftData(i).xMat_Ia = xMat_Ia;
        %lftData(i).yMat_Ia = yMat_Ia;
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
        lftData(i) = load(fpath);
    end
end