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
        intMat_Ia = NaN(nt,data(i).movieLength,nCh);
        startBuffer_Ia = NaN(nt,b,nCh);
        endBuffer_Ia = NaN(nt,b,nCh);
        sigma_r_Ia = NaN(nt,data(i).movieLength+2*b,nCh);
        gapMat_Ia = false(nt,data(i).movieLength);
        for k = 1:nt
            intMat_Ia(k,1:trackLengths(idx_Ia(k)),:) = tracks(k).A';
            startBuffer_Ia(k,:,:) = tracks(k).startBuffer.A';
            endBuffer_Ia(k,:,:) = tracks(k).endBuffer.A';
            sigma_r_Ia(k,1:trackLengths(idx_Ia(k))+2*b,:) = [tracks(k).startBuffer.sigma_r'; tracks(k).sigma_r'; tracks(k).endBuffer.sigma_r'];            
            gapMat_Ia(k,1:trackLengths(idx_Ia(k))) = tracks(k).gapVect';
        end
        
        lftData(i).intMat_Ia = intMat_Ia;
        lftData(i).startBuffer_Ia = startBuffer_Ia;
        lftData(i).endBuffer_Ia = endBuffer_Ia;
        lftData(i).sigma_r_Ia = sigma_r_Ia;
        lftData(i).gapMat_Ia = gapMat_Ia;

        [~,~] = mkdir([data(i).source 'Analysis']);
        iData = lftData(i);
        save(fpath, '-struct', 'iData');
    else
        lftData(i) = load(fpath);
    end
end