function lftData = getLifetimeData(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Overwrite', false, @islogical);
ip.parse(varargin{:});

nd = numel(data);
lftData(1:nd) = struct('lifetime_s', [], 'trackLengths', [], 'start', [], 'catIdx', [], 'intMat_Ia', []);
for i = 1:nd
    fpath = [data(i).source 'Analysis' filesep 'lifetimeData.mat'];
    if ~(exist(fpath, 'file')==2) || ip.Results.Overwrite
        
        tracks = loadTracks(data(i), 'Mask', true, 'Category', 'all', 'Cutoff_f', 0);
        
        % concatenate amplitudes of master channel into matrix
        lifetime_s = [tracks.lifetime_s];
        trackLengths = [tracks.end]-[tracks.start]+1;
        start = [tracks.start];
        catIdx = [tracks.catIdx];
        
        % store intensities of cat. Ia tracks
        idx_Ia = find([tracks.catIdx]==1);
        tracks = tracks(idx_Ia);
        
        nt = numel(tracks);
        intMat_Ia = NaN(nt, data(i).movieLength);
        mCh = strcmp(data(i).source, data(i).channels);
        for k = 1:nt
            intMat_Ia(k,1:trackLengths(idx_Ia(k))) = tracks(k).A(mCh,:);
        end
        
        lftData(i).lifetime_s = lifetime_s;
        lftData(i).trackLengths = trackLengths;
        lftData(i).start = start;
        lftData(i).catIdx = catIdx;
        lftData(i).intMat_Ia = intMat_Ia;
       
        [~,~] = mkdir([data(i).source 'Analysis']);
        save(fpath, 'lifetime_s', 'trackLengths', 'start', 'catIdx', 'intMat_Ia');
    else
        lftData(i) = load(fpath);
    end
end