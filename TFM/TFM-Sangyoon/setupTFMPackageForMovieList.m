function setupTFMPackageForMovieList(movieListPath, tfmParams)
% setupTFMPackageForMovieList
% MovieList ? ?? MD? ?? TFMPackage? ??? ????? ????.
%
% movieListPath: ".../movieList.mat"
% tfmParams: struct (?? ?? ??)

movieListPath = char(movieListPath);
assert(exist(movieListPath,'file')==2, 'movieListPath not found: %s', movieListPath);

ML = MovieList.load(movieListPath);
nMovies = numel(ML.movieDataFile_);

% ???? ???
p = defaultTFMParams();
if nargin >= 2 && ~isempty(tfmParams)
    fn = fieldnames(tfmParams);
    for i=1:numel(fn), p.(fn{i}) = tfmParams.(fn{i}); end
end

for i = 1:nMovies
    MD = ML.getMovie(i);
    % ---- Ensure frame rate exists (required by TFMPackage sanityCheck) ----
    if isempty(MD.timeInterval_) || ~(MD.timeInterval_ > 0)
        MD.timeInterval_ = p.timeInterval;   % seconds per frame
    end
    if isempty(MD.pixelSize_) || ~(MD.pixelSize_ > 0)
        MD.pixelSize_ = p.pixelSize;   % seconds per frame
    end
    MD.save;

    % reference image path (? MD output ?? ? reference/ref_beads_bestZ.tif)
    refImgPath = fullfile(MD.outputDirectory_, 'reference', 'ref_beads_bestZ.tif');
    assert(exist(refImgPath,'file')==2, 'Reference image missing: %s', refImgPath);

    % --- TFMPackage ?? ---
    iPack = MD.getPackageIndex('TFMPackage');
    if isempty(iPack)
        MD.addPackage(TFMPackage(MD));
        iPack = MD.getPackageIndex('TFMPackage');
    end
    tfmPack = MD.getPackage(iPack);

    % --- Process 1: EfficientSubpixelRegistration (Stage drift correction) ---
    if isempty(tfmPack.processes_{1})
        tfmPack.createDefaultProcess(1);
    end
    proc1 = tfmPack.getProcess(1);
    fp = proc1.funParams_;
    
    % Apply shift to ALL channels
    if isfield(fp,'ChannelIndex')
        fp.ChannelIndex = 1:numel(MD.channels_);
    end
    
    % Use bead channel to compute shift
    if isfield(fp,'iBeadChannel')
        fp.iBeadChannel = p.beadChan;   % p.beadChan = 2
    elseif isfield(fp,'BeadChannelIndex')
        fp.BeadChannelIndex = p.beadChan;
    elseif isfield(fp,'refChannel')
        fp.refChannel = p.beadChan;
    end
    
    % Reference image (bestZ slice)
    if isfield(fp,'referenceFramePath')
        fp.referenceFramePath = refImgPath;
    elseif isfield(fp,'referenceFrame')
        fp.referenceFrame = refImgPath;
    end
    
    % Registration accuracy (if exists)
    if isfield(fp,'usfac')
        fp.usfac = p.usfac;
    end
    
    proc1.setPara(fp);

    % --- Process 2: Displacement field ---
    if isempty(tfmPack.processes_{2}), tfmPack.createDefaultProcess(2); end
    proc2 = tfmPack.getProcess(2);
    fp = proc2.funParams_;
    fp.referenceFramePath = refImgPath;
    if isfield(fp,'ChannelIndex'); fp.ChannelIndex = p.beadChan; end

    % (?? ??? ? tfmBatch.m?? ??? ???)
    if isfield(fp,'alpha'); fp.alpha = p.alpha; end
    if isfield(fp,'minCorLength'); fp.minCorLength = p.minCorLength; end
    if isfield(fp,'addNonLocMaxBeads'); fp.addNonLocMaxBeads = p.addNonLocMaxBeads; end
    if isfield(fp,'maxFlowSpeed'); fp.maxFlowSpeed = p.maxFlowSpeed; end
    if isfield(fp,'highRes'); fp.highRes = p.highRes; end
    if isfield(fp,'useGrid'); fp.useGrid = p.useGrid; end
    if isfield(fp,'mode'); fp.mode = p.mode; end
    if isfield(fp,'noFlowOutwardOnBorder'); fp.noFlowOutwardOnBorder = p.noFlowOutwardOnBorder; end
    if isfield(fp,'trackSuccessively'); fp.trackSuccessively = p.trackSuccessively; end

    proc2.setPara(fp);

    % --- Process 3: Displacement post-processing ---
    % if isempty(tfmPack.processes_{3}), tfmPack.createDefaultProcess(3); end
    % proc3 = tfmPack.getProcess(3);
    % fp = proc3.funParams_;
    % if isfield(fp,'outlierThreshold'); fp.outlierThreshold = p.outlierThreshold; end
    % if isfield(fp,'fillVectors'); fp.fillVectors = p.fillVectors; end
    % proc3.setPara(fp);
    if ~isempty(tfmPack.processes_{3})
        MD.deleteProcess(tfmPack.processes_{3})
        MD.save
    end

    % --- Process 4: Force reconstruction ---
    if isempty(tfmPack.processes_{4}), tfmPack.createDefaultProcess(4); end
    proc4 = tfmPack.getProcess(4);
    fp = proc4.funParams_;

    if isfield(fp,'YoungModulus'); fp.YoungModulus = p.YoungModulus; end
    if isfield(fp,'regParam'); fp.regParam = p.regParam; end
    if isfield(fp,'method'); fp.method = p.method; end
    if isfield(fp,'solMethodBEM'); fp.solMethodBEM = p.solMethodBEM; end
    if isfield(fp,'useLcurve'); fp.useLcurve = p.useLcurve; end

    % basis table (FastBEM?? ??? ? ??)
    if isfield(fp,'basisClassTblPath') && ~isempty(p.basisClassTblPath)
        fp.basisClassTblPath = p.basisClassTblPath;
    end

    proc4.setPara(fp);

    % --- Process 5: Strain energy (???) ---
    if numel(tfmPack.processes_) >= 5
        if isempty(tfmPack.processes_{5}), tfmPack.createDefaultProcess(5); end
        proc5 = tfmPack.getProcess(5);
        fp = proc5.funParams_;
        proc5.setPara(fp);
    end

    MD.save;
    fprintf('[setupTFMPackage] %d/%d done: %s\n', i, nMovies, MD.outputDirectory_);
end

end

function p = defaultTFMParams()
% ? tfmBatch.m ??? ?? + ?? ????? ???

p = struct();
p.beadChan = 2;          % ? ?? ?? bead channel
p.usfac = 20;

% Displacement params
p.alpha = 0.05;
p.minCorLength = 7;
p.addNonLocMaxBeads = true;
p.maxFlowSpeed = 20;
p.highRes = true;
p.useGrid = true;
p.mode = 'accurate';
p.noFlowOutwardOnBorder = 1;
p.trackSuccessively = false;

% Post-proc
p.outlierThreshold = 2;
p.fillVectors = false;

% Force recon
p.YoungModulus = 5000;   % Pa (??? ?? ??)
p.regParam = 1e-5;       % ??? ??
p.method = 'FTTC';    % 'FastBEM' or 'FTTC'
p.solMethodBEM = 'QR';   % ?: 'QR' or '1NormReg'
p.useLcurve = false;
p.basisClassTblPath = ''; % FastBEM?? ?? ?? ?? (??? ''?)
p.timeInterval = 120;  % seconds per frame
p.pixelSize = 325;
end