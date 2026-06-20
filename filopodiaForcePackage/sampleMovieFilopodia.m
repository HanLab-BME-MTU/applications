function sampleMovieFilopodia(movieData)
%SAMPLEMOVIEFILOPODIA  Process 5. Sample traction force + talin intensity at
%filopodium tip / shaft / base, and along the shaft line.
%
% For each assembled filopodium (P4), at each frame:
%   - traction is read from the TFMPackage ForceFieldCalculationProcess
%     (forceField(t).pos / .vec, traction stress in Pa) by interpolation,
%   - talin-GFP intensity is read from the talin channel image,
% sampled at the tip, at each shaft adhesion, at the base, and at regular
% arc-length steps along the tip->base line. Traction is decomposed into an
% axial component (along the shaft, + = outward toward the tip) and a lateral
% component. Output is per-filopodium time series + shaft profiles.
% Sangyoon J. Han / 2026

%% process & params
iProc = movieData.getProcessIndex('FilopodiaSamplingProcess', 1, 0);
assert(~isempty(iProc), 'No FilopodiaSamplingProcess found.');
proc = movieData.processes_{iProc};
p = parseProcessParams(proc);
iChan = p.ChannelIndex;

% upstream: classification (P4)
iCls = getfielddef(p,'ClassProcessIndex',[]);
if isempty(iCls), iCls = movieData.getProcessIndex('FilopodiaClassificationProcess',1,0); end
assert(~isempty(iCls), 'Run FilopodiaClassificationProcess (P4) first.');
clsProc = movieData.processes_{iCls};
clsFile = clsProc.outFilePaths_{1,iChan};
if isempty(clsFile)||exist(clsFile,'file')~=2, clsFile = fullfile(clsProc.funParams_.OutputDirectory,'filoClassification.mat'); end
Scls = load(clsFile, 'filopodia', 'roleByTrack');
filopodia = Scls.filopodia;

% traction source: TFMPackage ForceFieldCalculationProcess
forceField = loadForceField(movieData, p);

%% I/O
inFilePaths = cell(1, numel(movieData.channels_)); inFilePaths{1,iChan} = clsFile;
proc.setInFilePaths(inFilePaths);
outDir = p.OutputDirectory; mkClrDir(outDir);
outFile = fullfile(outDir, 'filoSamples.mat');
outFilePaths = cell(1, numel(movieData.channels_)); outFilePaths{1,iChan} = outFile;
proc.setOutFilePaths(outFilePaths);

pix = movieData.pixelSize_; if isempty(pix), pix = 1; end
shaftStep = getfielddef(p,'ShaftSampleStep',3);   % px along shaft
sampR     = getfielddef(p,'SampleRadius',1);      % px; local averaging radius for talin
nNorm     = getfielddef(p,'NormProfileN',50);     % # points for length-normalized profile
sN        = linspace(0,1,nNorm)';                 % common normalized arc-length grid (0=tip,1=base)

%% build per-frame traction interpolants (lazy cache)
nF = movieData.nFrames_;
Fx = cell(1,nF); Fy = cell(1,nF); haveF = false(1,nF);
for t = 1:nF
    if ~isempty(forceField) && numel(forceField)>=t && ~isempty(forceField(t).pos)
        P = forceField(t).pos; V = forceField(t).vec;
        ok = all(isfinite(P),2) & all(isfinite(V),2);
        if nnz(ok) >= 3
            Fx{t} = scatteredInterpolant(P(ok,1),P(ok,2),V(ok,1),'linear','none');
            Fy{t} = scatteredInterpolant(P(ok,1),P(ok,2),V(ok,2),'linear','none');
            haveF(t) = true;
        end
    end
end

%% sample each filopodium
nFil = numel(filopodia);
filoSamples = struct('tipTrackId',{},'frames',{}, ...
    'tipForce',{},'tipForceAxial',{},'tipForceLateral',{},'tipTalin',{}, ...
    'baseForce',{},'baseTalin',{},'shaftProfile',{}, ...
    'normS',{},'normForce',{},'normForceAxial',{},'normTalin',{}, ...
    'normForceMean',{},'normForceAxialMean',{},'normTalinMean',{});

progressText(0,'Filopodia sampling','Filopodia sampling');
for f = 1:nFil
    fil = filopodia(f);
    nfr = numel(fil.frames);
    tipF=nan(1,nfr); tipFa=nan(1,nfr); tipFl=nan(1,nfr); tipI=nan(1,nfr);
    baseF=nan(1,nfr); baseI=nan(1,nfr);
    shaftProfile = cell(1,nfr);
    % length-normalized shaft profiles, one column per frame, resampled to sN grid
    nF_norm  = nan(nNorm, nfr);   % traction magnitude (Pa)
    nFa_norm = nan(nNorm, nfr);   % axial traction (+ = outward toward tip)
    nI_norm  = nan(nNorm, nfr);   % talin intensity
    for j = 1:nfr
        t = fil.frames(j);
        img = double(movieData.channels_(iChan).loadImage(t));
        tip = fil.tipPos(j,:); base = fil.basePos(j,:);
        axis_ = base - tip; aL = hypot(axis_(1),axis_(2));
        if aL < 1, continue; end
        u = axis_/aL;                       % tip->base unit; outward = -u
        % tip
        [tipF(j), tipFa(j), tipFl(j)] = sampleForce(tip, -u, Fx, Fy, haveF, t);
        tipI(j) = sampleInt(img, tip, sampR);
        % base
        [baseF(j), ~, ~] = sampleForce(base, -u, Fx, Fy, haveF, t);
        baseI(j) = sampleInt(img, base, sampR);
        % shaft profile along tip->base
        ss = (0:shaftStep:aL)';
        sp = struct('s_nm',[],'force',[],'forceAxial',[],'talin',[]);
        sp.s_nm = ss * pix;
        sp.force=nan(numel(ss),1); sp.forceAxial=nan(numel(ss),1); sp.talin=nan(numel(ss),1);
        for q = 1:numel(ss)
            pt = tip + ss(q)*u;
            [sp.force(q), sp.forceAxial(q), ~] = sampleForce(pt, -u, Fx, Fy, haveF, t);
            sp.talin(q) = sampleInt(img, pt, sampR);
        end
        shaftProfile{j} = sp;
        % length-normalized resampling (tip=0 -> base=1) onto common grid sN
        if numel(ss) >= 2
            sNorm_j = ss / aL;            % 0..1 along this frame's shaft
            [su, iu] = unique(sNorm_j(:));
            if numel(su) >= 2
                nF_norm(:,j)  = interp1(su, sp.force(iu),      sN, 'linear', NaN);
                nFa_norm(:,j) = interp1(su, sp.forceAxial(iu), sN, 'linear', NaN);
                nI_norm(:,j)  = interp1(su, sp.talin(iu),      sN, 'linear', NaN);
            end
        end
    end
    filoSamples(f).tipTrackId      = fil.tipTrackId;
    filoSamples(f).frames          = fil.frames;
    filoSamples(f).tipForce        = tipF;       % traction magnitude (Pa)
    filoSamples(f).tipForceAxial   = tipFa;      % + = outward (toward tip)
    filoSamples(f).tipForceLateral = tipFl;
    filoSamples(f).tipTalin        = tipI;
    filoSamples(f).baseForce       = baseF;
    filoSamples(f).baseTalin       = baseI;
    filoSamples(f).shaftProfile    = shaftProfile;
    % length-normalized shaft profiles (nNorm x nFrames) + per-filopodium time-average
    filoSamples(f).normS              = sN;
    filoSamples(f).normForce          = nF_norm;
    filoSamples(f).normForceAxial     = nFa_norm;
    filoSamples(f).normTalin          = nI_norm;
    filoSamples(f).normForceMean      = mean(nF_norm, 2,'omitnan');
    filoSamples(f).normForceAxialMean = mean(nFa_norm,2,'omitnan');
    filoSamples(f).normTalinMean      = mean(nI_norm, 2,'omitnan');
    progressText(f/max(1,nFil),'Filopodia sampling');
end

save(outFile, 'filoSamples', '-v7.3');
progressText(1,'Filopodia sampling');
fprintf('Sampling done: %d filopodia sampled (%d frames with traction).\n', nFil, nnz(haveF));
end

% ===================================================================
function [mag, axial, lateral] = sampleForce(pt, outDir, Fx, Fy, haveF, t)
mag=NaN; axial=NaN; lateral=NaN;
if t<1 || t>numel(haveF) || ~haveF(t), return; end
fx = Fx{t}(pt(1),pt(2)); fy = Fy{t}(pt(1),pt(2));
if isnan(fx)||isnan(fy), return; end
mag = hypot(fx,fy);
axial   = fx*outDir(1) + fy*outDir(2);          % + = outward (toward tip)
lateral = -fx*outDir(2) + fy*outDir(1);
end

% ===================================================================
function v = sampleInt(img, pt, R)
[H,W] = size(img);
x = round(pt(1)); y = round(pt(2));
if x<1||x>W||y<1||y>H, v=NaN; return; end
xs = max(1,x-R):min(W,x+R); ys = max(1,y-R):min(H,y+R);
v = mean(mean(img(ys,xs)));
end

% ===================================================================
function ff = loadForceField(movieData, p)
ff = [];
pkgName = getfielddef(p,'ForcePackageName','TFMPackage');
prcName = getfielddef(p,'ForceProcessName','ForceFieldCalculationProcess');
iPr = movieData.getProcessIndex(prcName, 1, 0);
if isempty(iPr)
    warning('No %s found; traction will be NaN. Run the TFMPackage first.', prcName);
    return;
end
fp = movieData.processes_{iPr};
% try the process file paths, then its OutputDirectory
cand = {};
try, cand{end+1} = fp.outFilePaths_{1}; end %#ok<TRYNC>
try, cand{end+1} = fullfile(fp.funParams_.OutputDirectory,'forceField.mat'); end %#ok<TRYNC>
for c = 1:numel(cand)
    fn = cand{c};
    if ~isempty(fn) && exist(fn,'file')==2
        S = load(fn);
        if isfield(S,'forceField'), ff = S.forceField; return; end
    end
end
warning('Could not load forceField from %s; traction will be NaN.', prcName);
end

% ===================================================================
function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
