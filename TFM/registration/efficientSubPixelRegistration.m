function efficientSubPixelRegistration(movieDataOrProcess, varargin)
%EFFICIENTSUBPIXELREGISTRATION Method for registering MovieData frames to correct for stage drift  
%
% efficientSubPixelRegistration 
%
% SYNOPSIS Process Wrapper to execute registering of MovieData frames to correct for stage drift
%          based. Core algorithm/code by :% [1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
%          "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
%          156-158 (2008). (s)
%
% INPUT   
%   MovieData - A MovieData object describing the movie to be processed
%                    
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
%  
%   
% See also: dftregistration.m 
%
% Andrew R. Jamieson Feb. 2017


% ----------- Input ----------- %%
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData') || isa(x,'Process'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
paramsIn=ip.Results.paramsIn;


%Get the indices of any previous stage drift processes  
[movieData, thisProcess, iProc] = getOwnerAndProcess(movieDataOrProcess,'EfficientSubpixelRegistrationProcess',true);

% assert(movieData==movieData1);

%Parse input, store in parameter structure
p = parseProcessParams(thisProcess, paramsIn);

%% Backup the original vectors to backup folder
% if exist(p.OutputDirectory,'dir')
%     display('Backing up the original data')
%     ii = 1;
%     backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
%     while exist(backupFolder,'dir')
%         backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
%         ii=ii+1;
%     end
%     try
%         mkdir(backupFolder);
%     catch
%         system(['mkdir -p ' backupFolder]);
%     end
%     copyfile(p.OutputDirectory, backupFolder,'f')
% end
mkClrDir(p.OutputDirectory);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name', thisProcess.getName());
else
    wtBar = -1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
imageFileNames = movieData.getImageFileNames;
nFrames = movieData.nFrames_;

% Set up the input directories (input images)
inFilePaths = cell(3,numel(movieData.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = imDirs{j};
end
thisProcess.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(3,numel(movieData.channels_));
mkClrDir(p.OutputDirectory);
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    mkClrDir(outFilePaths{1,i});
end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  [~,refName,refExt] = fileparts(thisProcess.funParams_.referenceFramePath);
else
  refName = ['refernceFrame' num2str(p.referenceFrameNum)];
  refExt  ='.tiff';
  % [~,refName,refExt] = fileparts(thisProcess.funParams_.referenceFramePath);
  % refFrame = double(imread(p.referenceFramePath));
end

outFilePaths{2,p.ChannelIndex(p.iBeadChannel)} = [p.OutputDirectory filesep refName refExt];
outFilePaths{3,p.ChannelIndex(p.iBeadChannel)} = [p.OutputDirectory filesep 'transformationParameters.mat'];

thisProcess.setOutFilePaths(outFilePaths);


%% --------------- Stage drift correction ---------------%%% 

disp('Starting correcting stage drift [EfficientSubpixelRegistrationbyCrossCorrelation]...')

% Anonymous functions for reading input/output
outFile=@(chan,frame) [outFilePaths{1,chan} filesep imageFileNames{chan}{frame}];


ImStack = zeros([movieData.imSize_ nFrames]);
beadsChannel = movieData.channels_(p.iBeadChannel);
for j = 1:nFrames, ImStack(:,:,j) = double(beadsChannel.loadImage(j)); end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  refFrame = double(imread(p.referenceFramePath));
else
  disp(['Using frame ' num2str(p.referenceFrameNum) ' as reference']);
  refFrame = ImStack(:,:,p.referenceFrameNum);
  % refFrame = double(imread(p.referenceFramePath));
end


disp('Calculating subpixel-wise registration...')
logMsg = @(t) ['Performing sub-pixel registration on frame: ' num2str(t)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;
 % Perform sub-pixel registration
if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg(1))); end


nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
Tout = zeros(nFrames, 4);

% ---- Global alias correction decision (vote across sample frames) ----
% Running detectAliasedShift per-frame is unreliable when Q scores are close.
% Instead, sample a few frames, tally which quadrant wins most often, and apply
% that correction GLOBALLY to all frames. Stage drift is systematic, so the
% winning quadrant should be consistent.
nSample   = min(10, nFrames);
sampleIdx = round(linspace(1, nFrames, nSample));
quadVotes = zeros(1, 4);   % votes for Q1, Q2, Q3, Q4

disp('Detecting alias correction from sample frames...');
for k = 1:nSample
    fk = sampleIdx(k);
    mv = ImStack(:,:,fk);
    [dft_s, Greg_s] = dftregistration(fft2(refFrame), fft2(mv), 1);
    I_s = abs(ifft2(Greg_s));
    [~,~,inf_s] = detectAliasedShift(refFrame, I_s, dft_s.row_shift, dft_s.col_shift, ...
        size(refFrame,1), size(refFrame,2));
    quadVotes(inf_s.best_quadrant) = quadVotes(inf_s.best_quadrant) + 1;
end

[~, globalBestQuad] = max(quadVotes);

% Compute median shift from sample frames to check plausibility
med_row_s = median(arrayfun(@(k) dftregistration(fft2(refFrame), fft2(ImStack(:,:,sampleIdx(k))), 1).row_shift, 1:nSample));
med_col_s = median(arrayfun(@(k) dftregistration(fft2(refFrame), fft2(ImStack(:,:,sampleIdx(k))), 1).col_shift, 1:nSample));
[nR_g, nC_g] = size(refFrame);
rowAliasPlausible_g = abs(med_row_s) > nR_g / 4;
colAliasPlausible_g = abs(med_col_s) > nC_g / 4;

global_row_corrected = rowAliasPlausible_g && (globalBestQuad == 3 || globalBestQuad == 4);
global_col_corrected = colAliasPlausible_g && (globalBestQuad == 2 || globalBestQuad == 4);
disp(['Global best quadrant: Q' num2str(globalBestQuad) ...
      '  median shift: row=' num2str(med_row_s) ' col=' num2str(med_col_s)]);
disp(['  row_corrected=' num2str(global_row_corrected) ...
      '  col_corrected=' num2str(global_col_corrected)]);

% ---- Per-frame processing ----
for frame_num = 1:nFrames

   moving = ImStack(:,:,frame_num);

   % ---- Shift estimation (dftregistration, unchanged) ----
   [DFT_output, Greg_beads] = dftregistration(fft2(refFrame), fft2(moving), p.usfac);
   row_shift = DFT_output.row_shift;
   col_shift = DFT_output.col_shift;

   % ---- Apply global alias correction decision ----
   [nR, nC] = size(refFrame);
   I_shifted = abs(ifft2(Greg_beads));
   
   row_shift_true = row_shift;
   col_shift_true = col_shift;
   if global_col_corrected
       col_shift_true = col_shift + sign(-col_shift) * nC;
   end
   if global_row_corrected
       row_shift_true = row_shift + sign(-row_shift) * nR;
   end
   
   quadInfo.row_corrected = global_row_corrected;
   quadInfo.col_corrected = global_col_corrected;
   
   if row_shift_true ~= row_shift || col_shift_true ~= col_shift
       disp(['  [alias corrected] row: ' num2str(row_shift) '->' num2str(row_shift_true) ...
             '  col: ' num2str(col_shift) '->' num2str(col_shift_true)]);
   end

   % Store corrected shifts + original shifts + flags in DFT_output
   DFT_output.row_shift       = row_shift_true;
   DFT_output.col_shift       = col_shift_true;
   DFT_output.row_shift_orig  = row_shift;       % original dftregistration value
   DFT_output.col_shift_orig  = col_shift;
   DFT_output.row_corrected   = quadInfo.row_corrected;
   DFT_output.col_corrected   = quadInfo.col_corrected;
   DFTout(frame_num) = DFT_output;

   disp(['-Frame: ' num2str(frame_num)]);
   disp(['row_shift: ' num2str(row_shift_true) '  (dft: ' num2str(row_shift) ')']);
   disp(['col_shift: ' num2str(col_shift_true) '  (dft: ' num2str(col_shift) ')']);

   % ---- Apply: keep ifft2 result (best visual quality) ----
   % The circularly-shifted image already has correct content in the
   % score-driving quadrant. We zero out the non-score-driving regions.
   fMask = computeValidMask(nR, nC, row_shift, col_shift, ...
       quadInfo.row_corrected, quadInfo.col_corrected);
   I_beads = I_shifted .* double(fMask);
   imwrite(uint16(I_beads), outFile(p.iBeadChannel, frame_num));

   % Apply same transform + mask to each additional channel
   for i = setdiff(p.ChannelIndex, p.iBeadChannel)
      iChan  = p.ChannelIndex(i);
      imChan = double(movieData.channels_(iChan).loadImage(frame_num));
      fft_imChan = DFT_apply(fft2(imChan), DFT_output);
      I2 = abs(ifft2(fft_imChan)) .* double(fMask);
      imwrite(uint16(I2), outFile(iChan, frame_num));
   end

   if ishandle(wtBar)
       tj = toc;
       nj = frame_num;
       waitbar(nj/nFrames, wtBar, sprintf([logMsg(frame_num) timeMsg(tj*nTot/nj-tj)]));
   end
end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  refFrame = double(imread(p.referenceFramePath));
else
  disp(['Using frame ' num2str(p.referenceFrameNum) ' as reference']);
  refFrame = ImStack(:,:,p.referenceFrameNum);
  % refFrame = double(imread(p.referenceFramePath));
end

T = [DFTout.row_shift; DFTout.col_shift]';

% Combined valid mask using corrected shifts
[nRows, nCols] = size(refFrame);
validMask = true(nRows, nCols);
for frame_num = 1:nFrames
    validMask = validMask & computeValidMask(nRows, nCols, ...
        DFTout(frame_num).row_shift_orig, DFTout(frame_num).col_shift_orig, ...
        DFTout(frame_num).row_corrected,  DFTout(frame_num).col_corrected);
end

% Save reference masked to valid region
refFrame_masked = refFrame .* double(validMask);
imwrite(uint16(refFrame_masked), outFilePaths{2, p.ChannelIndex(p.iBeadChannel)});
save(outFilePaths{3, p.ChannelIndex(p.iBeadChannel)}, 'DFTout', 'T', 'validMask');
if ishandle(wtBar), close(wtBar); end
disp('Finished correcting stage drift!')

end % efficientSubPixelRegistration


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask = computeValidMask(nR, nC, row_shift_orig, col_shift_orig, row_corrected, col_corrected)
% Valid region in the circular-shift (ifft2) output.
%
% The circular shift applies: output[r,c] = moving[r - row_shift, c - col_shift] (circular)
% This creates 4 quadrants. The "main" quadrant contains content from the
% non-wrapped portion of moving; the "wrapped" quadrant contains the aliased
% portion. If alias correction was applied (row/col_corrected = true), the
% score-driving region is the WRAPPED quadrant, not the main one.
%
% row_shift_orig, col_shift_orig : dftregistration output (possibly aliased)
% row_corrected, col_corrected   : true if that axis was alias-corrected

if nargin < 5, row_corrected = false; col_corrected = false; end

r_orig = round(row_shift_orig);
c_orig = round(col_shift_orig);

% ---- Row: determine valid rows ----
if r_orig > 0
    % circshift by -r_orig rows: top rows get wrapped content from bottom
    main_r    = [1 : nR - r_orig];         % output rows with non-wrapped content
    wrapped_r = [nR - r_orig + 1 : nR];    % output rows with wrapped content
elseif r_orig < 0
    main_r    = [-r_orig + 1 : nR];        % bottom portion
    wrapped_r = [1 : -r_orig];             % top portion (wrapped)
else
    main_r    = 1:nR;
    wrapped_r = 1:nR;
end

% ---- Col: determine valid cols ----
if c_orig < 0
    % circshift by +|c_orig| cols: RIGHT side gets wrapped content from LEFT
    main_c    = [1 : nC + c_orig];         % LEFT portion = main content
    wrapped_c = [nC + c_orig + 1 : nC];   % RIGHT portion = wrapped content
elseif c_orig > 0
    main_c    = [c_orig + 1 : nC];         % RIGHT portion = main
    wrapped_c = [1 : c_orig];              % LEFT portion = wrapped
else
    main_c    = 1:nC;
    wrapped_c = 1:nC;
end

% Select correct rows/cols based on alias correction
if row_corrected
    valid_r = wrapped_r;
else
    valid_r = main_r;
end

if col_corrected
    valid_c = wrapped_c;
else
    valid_c = main_c;
end

mask = false(nR, nC);
if ~isempty(valid_r) && ~isempty(valid_c)
    mask(valid_r, valid_c) = true;
end
end


function [row_true, col_true, info] = detectAliasedShift(ref, shifted, row_sh, col_sh, nR, nC)
% Detect and correct aliased shifts from circular FFT cross-correlation.
%
% dftregistration restricts shifts to [-N/2, N/2). If the true shift exceeds
% this range, the reported value is aliased: shift_true = shift +/- N.
%
% APPROACH: Score all 4 quadrants of the circularly-shifted output against the
% reference. The winning quadrant (highest corr2) directly tells us which axes
% need alias correction.
%
% Quadrant layout in the circular-shift output:
%   (output[r,c] = moving[r - row_sh, c - col_sh] circular)
%
%   For row_sh > 0: output rows [1:row_sh] = wrapped (from bottom of moving)
%                   output rows [row_sh+1:nR] = main
%   For col_sh < 0: output cols [1:nC+col_sh] = main
%                   output cols [nC+col_sh+1:nC] = wrapped (from left of moving)
%
% Winning quadrant -> correction:
%   Q1 (main  row x main  col): no correction
%   Q2 (main  row x wrap  col): col correction only
%   Q3 (wrap  row x main  col): row correction only
%   Q4 (wrap  row x wrap  col): both corrections

row_true = row_sh;
col_true = col_sh;

r = round(row_sh);
c = round(col_sh);

% Row quadrant boundaries
if r > 0
    r_main    = r+1 : nR;
    r_wrapped = 1   : r;
elseif r < 0
    r_main    = -r+1 : nR;
    r_wrapped = 1    : -r;
else
    r_main    = 1:nR;  r_wrapped = [];
end

% Col quadrant boundaries
if c < 0
    c_main    = 1      : nC+c;
    c_wrapped = nC+c+1 : nC;
elseif c > 0
    c_main    = c+1 : nC;
    c_wrapped = 1   : c;
else
    c_main    = 1:nC;  c_wrapped = [];
end

% Score all 4 quadrants
scores = [-Inf -Inf -Inf -Inf];  % [Q1 Q2 Q3 Q4]
if ~isempty(r_main) && ~isempty(c_main)
    scores(1) = quadScore(ref, shifted, r_main,    c_main);
end
if ~isempty(r_main) && ~isempty(c_wrapped)
    scores(2) = quadScore(ref, shifted, r_main,    c_wrapped);
end
if ~isempty(r_wrapped) && ~isempty(c_main)
    scores(3) = quadScore(ref, shifted, r_wrapped, c_main);
end
if ~isempty(r_wrapped) && ~isempty(c_wrapped)
    scores(4) = quadScore(ref, shifted, r_wrapped, c_wrapped);
end

[~, best] = max(scores);

% Alias correction plausibility guard:
% A true shift of ~N would leave almost no image content - physically implausible.
% Only apply correction if the wrapped region is large enough (|shift| > N/4).
% For small shifts (< N/4), the wrapped strip is tiny and corr2 there is
% unreliable - do NOT correct even if wrapped score happens to be higher.
rowAliasPlausible = abs(r) > nR / 4;
colAliasPlausible = abs(c) > nC / 4;

% Determine corrections from winning quadrant, with plausibility guard
row_corrected = rowAliasPlausible && (best == 3 || best == 4);
col_corrected = colAliasPlausible && (best == 2 || best == 4);

if col_corrected
    col_true = col_sh + sign(-col_sh) * nC;
end
if row_corrected
    row_true = row_sh + sign(-row_sh) * nR;
end

info.scores        = scores;
info.best_quadrant = best;
info.row_corrected = row_corrected;
info.col_corrected = col_corrected;
info.row_sh_orig   = row_sh;
info.col_sh_orig   = col_sh;
end


function s = quadScore(ref, shifted, rows, cols)
% corr2 score in a quadrant. Returns -Inf if region is too small or flat.
if numel(rows) < 10 || numel(cols) < 10
    s = -Inf; return;
end
R = double(ref(rows, cols));
S = double(shifted(rows, cols));
if std(R(:)) < eps || std(S(:)) < eps
    s = -Inf; return;
end
s = corr2(R, S);
end


function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
%
% Rewrote all code not authored by either Manuel Guizar or Jim Fienup
% Manuel Guizar - May 13, 2016
%
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
%
% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
%
% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
%
%
% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('usfac','var')
    usfac = 1;
end

[nr,nc]=size(buf2ft);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);

if usfac == 0
    % Simple computation of error and phase difference without registration
    CCmax = sum(buf1ft(:).*conj(buf2ft(:)));
    row_shift = 0;
    col_shift = 0;
elseif usfac == 1
    % Single pixel registration
    CC = ifft2(buf1ft.*conj(buf2ft));
    CCabs = abs(CC);
    [row_shift, col_shift] = find(CCabs == max(CCabs(:)));
    CCmax = CC(row_shift,col_shift)*nr*nc;
    % Now change shifts so that they represent relative shifts and not indices
    row_shift = Nr(row_shift);
    col_shift = Nc(col_shift);
elseif usfac > 1
    % Start with usfac == 2
    CC = ifft2(FTpad(buf1ft.*conj(buf2ft),[2*nr,2*nc]));
    CCabs = abs(CC);
    [row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');
    CCmax = CC(row_shift,col_shift)*nr*nc;
    % Now change shifts so that they represent relative shifts and not indices
    Nr2 = ifftshift(-fix(nr):ceil(nr)-1);
    Nc2 = ifftshift(-fix(nc):ceil(nc)-1);
    row_shift = Nr2(row_shift)/2;
    col_shift = Nc2(col_shift)/2;
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac));
        % Locate maximum and map back to original pixel grid 
        CCabs = abs(CC);
        [rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');
        CCmax = CC(rloc,cloc);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    
    end

    % If its only one row or column the shift along that dimension has no
    % effect. Set to zero.
    if nr == 1,
        row_shift = 0;
    end
    if nc == 1,
        col_shift = 0;
    end
    
end  

rg00 = sum(abs(buf1ft(:)).^2);
rf00 = sum(abs(buf2ft(:)).^2);
error = 1.0 - abs(CCmax).^2/(rg00*rf00);
error = sqrt(abs(error));
diffphase = angle(CCmax);

% output=[error,diffphase,row_shift,col_shift, Nr, Nc, nr, nc];

output.error = error;
output.diffphase = diffphase;
output.row_shift = row_shift;
output.col_shift = col_shift;
output.Nr = Nr;
output.Nc = Nc;
output.nr = nr;
output.nc = nc;
output.usfac = usfac;

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(1i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(1i*diffphase);
end
end


function Greg = DFT_apply(buf2ft, p)
     
   % Compute registered version of buf2ft
   if p.usfac > 0
       [p.Nc, p.Nr] = meshgrid(p.Nc, p.Nr);
       Greg = buf2ft.*exp(1i*2*pi*(-p.row_shift*p.Nr/p.nr-p.col_shift*p.Nc/p.nc));
       Greg = Greg*exp(1i*p.diffphase);
   elseif p.usfac == 0
       Greg = buf2ft*exp(1i*p.diffphase);
   end
end


function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff', 'var')~=1, roff=0;  end
if exist('coff', 'var')~=1, coff=0;  end
if exist('usfac','var')~=1, usfac=1; end
if exist('noc',  'var')~=1, noc=nc;  end
if exist('nor',  'var')~=1, nor=nr;  end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
kernr=exp((-1i*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
end


function [ imFTout ] = FTpad(imFT,outsize)
% imFTout = FTpad(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1]
% outsize   - Output size of array [ny nx] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1]
% Manuel Guizar - 2014.06.02

if ~ismatrix(imFT)
    error('Maximum number of array dimensions is 2')
end
Nout = outsize;
Nin = size(imFT);
imFT = fftshift(imFT);
center = floor(size(imFT)/2)+1;

imFTout = zeros(outsize);
centerout = floor(size(imFTout)/2)+1;

% imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
%     = imFT;
cenout_cen = centerout - center;
imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2))) ...
    = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)));

imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)/(Nin(1)*Nin(2));
end