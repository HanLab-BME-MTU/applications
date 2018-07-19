% Multiresolution statistical analysis for SpH signal
% in the olfactory bulb of mice
%
% Running example
%
% version 1.0
% October 2006
%
% Dimitri Van De Ville, Brice Bathellier
% Ecole Polytechnique Federale de Lausanne
%
%
% Only available for academic purposes. Please cite:
%
% B. Bathellier, D. Van De Ville, T. Blu, M. Unser, A. Carleton,
% "Wavelet-based multi-resolution statistics for optical imaging signals:
% application to automated detection of odour activated glomeruli in the
% mouse olfactory bulb," NeuroImage, submitted.
%
%
% For any other use, please contact the authors
%   Brice.Bathellier@epfl.ch
%   Dimitri.VanDeVille@epfl.ch
%   Thierry.Blu@epfl.ch
%   Michael.Unser@epfl.ch
%   Alan.Carleton@epfl.ch
%
%
% For documentation: please see ReadMe.txt
%

function MultiResWaveletStatistic(A)


%%%%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set confidence level
alpha=0.001;        % desired type I error level before Bonferroni correction for multiple testing

% Linear Models setup
CONST_T=39;        % start of stimulus response
% Signal
CONST_TAUF=29;     % FLUORESCENCE response time constant
% Transient dip
CONST_TAUx1=24;    % Time constant : first exponential
CONST_TAUx2=31.5;    % Time constant : second exponential
% Bleaching
CONST_TAUB=55;     % Time constant of bleaching 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the blank condition in the data
B=A(:,:,:,2);
% Get the stimulus condition in the data
A=A(:,:,:,1);


% Determine original size
n01=size(A,1);
n02=size(A,2);
Nt=size(A,3);

% Number of iterations of the wavelet transform
Jw=6;

% Extend data (mirror boundary conditions) to be able to perform the
% required number of iterations
n1=finddim(n01,2^Jw);
n2=finddim(n02,2^Jw);
A=extend(A,[n1 n2 Nt]);  % extend can be found in the accompagning Wavelet Processing folder


% Build the linear model with four regressors
X=zeros(Nt,4);    % empty matrix
% 1) stimulus response
Xidx=(1:Nt)';
model=zeros(Nt,1);
model(1:CONST_T-1)=zeros(CONST_T-1,1);
model(CONST_T:Nt)=1-exp(-(Xidx(CONST_T:Nt)-CONST_T)/CONST_TAUF);
X(:,1)=model;

% 2) constant
model=ones(Nt,1);
X(:,2)=model;

% 3) Bleach
model=zeros(Nt,1);
model=exp(-(Xidx-1)/CONST_TAUB)-1;
X(:,3)=model;

% 4) Dip
model=zeros(Nt,1);
model(CONST_T:Nt)=exp(-(Xidx(CONST_T:Nt)-CONST_T)/CONST_TAUx1)-exp(-(Xidx(CONST_T:Nt)-CONST_T)/CONST_TAUx2); % deep
X(:,4)=model;


% Least-square fitting (pseudo-inverse)
Xi=pinv(X);
XXi=pinv(X'*X);

% Degrees of freedom
J=Nt-rank(X);

% Contrast vectors:
% pick out constrast of stimulus response
c=[1; 0; 0; 0];
% Pick out constant
c1=[0; 1; 0; 0];
% Pick out bleach
c2=[0; 0; 1; 0];
% Pick out dip
c3=[0; 0; 0; 1];


% Type of wavelet transform
wav_type='*ortho';  % Symmetric orthogonal B-spline wavelets (Battle-Lemarie)
wav_degree=3.0;     % Degree 3 (cubic)

% Generate analysis and synthesis filters
[FA1,FS1]=FFTfractsplinefilters(n1,wav_degree,wav_type);
[FA2,FS2]=FFTfractsplinefilters(n2,wav_degree,wav_type);

% Perform spatial wavelet transform frame by frame
clear W;
fprintf('\nComputing DWT: frame %03d',0);
for iter=1:Nt,
    fprintf('\b\b\b%03d',iter);
    W(:,:,iter)=FFTwaveletanalysis2D(A(:,:,iter),FA1,FA2,Jw);
end;
fprintf('\n');


% Calculate threshold values according to spatio-wavelet statistical
% framework
a=alpha/(n01*n02);  % conservative (Bonferroni) correction for multiple testing

fprintf('Setting confidence level to %f (corrected)\nComputing thresholds... ',alpha);
[TW,TS]=threshold_search2d(Nt,J,a);
fprintf('[TW,TS]=[%f,%f]\n',TW,TS);

% Reshape the wavelet-temporal data
v=reshape(W(:,:,:),n1*n2,Nt)';

% Fit linear model (i.e find contrast values for every coefficient)
y=Xi*v;

% Pick out fluorescence increase parameter
g=c'*y;
% Pick out constant component
gc=c1'*y;
% Pick out bleach
gb=c2'*y;
% Pick out dip
gd=c3'*y;


% Residual error
e=v-X*y;
% Variance for fluorescense increase
s2=sum(e.^2,1)*(c'*XXi*c)/J;
% Variance dip component
s2d=sum(e.^2,1)*(c3'*XXi*c3)/J;


% Reshape data
g=reshape(g,n1,n2);
gc=reshape(gc,n1,n2);
gb=reshape(gb,n1,n2);
gd=reshape(gd,n1,n2);

% Reshape error and compute standard deviation (square root)
s2=reshape(s2,n1,n2); ss=sqrt(s2);
s2d=reshape(s2d,n1,n2); ssd=sqrt(s2d);

% T-values in wavelet domain
t=g./ss;        % fuorescence increase
td=gd./ssd;     % dip 

% Wavelet processing: mask contrast according to t-values that survive TW
% fluorescence
tmp=g;
tmp(find(abs(t)<TW))=0;
% dip
tmpd=gd;
tmpd(find(abs(td)<TW))=0;

% Reconstruct data (inverse DWT)
% unfiltered data
rc=FFTwaveletsynthesis2D(gc,FS1,FS2,Jw);    % Resting level 
rb=FFTwaveletsynthesis2D(gb,FS1,FS2,Jw);    % Bleach

% data thresholded in the wavelet domain
r=FFTwaveletsynthesis2D(tmp,FS1,FS2,Jw);    % Fluorescence increase
rd=FFTwaveletsynthesis2D(tmpd,FS1,FS2,Jw);  % transient dip


% Successively remove lowpass subbands for multiresolution
rhp=[];
for i1=0:Jw-2
    tmp1=tmp;
    tmp1(1:n1/2^(Jw-i1),1:n2/2^(Jw-i1))=0;  
    rhp(:,:,i1+1)=FFTwaveletsynthesis2D(tmp1,FS1,FS2,Jw);
end



% Special reconstruction of error estimate
rs=FFTwaveletsynthesis2Dabs(ss,FS1,FS2,Jw,ones(n1,n2));   % Fluorescence increase
rsd=FFTwaveletsynthesis2Dabs(ssd,FS1,FS2,Jw,ones(n1,n2)); % Dop

% Cropping to size of original data
r=r(1:n01,1:n02);       % Fluorescence increase
rc=rc(1:n01,1:n02);     % Resting (constant)
rb=rb(1:n01,1:n02);     % Bleaching
rd=rd(1:n01,1:n02);     % Transient dip
rhp=rhp(1:n01,1:n02,:); % Subbands selections
    
% Cropping error to size of original data
rs=rs(1:n01,1:n02);     % Fluorescence increase
rsd=rsd(1:n01,1:n02);


% Masks in spatial domain
% Mask for the full reconstruction of fluorescence increase
mask=ones(size(r));
mask(find((r(:,:)./rs)<TS))=0;  % positive activation

% Mask for transient dip 
maskd=ones(size(rd));
maskd(find((rd(:,:)./rsd)<TS))=0;

% Masks for subbands selections of fluorescence increase
mskhp=ones(size(mask));
maskhp=mask;
for i1=1:Jw-1
    idx=find((rhp(:,:,i1)./rs)<TS);
    maskhp(idx)=0;       % Consistency
    mskhp(:,:,i1)=maskhp;
end

% Normalize with constant component
r=r./rc;                            % Fluorescence increase
rhp=rhp./repmat(rc,[1 1 Jw-1]);     % Subbands selections
rd=rd./rc*min(X(:,4),[],1);         % Transient dip

% Masked images (including statistical final decision)
sol=r.*mask;                        % Fluorescence increase
sold=rd.*maskd;                     % Transient dip
solhp=rhp.*mskhp;                   % Subbands selections




%%%%%%%%%%%%%%%%%%%  Graphical Output  %%%%%%%%%%%%%%%%%%%%%

% Display results for multiresolution analysis
figure;
% Full reconstruction
subplot(2,3,1);
imagesc(sol); colorbar; axis image; title('Full reconstruction for activation');
% Low pass removed
subplot(2,3,2);
imagesc(solhp(:,:,1)); colorbar; axis image; title('Maximum width: 544 \mu m');
% First band pass removed 
subplot(2,3,3);
imagesc(solhp(:,:,2)); colorbar; axis image; title('Maximum width: 271 \mu m');
% Second band pass removed 
subplot(2,3,4);
imagesc(solhp(:,:,3)); colorbar; axis image; title('Maximum width: 135 \mu m');
% Third band pass removed 
subplot(2,3,5);
imagesc(solhp(:,:,4)); colorbar; axis image; title('Maximum width: 66 \mu m');
% Fourth band pass removed 
subplot(2,3,6);
imagesc(solhp(:,:,5));  colorbar; axis image; title('Maximum width: 30 \mu m');


% Display different signal component
figure;
% Full reconstruction for activation
subplot(2,2,1);
imagesc(sol); colorbar; axis image; title('Full reconstruction');
% Map for activation
subplot(2,2,2);
imagesc(solhp(:,:,4)); colorbar; axis image; title('Map for activation');
% Map of bleaching amplitudes
subplot(2,2,3);
imagesc(rb); colorbar; axis image; title('Map of bleaching amplitude (no stat)');
% Map of transient deep
subplot(2,2,4);
imagesc(sold); colorbar; axis image; title('Map of transient deep (no stat)');
