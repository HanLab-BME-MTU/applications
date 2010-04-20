function cands = fsmPrepTestLocalMaxima(img,cands,parameters,imgO)
% fspPrepTestLocalMaxima selects statistically significant local maxima using loc max and background info from analyzeSpeckles
%
% SYNOPSIS   cands=fsmPrepTestLocalMaxima(img,cands,parameters)
%
% INPUT      img        : image matrix
%            cands      : Imax and Imin for every speckles obtained 
%                         from analyzeSpeckles (Delaunay triangulation)
%            parameters : [k sigmaD PoissonNoise I0]
%                               k : gives the confidence interval Xavg+/-k*sigma 
%                                   for normally distributed data
%        sigmaD, PoissonNoise, I0 : from the noise calibration model 
%                                       sigma=sqrt(sigmaD^2+PoissonNoise*(I-I0))
%                                       I0 : mean intensity of the background images
%            imgO       : the original (filtered) image 
%                              
% OUTPUT     cands      : cands (input) augmented by additional (statistical) information 
%
%                         cands.Lmax       : Local maximum position - [y x]
%                              .Bkg1       : First local minimum position - [y x]
%                              .Bkg2       : Second local minimum position - [y x]
%                              .Bkg3       : Third local minimum position - [y x]
%                              .ILmax      : Local maximum intensity
%                              .IBkg       : Mean background intensity
%                              .deltaI     : Intensity difference: ILmax-IBkg
%                              .deltaICrit : Critical intensity difference as calculated with the noise model
%                              .sigmaLmax  : Error on local maximum intensity
%                              .sigmaBkg   : Error on background intensity 
%                              .status     : Significance of the local maximum: 1, speckle; 0, weak local maximum
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002

% Parameters
k=parameters(1);
sigmaD=parameters(2);
PoissonNoise=parameters(3);
I0=parameters(4);

% Correct intensity at the four edges 
% img([1 size(img,1) (size(img,2)-1)*size(img,1)+1 size(img,1)*size(img,2)])=mean(img(:));

% Calculate standard deviation of the non-black pixels as statistical value for comparison
nonBlackPixels=img(img>0);
stdImg=std(nonBlackPixels(:));

% Pad img to avoid border effects.
maskR = 5;
aImg = padarray(img, [maskR,maskR], 'replicate');

% SB: Remove A

% Create addressing matrix
A=zeros(size(img));

% Calculate difference Imax-Imin for every speckle
for c1=1:length(cands)
	
	% Read values for Imax and Imin from the image at the coordinates specified by info
	Imax=img(cands(c1).Lmax(1),cands(c1).Lmax(2)); % img is the current image
    
    % SB: This branch should be seldom since the Delaunay triangulation should
    % not failed anymore.
    
	if cands(c1).Bkg1(1)==-1 % This means that the Deulaunay triangulation has failed
        % imgO is the original filtered image
		[Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignifAux(aImg,imgO,[cands(c1).Lmax(1) cands(c1).Lmax(2)],maskR,stdImg,A,parameters);
	else
        
        % imgO is the original filtered image
        Imin=mean([imgO(cands(c1).Bkg1(1),cands(c1).Bkg1(2)) imgO(cands(c1).Bkg2(1),cands(c1).Bkg2(2)) imgO(cands(c1).Bkg3(1),cands(c1).Bkg3(2))]);
        [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignif([cands(c1).Lmax(1) cands(c1).Lmax(2)],Imax,Imin,k,sigmaD,PoissonNoise,I0,A);
	end
	
	% Complete cands
    cands(c1).ILmax=Imax;             % Local maximum intensity
    cands(c1).IBkg=Imin;              % Mean background intensity
    cands(c1).deltaI=deltaI;          % Intensity difference: ILmax-IBkg
    cands(c1).deltaICrit=k*sigmaDiff; % Critical intensity difference as calculated with the noise model
    cands(c1).sigmaLmax=sigmaMax;     % Error on local maximum intensity
    cands(c1).sigmaBkg=sigmaMin;      % Error on background intensity 
    cands(c1).status=status;          % Significance of the local maximum: 1, speckle; 0, weak local maximum
end

