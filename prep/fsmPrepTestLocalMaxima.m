function cands = fsmPrepTestLocalMaxima(img,cands,parameters,imgO)
% fspPrepTestLocalMaxima selects statistically significant local maxima using loc max and background info from analyzeSpeckles
%
% SYNOPSIS   cands=fsmPrepTestLocalMaxima(img,cands,parameters)
%
% INPUT      img        : image matrix
%            cands      : Imax and Imin for every speckles obtained 
%                         from analyzeSpeckles (Delaunay triangulation)
%            parameters : [k sigmaD PoissonNoise I0]
%                         k : gives the confidence interval Xavg+/-k*sigma 
%                         for normally distributed data
%                         sigmaD, PoissonNoise, I0 : from the noise
%                         calibration model
%                         sigma=sqrt(sigmaD^2+PoissonNoise*(I-I0))
%                         I0 : mean intensity of the background images
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
% Sylvain Berlemont, Jan 2010

% Parameters
k=parameters(1);
sigmaD=parameters(2);
PoissonNoise=parameters(3);
I0=parameters(4);

Lmax = vertcat(cands(:).Lmax);
Bkg1 = vertcat(cands(:).Bkg1);
Bkg2 = vertcat(cands(:).Bkg2);
Bkg3 = vertcat(cands(:).Bkg3);

indLmax = sub2ind(size(img), Lmax(:, 1), Lmax(:, 2));
indBkg1 = sub2ind(size(img), Bkg1(:, 1), Bkg1(:, 2));
indBkg2 = sub2ind(size(img), Bkg2(:, 1), Bkg2(:, 2));
indBkg3 = sub2ind(size(img), Bkg3(:, 1), Bkg3(:, 2));

% Calculate difference Imax-Imin for every speckle
for c1=1:length(cands)
	
    % Read values for Imax and Imin from the image at the coordinates
    % specified by info 
    Imax = img(indLmax(c1));
    
    % imgO is the original filtered image
    Imin = mean(nonzeros([imgO(indBkg1(c1)), imgO(indBkg2(c1)), imgO(indBkg3(c1))]));
    
    [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=...
        fsmPrepTestSpeckleSignif(Imax,Imin,k,sigmaD,PoissonNoise,I0);
	
	% Complete cands
    cands(c1).ILmax=Imax;             % Local maximum intensity
    cands(c1).IBkg=Imin;              % Mean background intensity
    cands(c1).deltaI=deltaI;          % Intensity difference: ILmax-IBkg
    cands(c1).deltaICrit=k*sigmaDiff; % Critical intensity difference as calculated with the noise model
    cands(c1).sigmaLmax=sigmaMax;     % Error on local maximum intensity
    cands(c1).sigmaBkg=sigmaMin;      % Error on background intensity 
    cands(c1).status=status;          % Significance of the local maximum: 1, speckle; 0, weak local maximum
end

