function fsmPrepStatSpeckleSelection(img,enhTriang,noiseParam,strg,counter)
% fsmPrepStatSpeckleSelection performs a statistical selection on local maxima based on a noise model
%
% SYNOPSIS   fsmPrepStatSpeckleSelection(img,enhTriang,noiseParam,str)
%
% INPUT      img            : loaded image
%            enhTriang      : if set to 1, an enhanced Delaunay triangulation is performed
%            noiseParam     : parameters of the noise model
%            str            : string format to save the result with the correct index
% 
% OUTPUT     I              : speckle mask
%            The structure 'info' is saved to disk
%
% DEPENDENCES   fsmPrepStatSpeckleSelection uses { locMax2d.m
%                                                         locMin2d.m
%                                                         fsmPrepBkgEstimationDelaunay.m
%                                                         fsmPrepTestLocalMaxima.m }
%               fsmPrepStatSpeckleSelection is used by { fsmPrepMain.m }
%
% Aaron Ponti, October 4th, 2002

if nargin~=5
    error('5 input parameters expected');
end

% Assure class == double
img=double(img);

% Loc max maps
locMax=locmax2d(img,[5 5]);
	
% Loc min maps - [3 3]!
locMin=locmin2d(img,[3 3]);

% Use Delaunay triangulation to assign 3 local minima to each local maximum
matlabVersion=ver('MATLAB');
if str2num(matlabVersion.Version)<6.1
    % Use enhanced triangulation
    cands=fsmPrepBkgEstimationDelaunay(size(img),locMax,locMin,enhTriang); % Finds 3 loc min around each loc max
else
    % Use delaunay triangulation of Matlab version >=6.1
    cands=fsmPrepBkgEstimDelauNoEnh(size(img),locMax,locMin); % Finds 3 loc min around each loc max
end
% Test local maxima against background using using a noise model (updates locMax and cands)
[locMax,cands]=fsmPrepTestLocalMaxima(img,locMax,cands,noiseParam);

% Save speckle information (cands and locMax) to disk for future use
indxStr=sprintf(strg,counter);
eval(strcat('save cands',filesep,'cands',indxStr,'.mat cands;')); % Save speckle info
eval(strcat('save locMax',filesep,'locMax',indxStr,'.mat locMax;')); % Save loc max positions
