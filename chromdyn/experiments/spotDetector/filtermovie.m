function fMov=filtermovie(mov,par,verbose)
%FILTERMOVIE 3D gaussian filter the movie
%
% SYNOPSIS  fMov=filtermovie(mov,par,verbose)
%
% INPUT mov  : full  movie stack
%           par    : vector with gaussian sigma parameter in [sx sy sz]
%           verbose: (opt) if 1 (default), waitbar is shown. Verbose 0
%                   makes no output to the screen. If 2 or if string,
%                   progressText is shown.
%
% OUTPUT fMov : filtered movie

% c: 14/04/01

% test optional input
if nargin < 3 || isempty(verbose)
    verbose = 1;
end
if verbose == 2
    progTxt = '';
end
if ischar(verbose)
    progTxt = verbose;
    verbose = 2;
end


% init variables
fMov=zeros(size(mov));
%h=waitbar(0,'Filtering movie...');
tsteps=size(mov,5);
switch verbose
    case 1
        h= mywaitbar(0,[],tsteps,'filtering movie...');
    case 2
        progressText(0,progTxt);
end

% loop through all time points
for t=1:tsteps
    % use old border correction for reasons of consistency
    fMov(:,:,:,1,t)=fastGauss3D(mov(:,:,:,1,t),par(1:3),par(4:6),2);

    %remove  bg
    curT=fMov(:,:,:,1,t);
    fMov(:,:,:,1,t)=curT-min(curT(:));
    clear curT;
    switch verbose
        case 1
            mywaitbar(t/tsteps,h,tsteps);
        case 2
            progressText(t/tsteps);
    end

end;
if verbose == 1
    close(h);
end