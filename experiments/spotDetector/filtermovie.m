function fMov=filtermovie(mov,par,verbose)
%FILTERMOVIE 3D gaussian filter the movie
%
% SYNOPSIS  fMov=filtermovie(mov,par)
%
% INPUT mov  : full  movie stack
%           par    : vector with gaussian sigma parameter in [sx sy sz]
%           verbose: (opt) if 1 (default), waitbar is shown. Verbose 0 
%           makes no output to the screen.
%
% OUTPUT fMov : filtered movie

% c: 14/04/01

% test optional input
if nargin < 3 || isempty(verbose)
    verbose = 1;
end


% init variables
fMov=zeros(size(mov));
%h=waitbar(0,'Filtering movie...');
tsteps=size(mov,5);
if verbose
    h= mywaitbar(0,[],tsteps,'filtering movie...');
end

% loop through all time points
for t=1:tsteps
    fMov(:,:,:,1,t)=fastGauss3D(mov(:,:,:,1,t),par(1:3),par(4:6),1);

    %remove  bg
    curT=fMov(:,:,:,1,t);
    fMov(:,:,:,1,t)=curT-min(curT(:));
    clear curT;
    if verbose
        mywaitbar(t/tsteps,h,tsteps);
    end

end;
if verbose
    close(h);
end