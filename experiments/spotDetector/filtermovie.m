function fMov=filtermovie(mov,par)
%FILTERMOVIE 3D gaussian filter the movie
%
% SYNOPSIS  fMov=filtermovie(mov,par)
%
% INPUT mov  : full  movie stack 
%           par    : vector with gaussian sigma parameter in [sx sy sz] 
%       
% OUTPUT fMov : filtered movie

% c: 14/04/01

% init variables
fMov=zeros(size(mov));
%h=waitbar(0,'Filtering movie...');
tsteps=size(mov,5);
h= mywaitbar(0,[],tsteps,'filtering movie...');

% loop through all time points
for t=1:tsteps
    fMov(:,:,:,1,t)=fastGauss3D(mov(:,:,:,1,t),par(1:3),par(4:6),1);
    
    %correct filter effect in border slices
    % equalize mean intensity
    % midSl = mid slice
%      midSl=fMov(:,:,ceil(size(mov,3)/2),1,t);
%      for k=1:size(mov,3)
%          curSl=fMov(:,:,k,1,t);
%          %scale intensity of curr-slice according to mid slice
%          fMov(:,:,k,1,t)=mean(midSl(:))/mean(curSl(:))*curSl;
%      end;

    %remove  bg
    curT=fMov(:,:,:,1,t);
    fMov(:,:,:,1,t)=curT-min(curT(:));
    clear curT;
    
    mywaitbar(t/tsteps,h,tsteps);
    %waitbar(t/tsteps,h);
end;
close(h);