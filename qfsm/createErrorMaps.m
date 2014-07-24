function [img3C_map img3C_SNR] =createErrorMaps(stack,E,S)
% CREATEERRORMAPS creates error maps from the error vectors
%
% createErrorMaps goes through the whole M (or Md) stack of s matrices (each
%  matrix corresponds to the    matches returned by the tracker for frames
%  2 consecutive frames) and creates t<=s speed maps  each of which is the
%  average over n frames [j-n/2:j+n/2] around frame j, j in
%  [fix(n/2)+1:s-fix(n/2)]
%
% SYNOPSIS      [img3C_map img3C_SNR]=createErrorMaps(E,S,,imSize)
%
% INPUT
%            E         : output from analyzeFlow
%            S         : output from analyzeFlow
%            imgSize   : pixel size in the image domain (pixels)
%            mask      : [ 0 | 1 ] overlays vector field to speed map
%
% OUTPUT
%            img3C_map  : a cell array of matrices of size imgSize
%            img3C_SNR  : a cell array of matrices of size imgSize
%
%
% Aaron Ponti, May 15th, 2003
% Sebastien Besson, June 2011
% Adapted from fsmVectorAnalysis

% Calculate average vector fields
nMaps=size(stack,3);
img3C_map=cell(1,nMaps);
img3C_SNR=cell(1,nMaps);

imSize = size(stack(:,:,1));

for i=1:nMaps
    % Create empty maps
    match=zeros(imSize);
    SNRmap=zeros(imSize);
    
    % Put errors
    ind=sub2ind(imSize,E{i}(:,1),E{i}(:,2));
    match(ind)=31./E{i}(:,3);
    
    ind=sub2ind(imSize,S{i}(:,1),S{i}(:,2));
    SNRmap(ind)=31.*S{i}(:,3);
    
    % Calculate the approximate mean error distance and sigma for the gaussian kernel
    sigma=sqrt(numel(stack(:,:,i))/size(E{i},1))/3; % Doesn't take into accout whether the error
    % distribution occupies the whole image surface or not

    % Filter result with a gaussian kernel and the calculated sigma
    match=filterGauss2D(match,sigma);
    SNRmap=filterGauss2D(SNRmap,sigma);
    
    % Apply color map
    img3C_map{i}=applyColorMap(stack(:,:,i),match,[0 30],jet(31),16);
    img3C_SNR{i}=applyColorMap(stack(:,:,i),SNRmap,[0 30],jet(31),16);

end

