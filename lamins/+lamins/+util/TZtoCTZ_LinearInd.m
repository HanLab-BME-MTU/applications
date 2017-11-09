function [ ctz ] = TZtoCTZ_LinearInd( MD, tz )
%TZtoCTZ_LinearInd Expands tz to full ctz linear index
% MD is a MovieData object
% tz is a TxZ linear index

% ctz dimensions
nChannels = length(MD.channels_);
ctzDim = [nChannels MD.nFrames_ MD.zSize_];

% repeat channel for each tz element
ch = 1:nChannels;
ch = repmat(ch,length(tz),1);

tz = repmat(tz,nChannels,1);

ctz = sub2ind(ctzDim,ch(:),tz);



end

