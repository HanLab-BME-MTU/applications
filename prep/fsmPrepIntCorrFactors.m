function factors=fsmPrepIntCorrFactors(FileNameList,imageNumber,normValues)
% fsmPrepIntCorrFactors returns an intensity correction factor for bleaching
%
% SYNOPSIS   factors=fsmPrepIntCorrFactors(FileNameList,imageNumber,normValues)
%
% INPUT      FileNameList: cell array containing all file names
%            imageNumber : number of images to be analyzed
%            LTParam     : [xmin xmax] minimum and maximum values for intensity normalization
%
% OUTPUT     factors     : [1 x imageNumber] vector of intensity correction factors
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002

% Check that n is not larger than the number of files
if imageNumber>size(FileNameList,1);
    imageNumber=size(FileNameList,1);
end

% Initializing progress bar
h = waitbar(0,'Estimating intensity correction factors...');

% Build a vector with mean image intensities
m=0;
for c1=1:imageNumber
	a=double(imread(FileNameList(c1,:)));
	a=(a-normValues(1))/(normValues(2)-normValues(1));
	m(c1)=mean(a(:));
	% Update wait bar
	waitbar(c1/imageNumber,h);
end

% Close waitbar
close(h);

% Find the max mean intensity 
maxI(1:imageNumber)=max(m);

% Calculate for every image the ratio maxI/mean intensity
factors=maxI./m;

