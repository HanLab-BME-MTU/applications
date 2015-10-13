function printMIP(MD)
% Title says it all -- PR, augmented from Meghan D. 2015

ip = inputParser;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.parse(MD);



% turn a specific warning off
warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');

if(MD.zSize_==1)
    error('This seems to be a 2D movie, printing a MIP is not smart.');
end

savePath=[MD.outputDirectory_ filesep 'MIP'];
if ~isdir(savePath) || ~isdir([savePath filesep 'XY']) || ~isdir([savePath filesep 'ZY']) || ~isdir([savePath filesep 'ZX']) || ~isdir([savePath filesep 'Three'])
    mkdir(savePath)
    mkdir([savePath filesep 'XY'])
    mkdir([savePath filesep 'ZY'])
    mkdir([savePath filesep 'ZX'])
    mkdir([savePath filesep 'Three'])
end
nameCells=MD.getChannel(1).getImageFileNames;
% fprintf('printing MIP:');
% fprintf(['\n' repmat('.',1,MD.nFrames_) '\n\n']);
ZXRatio=MD.pixelSizeZ_/MD.pixelSize_;
for frameIdx=1:MD.nFrames_
    maxXY=[];maxZY=[];maxZX=[];three=[];
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(frameIdx);  
        [cmaxXY,cmaxZY,cmaxZX,cthree]=computeMIPs(vol,ZXRatio);
        maxXY=[ maxXY cmaxXY];
        maxZY=[ maxZY cmaxZY];
        maxZX=[ maxZX cmaxZX];
        three=[ three cthree];        
    end
    % save the maximum intensity projections
    imwrite(maxXY, [savePath filesep 'XY' filesep 'XY_' nameCells{frameIdx} ], 'Compression', 'none');
    imwrite(maxZY, [savePath filesep 'ZY' filesep 'ZY_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(maxZX, [savePath filesep 'ZX' filesep 'ZX_'  nameCells{frameIdx}], 'Compression', 'none');
    imwrite(three, [savePath filesep 'Three' filesep nameCells{frameIdx}], 'Compression', 'none');
%     fprintf('\b|\n');
end


function [maxXY,maxZY,maxZX,three]=computeMIPs(vol,ZXRatio)

% set other parameters
stripeSize = 8; % the width of the stripes in the image that combines all three maximum intensity projections
stripeColor = 0; %the stripe color, a number between 0 (black) and 1 (white).  (If you're not using all of the bit depth then white might be much lower than 1, and black might be higher than 0.)

ScaledZ=ceil(size(vol,3)*ZXRatio);
% find the maximum intensity projections
maxXY = (max(vol, [], 3));
maxZY = imresize((squeeze(max(vol, [], 2))),[size(vol,1) ScaledZ]);
maxZX = imresize((squeeze(max(vol, [], 1))),[size(vol,2) ScaledZ]);

% generate a single image with all three projections
threeTop = [maxXY, stripeColor*ones(size(vol,1), stripeSize), maxZY];
threeBottom = [maxZX', stripeColor*ones(ScaledZ, ScaledZ+stripeSize)];
three = [threeTop; stripeColor*ones(stripeSize, size(vol,2)+ScaledZ+stripeSize); threeBottom];

maxXY = uint8((2^8-1)*mat2gray(maxXY));
maxZY = uint8((2^8-1)*mat2gray(maxZY));
maxZX = uint8((2^8-1)*mat2gray(maxZX));
three = uint8((2^8-1)*mat2gray(three));
