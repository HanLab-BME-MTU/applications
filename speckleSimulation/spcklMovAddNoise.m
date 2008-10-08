function spcklMovAddNoise(SNR)
% add noise to 16-bit images

% query user for tifs directory and make new one under it for noisy images
tifDir=uigetdir(pwd,'Choose tifs directory to which noise should be added');
outputDirTifsWithNoise=[tifDir filesep 'tifs_SNR' int2str(SNR)];
if ~isdir(outputDirTifsWithNoise)
    mkdir(outputDirTifsWithNoise);
end

% get number of images
listOfFiles = searchFiles('image',[],tifDir,0);
nFrames=length(listOfFiles);

for fm=1:nFrames
    % read image
    newIm=double(imread([char(listOfFiles(fm,2)) filesep char(listOfFiles(fm,1))]));
    % std of noise is approximated  by 1/SNR * std(image)
    noiseSD=(1/SNR)*std(newIm(:));
    noise=noiseSD*randn(size(newIm));
    % add noise to the image
    newIm=newIm+noise;
    % cut gray levels to min and max for 16-bit image
    newIm(newIm>2^16)=2^16;
    newIm(newIm<0)=0;
    newIm=round(newIm);
    % convert from double to 16-bit
    newIm=uint16(newIm);
    % save image
    imwrite(newIm,[outputDirTifsWithNoise filesep char(listOfFiles(fm,1))]);
end