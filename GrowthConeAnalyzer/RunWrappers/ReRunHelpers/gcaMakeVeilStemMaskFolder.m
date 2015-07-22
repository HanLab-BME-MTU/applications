function [ output_args ] = gcaMakeVeilStemMaskFolder( movieData,movefile)
%GCAMakeVeilStemMask
% was makeNeuriteBodyMaskFolder until 2015067
%
%
% INPUT:
% movieData
%
%
if nargin<2
    movefile=true; 
end
%% Move the old file
source = [movieData.outputDirectory_ filesep 'masks' filesep 'masks_for_channel_1'];
if movefile == true

Global = [movieData.outputDirectory_ filesep 'GlobalThresholding'];
if ~isdir(Global)
    mkdir(Global)
end
dest = [Global filesep 'masks_for_channel_1'];
copyfile(source, dest);

source1 = [movieData.outputDirectory_ filesep 'masks' filesep 'threshold_values_for_channel_1.mat'];
dest1 = [Global filesep 'threshold_values_for_channel_1.mat'];
copyfile(source1,dest1);

rmdir(source,'s');
mkdir(source);
end 
%% Load the veilStemMasks and put them in the original thereshold folder
%  until can re-adjust the input-output path workflow easier to keep them in
%  the same format as the default - input output.

fileToLoad = [movieData.outputDirectory_ filesep 'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
    'III_veilStem_reconstruction' filesep 'Channel_1' filesep 'veilStem.mat'];
load(fileToLoad);
nImTot = numel(veilStem);
fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];

for iFrame = 1:nImTot
    mask  = veilStem(iFrame).finalMask;
    
    % should have already been done but just to make sure
    % get largest CC
    mask = double(getLargestCC(mask));
    
    % set the border to zero - this just helps make sure that the windows do not
    % error
    dims = size(mask);
    mask(1:dims(1),1) =0;
    mask(1:dims(1),dims(2))=0;
    mask(1,1:dims(2))= 0;
    mask(dims(1),1:dims(2)) =0;
    
    % should have been filled but again just in case
    
    mask = imfill(mask,'holes');
    imwrite(mask,[ source filesep 'veilStemMask' num2str(iFrame,fmt) '.tif']);
   % if overlaymovie ==1
        % img = double(imread(analInfo(iFrame).imgPointer));
        img =  movieData.channels_(1).loadImage(iFrame);
        img = double(img); 
        [ny,nx] = size(img);
        setFigure(nx,ny,'off');
        imshow(-img,[]);
        hold on
        roiYX = bwboundaries(mask);
        text(5,10,'VeilStem Outline', 'Color','g');
        cellfun(@(x) plot(x(:,2),x(:,1),'color','g'),roiYX);
        if ~isdir([source filesep 'Overlays']);
            mkdir([source filesep 'Overlays']) ;
        end
        saveas(gcf,[source filesep 'Overlays' filesep 'VeilStemOverlay' num2str(iFrame,fmt) '.png']);
        close gcf
        
    %end
    clear mask
end
% if overlaymovie ==1
%     cd([saveDir filesep 'Overlays'])
%
%     execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
%     system(execute);
% end

