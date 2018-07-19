function [] = createMaskCellAnterior(pathMD,erodeDist)
% createMaskCellAnterior creates masks that is
% only anterior part of a cell using user input)
% input:    pathMD:    path to the movieData file
%               erodeDis:   distance of erosion from the cell edge.
%                               This is usually in the case of flows from
%                               non-speclized cell where flows at the very
%                               edge is not accurate.
% output:   stores filtered flow vector mat file to analysis folder after
% backing up the original one
% flow_mag is in unit of nm/min.

if nargin <2
    erodeDist=0;
end
% Load movieData
movieDataPath = [pathMD '/movieData.mat'];
movieData = MovieData.load(movieDataPath);
nFrames = movieData.nFrames_;
nChannels = length(movieData.channels_);

iMask = movieData.getProcessIndex('MaskRefinementProcess');
maskProcess = movieData.getProcess(iMask);

% Backup the original vectors to backup folder
display('Backing up the original data')
[pathstr,~,~] = fileparts(maskProcess.outFilePaths_{1});
[pathstr1,analysis,~] = fileparts(pathstr);
backupFolder = fullfile(pathstr1, [analysis ' Backup']); % name]);
if ~exist(backupFolder,'dir')
    mkdir(backupFolder);
    copyfile(pathstr, backupFolder)
end

% channelName = @(x)movieData.getChannelPaths{x}(max(regexp(movieData.getChannelPaths{x},filesep))+1:end);   
% outputDir = cell(nChannels,1);
% for k = 1:nChannels    
%     %Create string for current directory
%     [pathstr,~,~] = fileparts(maskProcess.outFilePaths_{k});
%     outputDir{k} = fullfile(pathstr,channelName(k));
%     mkClrDir(outputDir{k});
% end
%Format string for zero-padding file names
% fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
% numStr = @(frame) num2str(frame,fString);
% 
% outFile=@(chan,frame) [outputDir{chan} filesep 'anteriorMask_' numStr(frame) '.tif'];
% mask = maskProcess.loadChannelOutput(1, 1);
% [nrow,ncol] = size(mask);
% roi = false(nrow,ncol);
for iChan=1:nChannels
    maskNames = maskProcess.getOutMaskFileNames(iChan);
    for jj=1:nFrames
        if iChan == 1 && jj==1
            I=double(movieData.channels_(iChan).loadImage(jj));
            dI = double(I)/max(I(:));
            % Load segmented masks
            maskf = maskProcess.loadChannelOutput(1, jj);
            maskf = bwmorph(maskf,'erode',erodeDist);
            maskl = maskProcess.loadChannelOutput(1, nFrames);
            maskl = bwmorph(maskl,'erode',erodeDist);
            % showing    
            mask_com(:,:,1) = dI+0.2*double(maskf);
            mask_com(:,:,2) = dI+0.2*double(maskl);
            mask_com(:,:,3) = dI;
            h=figure;
            imshow(mask_com,[]);
            % Get user input by impoly that will be excluded
            display('Please draw a polygon where you want to show the vector field. Finish with right click after a final point')
            hpoly = impoly;
            polyMask_first = createMask(hpoly);
            delete(hpoly);
            close(h)
        end
        mask = maskProcess.loadChannelOutput(iChan, jj);
        mask = bwmorph(mask,'erode',erodeDist);
        % subtract the user input roi from the segmented mask
        roi = ~polyMask_first & mask;
        imwrite(roi,[maskProcess.outFilePaths_{iChan} filesep maskNames{1}{jj}]);
    end
end

