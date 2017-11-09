function [flowMagAllCh,flowMagFrame,flowMagStruct] = filterVectorsSegmentedMask(pathMD,erodeDist)
% filterVectorsSegmentedMask filters flow vectors with a mask that is
% modified from segmented mask (e.g. showing only anterior part of a cell)
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

% Load flow vectors
iFlow = movieData.getProcessIndex('FlowAnalysisProcess');
if isempty(iFlow)
    iFlow = movieData.getProcessIndex('FlowTrackingProcess');
    if isempty(iFlow)
        error('Flow tracking has to be run.')
    else
        display('Flow analysis process has not been run, flow tracking process is used instead.')
    end
end
flowProcess = movieData.getProcess(iFlow);
flow = flowProcess.loadChannelOutput(1,'output','Md');
flow_filtered = cell(size(flow));
% Load segmented masks
iMask = movieData.getProcessIndex('MaskRefinementProcess');
maskProcess = movieData.getProcess(iMask);
dt = movieData.timeInterval_; 
res = movieData.pixelSize_;

% Backup the original vectors to backup folder
display('Backing up the original data')
[pathstr,~,~] = fileparts(flowProcess.outFilePaths_{1});
[pathstr1,analysis,~] = fileparts(pathstr);
backupFolder = fullfile(pathstr1, [analysis ' Backup']); % name]);
if ~exist(backupFolder,'dir')
    mkdir(backupFolder);
    copyfile(pathstr, backupFolder)
end

channelName = @(x)movieData.getChannelPaths{x}(max(regexp(movieData.getChannelPaths{x},filesep))+1:end);   
for k = 1:nChannels    
    %Create string for current directory
    [pathstr,~,~] = fileparts(flowProcess.outFilePaths_{k});
    outputDir{k} = fullfile(pathstr,channelName(k));
%     mkClrDir(outputDir{k});
end
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);

outFile=@(chan,frame) [outputDir{chan} filesep 'flowMaps_' numStr(frame) '.mat'];
outFileMag=@(chan) [outputDir{chan} filesep 'sflowMag.mat'];
outFileMagFrame=@(chan) [outputDir{chan} filesep 'sflowMagFrame.mat'];
% mask = maskProcess.loadChannelOutput(1, 1);
% [nrow,ncol] = size(mask);
% roi = false(nrow,ncol);
flowMagAllCh = cell(nChannels,1);
flowMagFrame = cell(nChannels,nFrames);
for iChan=1:nChannels
    flow = flowProcess.loadChannelOutput(iChan,'output','Md');
    flowMagAll = [];
    for jj=1:nFrames
        if iChan == 1 && jj==1
            % Load segmented masks
            maskf = maskProcess.loadChannelOutput(1, jj);
            maskf = bwmorph(maskf,'erode',erodeDist);
            maskl = maskProcess.loadChannelOutput(1, nFrames);
            maskl = bwmorph(maskl,'erode',erodeDist);
            % showing    
            mask_com(:,:,1) = double(maskf);
            mask_com(:,:,2) = double(maskl);
            mask_com(:,:,3) = 0;
            h=figure;
            imshow(mask_com,[]);
            hold on
            mag = 100/dt;
            flowf = flow{jj};
            flowl = flow{end};

            quiver(flowf(:,2),flowf(:,1),mag*(flowf(:,4)-flowf(:,2)),mag*(flowf(:,3)-flowf(:,1)),0,'k','Linewidth',0.5);
            quiver(flowl(:,2),flowl(:,1),mag*(flowl(:,4)-flowl(:,2)),mag*(flowl(:,3)-flowl(:,1)),0,'b','Linewidth',0.5);
            if nChannels>1
                flow2 = flowProcess.loadChannelOutput(2,'output','Md');
                flow2f = flow2{jj};
                flow2l = flow2{end};
                quiver(flow2f(:,2),flow2f(:,1),mag*(flow2f(:,4)-flow2f(:,2)),mag*(flow2f(:,3)-flow2f(:,1)),0,'r','Linewidth',0.5);
                quiver(flow2l(:,2),flow2l(:,1),mag*(flow2l(:,4)-flow2l(:,2)),mag*(flow2l(:,3)-flow2l(:,1)),0,'m','Linewidth',0.5);
            end
            % Get user input by impoly that will be excluded
            display('Please draw a polygon where you want to show the vector field. Finish with right click after a final point')
            hpoly = impoly;
            polyMask = createMask(hpoly);
            delete(hpoly);
            close(h)
        end
        mask = maskProcess.loadChannelOutput(1, jj);
        mask = bwmorph(mask,'erode',erodeDist);
        % subtract the user input roi from the segmented mask
        roi = ~polyMask & mask;

        % apply the mask to the vectors
        curFlow = flow{jj};
        % Outlier filtering
        idx = true(1,numel(curFlow(:,1)))'; %index for filtering
        dispMat = [curFlow(:,2),curFlow(:,1),curFlow(:,4)-curFlow(:,2),curFlow(:,3)-curFlow(:,1)];
        outlierThreshold = 2;
        outlierIndex = detectVectorFieldOutliers(dispMat,outlierThreshold,1);
        idx(outlierIndex) = false;
        curFlow = curFlow(idx,:);
        inMaskIdx = arrayfun(@(i,j) maskVectors(i,j,roi),curFlow(:,2),curFlow(:,1));
        curFlow = curFlow(inMaskIdx,:);
        flow_filtered{jj} = curFlow;
        % save flow magnitudes:
        dispMat = [curFlow(:,2),curFlow(:,1),curFlow(:,4)-curFlow(:,2),curFlow(:,3)-curFlow(:,1)];
        flow_mag = (dispMat(:,3).^2+dispMat(:,4).^2).^0.5*res*(60/dt); % now unit is nm/min
        flowMagAll = [flowMagAll; flow_mag];
        flowMagFrame{iChan,jj} = flow_mag;
        flowMagStruct.t(jj) = (jj-1)*movieData.timeInterval_;
        flowMagStruct.speed(jj) = mean(flow_mag);
        flowMagStruct.speed_std(jj) = std(flow_mag);
    end
    
    % Save the vectors to the original anlaysis folder
        % Fill output structure for each frame and save it
    Ms = flowProcess.loadChannelOutput(iChan,'output','Ms');
    E = flowProcess.loadChannelOutput(iChan,'output','E');
    S = flowProcess.loadChannelOutput(iChan,'output','S');
    speedMap = flowProcess.loadChannelOutput(iChan,'output','speedMap');
    img3C_map = flowProcess.loadChannelOutput(iChan,'output','img3C_map');
    img3C_SNR = flowProcess.loadChannelOutput(iChan,'output','img3C_SNR');

    disp('Results will be saved under:')
    disp(flowProcess.outFilePaths_{1,iChan});
    for jj=1:nFrames
        s.Md=flow_filtered{jj};
        s.Ms=Ms{jj};
        s.E=E{jj};
        s.S=S{jj};
        s.speedMap=speedMap{jj};
        s.img3C_map=img3C_map{jj};
        s.img3C_SNR=img3C_SNR{jj};

        save(outFile(iChan,jj),'-struct','s');
    end
    flowMagAllCh{iChan} = flowMagAll;
    save(outFileMag(iChan),'flowMagAll');
    flowMagFrameChan = flowMagFrame{iChan,:};
    save(outFileMagFrame(iChan),'flowMagFrameChan');
end
