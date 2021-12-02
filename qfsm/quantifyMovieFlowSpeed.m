function [speedCell]=quantifyMovieFlowSpeed(MD)
%function [speedCell,speedPeri]=quantifyMovieFlowSpeed(MD) quantifies mean
%speeds from speed map for cell mask and periphery. For better
%quantification, running windowing package is needed.
% input:
%   MD          MovieData object. Segmentation, QFSM package, and Windowing
%               package should have been run.
% output:
%   speedCell:  cell array containing mean speeds per layer upto layer 5
%               It will have a format of speedCell(windows:layers:frame). Speed at
%               cell periphery will be fram layer 1 or 2.
% Unit is in nm/min
% Sangyoon Han Nov, 2021

%% Load QFSM package
nFrames = MD.nFrames_;
nChannels = length(MD.channels_);
% Load flow vectors
iFlow = MD.getProcessIndex('FlowAnalysisProcess');
if isempty(iFlow)
    iFlow = MD.getProcessIndex('FlowTrackingProcess');
    if isempty(iFlow)
        error('Flow tracking has to be run.')
    else
        disp('Flow analysis process has not been run, flow tracking process is used instead.')
    end
end
flowProcess = MD.getProcess(iFlow);
iChan = find(flowProcess.checkChannelOutput);
%% Get the flow
flow1 = flowProcess.loadChannelOutput(iChan,'output','Md');
% Load segmented masks
dt = MD.timeInterval_; 
res = MD.pixelSize_;
%% Get the window package
iWinPack = MD.getPackageIndex('WindowingPackage');
winPack = MD.getPackage(iWinPack);
%% Get samples from windows
wsProc = winPack.getProcess(4);
for ii=1:numel(wsProc.outFilePaths_(:,iChan))
    [~,outName{ii}]=fileparts(wsProc.outFilePaths_{ii,iChan});
end
%Choose the name that contains Speed
iOutput = find(cellfun(@(x) contains(x,'Speed'),outName));
sampleStr = wsProc.loadChannelOutput(iChan,iOutput);
%% Need to interpret them
% e.g., if sampleStr is 223x22x31 in size, final one (33) is the number of
% frames, middle one (22) is depth, and first one 223 is the number of
% windows per layer? Yes I confirm that.
% figure, imagesc(squeeze(sampleStr.avg(:,1,:)))
maxNumLayers = 5;

speedCell = sampleStr.avg(:,1:maxNumLayers,:);

end

