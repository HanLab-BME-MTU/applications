function [speedCell]=quantifyMovieFlowSpeed(MD)
%function [speedCell,speedPeri]=quantifyMovieFlowSpeed(MD) quantifies mean
%speeds from speed map for cell mask and periphery. For better
%quantification, running windowing package is needed.
% input:
%   MD          MovieData object. Segmentation, QFSM package, and Windowing
%               package should have been run.
% output:
%   speedCell:  cell array containing mean speeds per layer upto layer 5
%               It will have a format of speedCell{frame}{layer} containing
%               mean speeds of all windows per layer per frame. Speed at
%               cell periphery will be fram layer 1 or 2.
% Unit is in ??
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
iChan = flowProcess.checkChannelOutput;
%% Get the flow
flow1 = flowProcess.loadChannelOutput(iChan,'output','Md');
% Load segmented masks
dt = movieData.timeInterval_; 
res = movieData.pixelSize_;
%% Get the window


end

