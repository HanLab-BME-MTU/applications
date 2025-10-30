function [tMapUnshifted]=quantifyTractionWindows(MD)
%function [tMapUnshifted]=quantifyTractionWindows(MD) quantifies mean
%traction from unshifted traction map for cell mask and periphery. For better
%quantification, running windowing package is needed.
% input:
%   MD              MovieData object. Segmentation, QFSM package, and Windowing
%                   package should have been run.
% output:
%   tMapUnshifted:  cell array containing mean traction per layer upto layer 5
%                   It will have a format of tractionCell(windows:layers:frame). Traction at
%                   cell periphery will be from layer 1 or 2.
% Unit is in 
% 
% Sangyoon Han Nov, 2021
% Etienne Michels Jul, 2022

%% Load QFSM package
% nFrames = MD.nFrames_;
% nChannels = length(MD.channels_);
% Load flow vectors
iTract = MD.getProcessIndex('ForceFieldCalculationProcess');
if isempty(iTract)
       disp('Force Field Calculation Process has not been run.');
       throw(MException('quantifyTractionWindows:runForceField','Force Field Calculation Process has not been run.'));
end
flowProcess = MD.getProcess(iTract);
iChan = find(flowProcess.checkChannelOutput);
%% Get the flow
% flow1 = flowProcess.loadChannelOutput(iChan,'output','Md');
% % Load segmented masks
% dt = MD.timeInterval_; 
% res = MD.pixelSize_;
%% Get the window package
iWinPack = MD.getPackageIndex('WindowingPackage');
winPack = MD.getPackage(iWinPack);
%% Get samples from windows
wsProc = winPack.getProcess(4);
for ii=1:numel(wsProc.outFilePaths_(:,iChan))
    [~,outName{ii}]=fileparts(wsProc.outFilePaths_{ii,iChan});
end
%Choose the name that contains Speed
maxNumLayers = 5;
iOutput = find(cellfun(@(x) contains(x,'Traction map unshifted'),outName));
sampleStr = wsProc.loadChannelOutput(iChan,iOutput);
tMapUnshifted = sampleStr.avg(:,1:maxNumLayers,:);            


end

