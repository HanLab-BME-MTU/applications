function fsmCenterCB_showCands
% fsmCenterCB_showCands
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% This function allows the user to overlay the result of a speckle selection
% step on the currently displayed image.
%
% INPUT   None
%
% OUTPUT  None
%
% Aaron Ponti, 11/18/2003

% The user should select the cands file corresponding to the current image
[fileName,dirName] = uigetfile(...
    {'*.mat;','MATLAB .mat files (*.mat)'
    '*.*','All Files (*.*)'},...
    'Select cands###.mat file');
if(isa(fileName,'char') & isa(dirName,'char'))
    load([dirName,fileName]);
else
    return
end

% Plot positions
hold on;
for i=1:length(cands)
    if cands(i).status==1
        switch cands(i).speckleType
            case 1, plot(cands(i).Lmax(2),cands(i).Lmax(1),'r.');
            case 2, plot(cands(i).Lmax(2),cands(i).Lmax(1),'y.');
            case 3, plot(cands(i).Lmax(2),cands(i).Lmax(1),'g.');
            otherwise, plot(cands(i).Lmax(2),cands(i).Lmax(1),'c.');
        end
    end
end

% Title
title('Speckle order - Red: 1st, Yellow: 2nd, Green: 3rd, Cyan: 4th and above'); 

% Return loaded cands to MATLAB base workspace
assignin('base','cands',cands);
