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

% Extract speckle classes
primary=find([cands.status]==1 & [cands.speckleType]==1);
secondary=find([cands.status]==1 & [cands.speckleType]==2);
tertiary=find([cands.status]==1 & [cands.speckleType]==3);
higher=find([cands.status]==1 & [cands.speckleType]>3);

% Extract speckle positions
pPos=reshape([cands(primary).Lmax],2,length([cands(primary).Lmax])/2)';
sPos=reshape([cands(secondary).Lmax],2,length([cands(secondary).Lmax])/2)';
tPos=reshape([cands(tertiary).Lmax],2,length([cands(tertiary).Lmax])/2)';
hPos=reshape([cands(higher).Lmax],2,length([cands(higher).Lmax])/2)';

% Plot speckles 
%    All speckles of a certain type are in one plot -   
%    this allows to easily change their properties 
plot(pPos(:,2),pPos(:,1),'ro'); % Primary speckles
plot(sPos(:,2),sPos(:,1),'yo'); % Secondary peckles
plot(tPos(:,2),tPos(:,1),'go'); % Tertiary speckles
plot(hPos(:,2),hPos(:,1),'co'); % Higher-order speckles

% Title
title('Speckle order - Red: 1st, Yellow: 2nd, Green: 3rd, Cyan: 4th and above'); 

% Return loaded cands to MATLAB base workspace
assignin('base','cands',cands);
