function fImg=fsmCentercb_filterImage
% fImg=fsmCentercb_filterImage
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% INPUT   none
%
% OUTPUT  none
%
% REMARK  The image is low-pass filtered with a user-defined sigma
%         and returned into the MATLAB bae workspace. The image in the 
%         figure will be updated.
%
% Aaron Ponti, 11/18/2003

% Ask the user to specify the sigma for filtering
prompt={'Please specify sigma (pixels) for low-pass filtering:'};
def={'1'};
dlgTitle='User input requested';
lineNo=1;
sigma=char(inputdlg(prompt,dlgTitle,lineNo,def));

% Retrieve image from figure
img=get(get(gca,'Children'),'CData');

% Filter
fImg=Gauss2D(img,str2num(sigma));

% Update figure
set(get(gca,'Children'),'CData',fImg);
refresh;

% Export to base workspace
assignin('base','fImg',fImg);
