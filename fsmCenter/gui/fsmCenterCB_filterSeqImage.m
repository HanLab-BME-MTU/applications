function fImg=fsmCentercb_filterSeqImage
% fImg=fsmCentercb_filterSeqImage
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
% Andre Kerstens, 12/15/2004

% Get handle to the fsmCenter gui
hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
handles = guidata(hFsm);

% Get handle to the filter menu
hMenu = findobj('Label','Filter image');

if strcmp(get(hMenu, 'Checked'),'on')
    set(hMenu, 'Checked', 'off');
    assignin('base','sigmaSeq',[]);
else 
    set(hMenu, 'Checked', 'on');

    % Check for the existence of sigma
    if ~exist('sigmaSeq','var') | ~isempty('sigmaSeq')
        % Ask the user to specify the sigma for filtering
        prompt={'Please specify sigma (pixels) for low-pass filtering:'};
        def={'1'};
        dlgTitle='User input requested';
        lineNo=1;
        sigma=char(inputdlg(prompt,dlgTitle,lineNo,def));

        % Check returned value
        if isempty(sigma)
            disp('Aborted by user.');
            set(hMenu, 'Checked', 'off');
            assignin('base','sigmaSeq',[]);
            return
        end
    end

    % Retrieve image from figure
    hImg=findall(gca,'Type','Image');
    img=get(hImg,'CData');

    % Filter
    fImg=Gauss2D(img,str2num(sigma));
    
    % Update figure
    set(hImg,'CData',fImg);
    %set(gca,'CLimMode','manual');
    %set(gca,'CLim',[0 1]);
    refresh;

    % Export to base workspace
    assignin('base','fImg',fImg);
    assignin('base','sigmaSeq',sigma);
    
    % Assign to the handles struct
    handles.imageSeq.sigma = sigma;
    guidata(hFsm, handles);
end
