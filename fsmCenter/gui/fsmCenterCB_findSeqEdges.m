function fsmCentercb_findEdges
% fsmCentercb_findEdges
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  A black and white mask (bwMask) and the drawn polygon are
%         returned into Matlab base workspace
%
% Aaron Ponti, 11/18/2003

% Get handle to the invert menu
hMenu = findobj('Label','Find edges');

if strcmp(get(hMenu, 'Checked'),'on')
    set(hMenu, 'Checked', 'off');
    
    % Delete the old edge
    currentH=findall(gca,'Tag','edge');
    if ~isempty(currentH)
        delete(currentH);
    end
else 
    set(hMenu, 'Checked', 'on');

    % Read current project path from fsmCenter
    hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
    if ~isempty(hFsm)
        handles=guidata(hFsm);
    end
    
    % Retrieve image from figure
    hImg=findall(gca,'Type','Image');
    img=get(hImg,'CData');

    % Ask the user to specify the bitdepth
    prompt={'Please specify the bit-depth of your image:'};
    def={'14'};
    dlgTitle='User input requested';
    lineNo=1;
    bitDepth=char(inputdlg(prompt,dlgTitle,lineNo,def));

    if isempty(bitDepth)
        uiwait(msgbox('Edge extraction interrupted by user.','Warning','modal'));
        return
    end

    % Change bitDepth to double
    bitDepth=str2num(bitDepth);

    % Try to extract edges
    try
        img=double(img);

        % Normalize the image
        img=img/(2^bitDepth-1);

        [ans,img_edge,img_bg,edge_pixel,length_edge,frame_pos]=imFindCellEdge(double(img),'',0,'filter_image',1,'img_sigma',1,'bit_depth',2^bitDepth-1);
        figure(gcf);
        hold on
        e1=plot(edge_pixel(:,1),edge_pixel(:,2),'y-'); % Coordinates in edge_pixel are returned a [x y]n
        set(e1,'Tag','edge');
        hold off
        
        % Assign cands data to handles struct
        handles.imageSeq.bitDepth = bitDepth;

        % Update the handles struct
        guidata(hFsm, handles);

        % Swap columns in img_edge from [x y] to [y x] (for compatibility)
        edge_pixel=edge_pixel(:,[2 1]);
        assignin('base','edgePixels',edge_pixel);
        assignin('base','imgEdge',img_edge);
        assignin('base','bwMask',img_bg);
    catch
        uiwait(errordlg('The function failed to retrieve edges.','Error'));
    end
end