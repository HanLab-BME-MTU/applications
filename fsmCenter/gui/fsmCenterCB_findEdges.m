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

% Retrieve image from figure
children=get(gca,'Children');
if length(children)>1
    children=children(end);
end
img=get(children,'CData');

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
    plot(edge_pixel(:,1),edge_pixel(:,2),'y-'); % Coordinates in edge_pixel are returned a [x y]n
    % Swap columns in img_edge from [x y] to [y x] (for compatibility)
    edge_pixel=edge_pixel(:,[2 1]);
    assignin('base','edgePixels',edge_pixel);
    assignin('base','imgEdge',img_edge);
    assignin('base','bwMask',img_bg);
   
catch
    
    uiwait(errordlg('The function failed to retrieve edges.','Error'));
   
end
