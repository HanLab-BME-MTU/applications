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
img=get(get(gca,'Children'),'CData');

try
    img=double(img);
    
    % Normalize the image
    bitdepth=14;
    img=img/(2^bitdepth-1);
    
    [ans,img_edge,img_bg,edge_pixel,length_edge,frame_pos]=imFindCellEdge(double(img),'',0,'bit_depth',bitdepth);
    figure(gcf);
    hold on
    plot(edge_pixel(:,1),edge_pixel(:,2),'y-');
    assignin('base','imgEdge',img_edge);
    assignin('base','bwMask',img_bg);
    assignin('base','edgePixels',edge_pixel);
   
catch
    
    uiwait(errordlg('The function failed to retrieve edges.','Error'));
   
end
