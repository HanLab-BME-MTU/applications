function fsmCentercb_drawPolygon
% fsmCentercb_loadPolygon
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% INPUT   None
%
% OUTPUT  None
%
% REMARK  The user is asked to load a polygon file. The polygon must be stored as a series 
%         of coordinates [y x]n 
%
% Aaron Ponti, 04/23/2004

% The user should pick a polygon file
[fName,dirName] = uigetfile(...
    {'*.mat','*.mat'},...
    'Select a polygon file');
if(isa(fName,'char') & isa(dirName,'char'))
    polygon=load([dirName,fName]);
else
    return 
end

% Plot the polygon
polygonFields=fieldnames(polygon);
if length(polygonFields)>1
    error('Invalid polygon');
end
coords=getfield(polygon,char(polygonFields));
if size(coords,2)~=2
    if size(coords,1)==2
        coords=coords';
    else
        error('The coordinates must be stored in an [y x]n matrix.');
    end
end

% Plot coordinates on current image
handle=findall(0,'Tag','ViewPanel');
if isempty(handle)
    return;
end
figure(handle);
hold on;
plot(coords(:,2),coords(:,1),'r-');

% Return polygon to MATLAB base workspace
assignin('base','polygon',coords);