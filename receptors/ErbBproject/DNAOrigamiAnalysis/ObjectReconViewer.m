function handles = ObjectReconViewer(ReconObject, varargin)
% ObjectReconViewer creates an instance of an interactive GUI for examning the results of fourEdgeWrapepr or lineWrapper
%
% inputs: 
%       ReconObject  -> a cell array containing structures
%                       with a .counts field that will be displayed.
%
% outputs: 
%       handles  -> a structure that contains all of the object handles of
%                   the gui
%
% Last Update: 2013_01_21
% Jeffrey L Werbin
% Harvard Medical School

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('ReconObject',@iscell);
ip.addOptional('field','counts',@ischar);
ip.parse(ReconObject,varargin{:});

field = ip.Results.field;

% makes object for figure
f = figure;
a = axes('DataAspectRatio',[1,1,1]);
a = get(f,'Children');
s = uicontrol('Style','slider','Min',1,'Max',numel(ReconObject), ...
              'Value',1,'Position',[10,10,150,30], ...
              'Callback',{@slideradjust}, 'SliderStep',[1/numel(ReconObject),10/numel(ReconObject)]);
t = uicontrol('Style','text','Position',[10,40,80,20],'String', 'N = 1');

%set state variables
state = 1; %starting position for display

%set gui data for object usage
slider = struct('a',a,'t',t,'i',state);
slider.data = ReconObject;
setappdata(s,'slider',slider);
setappdata(s,'field',field);

%creates output structure;
handles = struct('figure', f,'axes',a,'slider',s, 'text',t);

%itialize plot
figure(f);
    if strcmp(field,'img')
        imagesc(slider.data{state}.(field)(:,:,2),[0,8]);
        colormap('pink');
        imOverlayPointMap(slider.data{state}.(field)(:,:,1),'g+');
    else
        imshow(slider.data{state}.(field),[0,8]);
        colormap('hot');
    end
end


function slideradjust(hObject, eventData, handles)
    %get data stored in slider
    st = getappdata(hObject,'slider'); 
    field = getappdata(hObject,'field');
    st.i = fix(get(hObject,'Value'));
    axes(st.a);
    
    if strcmp(field,'img')
        imagesc(st.data{st.i}.(field)(:,:,2),[0,8]);
        colormap('pink');
        imOverlayPointMap(st.data{st.i}.(field)(:,:,1),'g+');
    else
        imshow(st.data{st.i}.(field),[0,8]);
        colormap('hot');
    end
    
    set(st.t,'String',['N = ',num2str(st.i,'%u')]);
    
    %update data stored in slider
    setappdata(hObject,'slider',st);
end