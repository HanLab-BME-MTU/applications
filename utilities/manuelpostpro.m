function manuelpostpro(hObject)
% manuelpostpro opens a figure and puts a slider with callback changeframe
%               into it
%
% SYNOPSIS       manuelpostpro(hObject)
%
% INPUT          hObject : handle of an object of the GUI calling manuelpostpro
%
% OUTPUT         none
%
% DEPENDENCIES   manuelpostpro uses {nothing}
%                                  
%                manuelpostpro is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04         




% A routine that, together with it's subprograms, allows the user to
% manually postprocess the results of the main analysis. 
% This programm only opens a figure and puts a slider into it. 
% The sliders callback is changeframe, so better look there if you want to
% know more.

%if the slider (and less important, the little windos with the number of
%the current frame) already exist, delete them now.
picturecountH =findall(0,'Style','text','Tag','picturecount');
pictureslideH =findall(0,'Style','slider','Tag','pictureslide');

if ~isempty(picturecountH)
    delete(picturecountH);
end
if ~isempty(pictureslideH)
    delete(pictureslideH);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%fetch and extract values

handles = guidata(hObject)

handles.listofcells='...';


%the good old jobvalues... we need them once more
imagedirectory=handles.jobvalues.imagedirectory;
imagename=handles.jobvalues.imagename;
firstimage=handles.jobvalues.firstimage;
lastimage=handles.jobvalues.lastimage;
increment=handles.jobvalues.increment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




ma=floor((lastimage-firstimage)/increment+0.001);
handles.ma=ma;
slider_step(1) = 1/(ma-1);
slider_step(2) = 3/(ma-1);


guidata(hObject, handles);

figure,

picturecountH=uicontrol('Style','text',...
                        'Units','normalized',...
                        'Tag','picturecount',...
                        'Position',[0.02,0.93,0.05,0.06]);
                     
set(picturecountH,'String',num2str(firstimage));



pictureslideH=uicontrol('Style','slider',...
                        'Units','normalized',... 
                        'Value',0.00000001,...
                        'Min',0,...
                        'Max',1,...
                        'SliderStep',slider_step,...
                        'Callback','changeframe',...
                        'Tag','pictureslide',...
                        'Position',[0.02,0.02,0.05,0.9]);
