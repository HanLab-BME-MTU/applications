function varargout = spotview(varargin)
% SPOTVIEW Application M-file for spotview.fig
%    FIG = SPOTVIEW launch spotview GUI.
%    SPOTVIEW('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 09-Apr-2001 13:51:54

if isa(varargin{1},'struct')  % LAUNCH GUI
    spotcords=varargin{1};
    mov=varargin{2};
	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    
    timeslideH=handles.slider1;
    textH=handles.timestep;
    setappdata(fig,'slist',spotcords);
    setappdata(fig,'movie',mov);
   
    nSlice=size(spotcords,2);
    %set timeslider and edit values
    set(timeslideH,'Min',0.9999999);
    set(timeslideH,'Max',nSlice);
    set(timeslideH,'SliderStep',[1/nSlice 10/nSlice]);
    set(timeslideH,'Value',1);
    %plot first timestep    
    dispspot3D(mov(:,:,:,1,1),spotcords(1).sp);
%     cla;
%     hold on;
%     grid on;
%     axis ij;
%     for i=1:size(spotcords(1).sp,2)
%         if strcmp(spotcords(1).sp(i).type,'spb')
%             mark='-.';
%         elseif strcmp(spotcords(1).sp(i).type,'kin')
%             mark='+:';
%         else
%             mark='r+';
%         end;
%         stem3(spotcords(1).sp(i).cord(1),spotcords(1).sp(i).cord(2),spotcords(1).sp(i).cord(3),mark);
%     end;
%     rotate3d on;
	if nargout > 0  
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
textH=handles.timestep;
spotcords=getappdata(gcbf,'slist');
mov=getappdata(gcbf,'movie');
val=round(get(h,'Value'));
set(h,'Value',val);
set(textH,'String',num2str(val));
dispspot3D(mov(:,:,:,1,val),spotcords(val).sp);

% spotsfigh=findobj(0,'Tag','spotdisp');
% %plot spots
% figure(spotsfigh)
% cla;
% hold on;
% grid on;
% axis ij;
% for i=1:size(spotcord(val).sp,2)
%     if strcmp(spotcord(val).sp(i).type,'spb')
%         mark='-.';
%     elseif strcmp(spotcord(val).sp(i).type,'kin')
%         mark='+:';
%     else
%         mark='r+';
%     end;
%    stem3(spotcord(val).sp(i).cord(1),spotcord(val).sp(i).cord(2),spotcord(val).sp(i).cord(3),mark);
% end;
% daspect([19.6 19.6 5]);
% rotate3d on;