function h = setupGUI(obj)

% ---------------------------------------
% ### SETUP GUI #########################
% ---------------------------------------

% Check if ROI was already set
if isempty(obj.roiPos)
    obj.roiPos = [100 100 0]; %nm
end

if isempty(obj.roiSize)
    obj.roiSize = [500 500 0]; %nm
end

roiRectPixel = [obj.roiPos(1:2) obj.roiSize(1:2)]*obj.previewScaleFactor;

% Init window
h.fig_SelectROI = figure(1);
windowSizeX = 950;
windowSizeY = min([1000 size(obj.preview,1)*5+200]);
set(h.fig_SelectROI,'Name','Select ROI',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'Position',[960,40,windowSizeX,windowSizeY],...
    'Color',[212 208 200]/255);

% Init "Ok" pushbutton
h.pushbutton_Ok = uicontrol('Style','PushButton',...
    'Units','pixel',...
    'String','Ok',...
    'Position',[10 10 100 20],...
    'Callback',@pushbutton_Ok_CB);

% Init preview image
h.image_Preview = imagesc(obj.preview);
% colormap(hot)
colormap(bone)
axis image;
axis ij;
set(gca,'Box','on',...
    'XTick',[],...
    'YTick',[]);

% Init ROI selection tool
h.imrect_PreviewSelection = imrect(gca,roiRectPixel); % gca: Current axis handler
addNewPositionCallback(h.imrect_PreviewSelection,@imrect_PreviewSelection_CB);

% Add all the editable text boxes and the corresponding labels
h.text_x = uicontrol('Style','text',...
    'Position',[230 30 100 20],...
    'String','X [nm]:');
h.edit_x = uicontrol('Style','edit',...
    'Callback',@edit_x_CB,...
    'Position',[230 10 100 20],...
    'String',num2str(obj.roiPos(1)));

h.text_width = uicontrol('Style','text',...
    'Position',[340 30 100 20],...
    'String','Width [nm]:');
h.edit_width = uicontrol('Style','edit',...
    'Callback',@edit_width_CB,...
    'Position',[340 10 100 20],...
    'String',num2str(obj.roiSize(1)));

h.text_y = uicontrol('Style','text',...
    'Position',[450 30 100 20],...
    'String','Y [nm]:');
h.edit_y = uicontrol('Style','edit',...
    'Callback',@edit_y_CB,...
    'Position',[450 10 100 20],...
    'String',num2str(obj.roiPos(2)));

h.text_height = uicontrol('Style','text',...
    'Position',[560 30 100 20],...
    'String','Height [nm]:');
h.edit_height = uicontrol('Style','edit',...
    'Callback',@edit_height_CB,...
    'Position',[560 10 100 20],...
    'String',num2str(obj.roiSize(2)));

h.text_z = uicontrol('Style','text',...
    'Position',[670 30 100 20],...
    'String','Z [nm]:');
h.edit_z = uicontrol('Style','edit',...
    'Callback',@edit_z_CB,...
    'Position',[670 10 100 20],...
    'String',num2str(obj.roiPos(3)));

h.text_depth = uicontrol('Style','text',...
    'Position',[780 30 100 20],...
    'String','Depth [nm]:');
h.edit_depth = uicontrol('Style','edit',...
    'Callback',@edit_depth_CB,...
    'Position',[780 10 100 20],...
    'String',num2str(obj.roiSize(3)));

% Wait with program execution until the user clicks "Ok"
uiwait(h.fig_SelectROI);

% ---------------------------------------
% ### GUI CALLBACKS #####################
% ---------------------------------------
    function pushbutton_Ok_CB(hObject, eventdata, handles)
        obj.roiPos = obj.roiPos;
        obj.roiSize = obj.roiSize;
        close(h.fig_SelectROI);
    end

    function imrect_PreviewSelection_CB(position)
        % Update the diplayed region if it was specified manually
        roiRectPixel = position;
        obj.roiPos(1:2) = position(1:2)/obj.previewScaleFactor;
        obj.roiSize(1:2) = position(3:4)/obj.previewScaleFactor;
        set(h.edit_x,'String',num2str(obj.roiPos(1)));
        set(h.edit_y,'String',num2str(obj.roiPos(2)));
        set(h.edit_width,'String',num2str(obj.roiSize(1)));
        set(h.edit_height,'String',num2str(obj.roiSize(2)));
    end

    function edit_x_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            roiRectPixel(1) = str2double(get(hObject,'String'))*obj.previewScaleFactor;
            obj.roiPos(1) = str2double(get(hObject,'String'));
            setConstrainedPosition(h.imrect_PreviewSelection,roiRectPixel);
        end
    end

    function edit_y_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            roiRectPixel(2) = str2double(get(hObject,'String'))*obj.previewScaleFactor;
            obj.roiPos(2) = str2double(get(hObject,'String'));
            setConstrainedPosition(h.imrect_PreviewSelection,roiRectPixel);
        end
    end

    function edit_z_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            obj.roiPos(3) = str2double(get(hObject,'String'));
        end
    end

    function edit_width_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            roiRectPixel(3) = str2double(get(hObject,'String'))*obj.previewScaleFactor;
            obj.roiSize(1) = str2double(get(hObject,'String'));
            setConstrainedPosition(h.imrect_PreviewSelection,roiRectPixel);
        end
    end

    function edit_height_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            roiRectPixel(4) = str2double(get(hObject,'String'))*obj.previewScaleFactor;
            obj.roiSize(2) = str2double(get(hObject,'String'));
            setConstrainedPosition(h.imrect_PreviewSelection,roiRectPixel);
        end
    end

    function edit_depth_CB(hObject, eventdata, handles)
        if ~isnan(str2double(get(hObject,'String')))
            obj.roiSize(3) = str2double(get(hObject,'String'));
        end
    end

end



