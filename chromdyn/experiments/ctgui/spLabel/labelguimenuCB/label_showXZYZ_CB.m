function zFigH = label_showXZYZ_CB(hObject, eventdata, handles, refreshCall)
% function to show xz/yz maximum projection of 3D movie stacks in labelgui

%is item checked
isChecked = get(hObject,'Checked');



if strcmp(isChecked,'on') %uncheck, close all windows
	set(hObject,'Checked','off');
	zFigH = findall(0,'Tag','XZYZFigure');
	close(zFigH);
    %delete all the saved zFigH
    labelPanelList = findall(0,'Tag','LabelPanel');
    for i = 1:length(labelPanelList)
        SetUserData(labelPanelList(i),[],1,'XZYZFigureH');
    end
	return
	
else %check and open figure
	set(hObject,'Checked','on');

    %find data
    imgFigH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
    if isempty(imgFigH)
        disp('no movie loaded')
        return
        %no movie loaded, keep checked though
    end
	dataProperties = GetUserData(imgFigH,'dataProperties');
	
	if isempty(dataProperties)
		set(hObject,'Checked','on'); %uncheck
        h = errordlg('can''t do this without dataProperties: Please load an idlist first!');
        uiwait(h)
        return
	end
	
	%calc ratio of z vs. xy. Has to be an integer
	pixelRatio = dataProperties.PIXELSIZE_Z/dataProperties.PIXELSIZE_XY;
	pixelRatio = round(pixelRatio);
    if pixelRatio < 1
        pixelRatio = 1;
    end
	
	%get movieSize (x/y/z)
	movieSize = dataProperties.movieSize;
    
    if isfield(dataProperties,'crop')
        if ~isempty(dataProperties.crop)
            movieSize(1:2) = dataProperties.crop(2,1:2);
        end
    end
	
	%create vector with which to blow up z-size
	[dummy, blowUpVector] = ndgrid(1:pixelRatio,movieSize(3):-1:1);
	%reshape, so that blowUpVector becomes [...,2;2;2;2;1;1;1;1] with a ratio of four
	blowUpVector=blowUpVector(:); 
	
% 	%open figure
    zFigH = findall(0,'Tag','XZYZFigure');
    if ishandle(zFigH)
        close(zFigH);
    end
    zFigH = uiViewPanel;
    dataName = get(imgFigH,'Name');
    set(zFigH,'Tag','XZYZFigure','Name',['X/Z&Y/Z-view ',dataName],'HandleVisibility','Callback','NumberTitle','off');
	%turn off menu 'Panels' in uiviewpanel
    uiCh = get(zFigH,'Children');
    set(uiCh(4),'Visible','off');
    
	%------create axes
    
    %calculate sizes for nice plots
    xyRatio = 0.85/sum(movieSize(1:2)); %normalized Units/#of pixels
    zRatio = 0.9/(pixelRatio*movieSize(3));
    
    %set figure size so that it has the right proportions
    figRatio = xyRatio/zRatio;
    figPos = get(zFigH,'Position');
    figPos(4) = figPos(3)*figRatio;
    %make sure projection figure is not too small
    resizeFactor = 800/sum(figPos(3:4));
    if resizeFactor > 1
        figPos(3:4) = figPos(3:4)*resizeFactor;
    end
    set(zFigH,'Position',figPos);
    
    %axes box sizes relative to figure size
    pos1 = [0.05, 0.05, movieSize(1)*xyRatio, 0.9];
    pos2 = [0.1+movieSize(1)*xyRatio, 0.05, movieSize(2)*xyRatio, 0.9];
   
	xzPlotH = subplot(1,2,1); 
	yzPlotH = subplot(1,2,2);
    axesH = [xzPlotH,yzPlotH];
    set(xzPlotH,'Position',pos1);
    set(yzPlotH,'Position',pos2);
  
	%save data
	
	SetUserData(zFigH, axesH,1);
	SetUserData(zFigH, blowUpVector,1);
	SetUserData(zFigH, pixelRatio,1);
	pixelSize=[dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z];
	SetUserData(zFigH, pixelSize,1);
	SetUserData(zFigH, movieSize,1);
    
	%save figureHandle
    SetUserData(imgFigH, zFigH, 1, 'XZYZFigureH');
    
	%draw images;
	if nargin==3 | refreshCall~=1
		label_refreshXZYZ;
	end
	
    

    
end