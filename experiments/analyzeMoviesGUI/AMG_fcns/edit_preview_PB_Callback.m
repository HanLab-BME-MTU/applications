function edit_preview_PB_Callback(hObject, eventdata, handles)
% --- Executes on button press in edit_preview_PB.
% hObject    handle to edit_preview_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%get movie length
movieLength = str2double(get(handles.edit_movieLength_txt,'String'));

%if no filterprm, no filtering/detecting is possible
if ~isfield(handles,'FILTERPRM')
    errordlg('Set NA first!','Insufficient Information');
    return
end

%try to load frame selection. if not there, calculate
if isfield(handles,'previewDefault')
    previewDefault = handles.previewDefault;
else
    stepSize = floor(movieLength/9);
    if (10*stepSize+1) <= movieLength
        stepSize = stepSize + 1; %make sure there are not 11 frames
    end
    previewDefault = {['1:',num2str(stepSize),':',num2str(movieLength)]};
end
    

%select frames to preview
done=0;
while ~done
    inputCell = inputdlg('select a vector of max. 10 frames to preview (e.g. 1:2,5,22)','Preview',1,previewDefault);
    inputString = char(inputCell);
    %test input
    try
        eval(['selFrames=[',inputString,'];']);
        test =  isempty(selFrames)+2*(length(selFrames)<=10)+4*any((max(selFrames)<=movieLength));
        switch test %possible outcomes: 1+2: cancel, 2+4: good, rerun otherwise 
            case 6 %good
                done=1;
            case 3 %canceled
                return
            otherwise %not good
                %rerun
        end
    catch
        %not good
    end %try
end %while

handles.previewDefault = inputCell;
guidata(hObject, handles);

%set colordef to black because it looks cool
colordef black;

%if there is a filtered movie, read selected timepoints, else read from original
%movie and filter. Calculate histogram and plot into figure
previewFig = figure;
set(previewFig,'Name',['Preview for movie ',handles.previewData.movieName]);

numPlots = length(selFrames);
plotCols = min(4,2*numPlots);
plotRows = ceil(numPlots/2); %=2*numPlots/4
%calc sequence in which to plot so that first the first col is filled first
movieSequence = [1:4:2*numPlots-1,3:4:2*numPlots-1];
histoSequence = [2:4:2*numPlots,4:4:2*numPlots];

dataProperties.PATCHSIZE = 7;
dataProperties.CH_MAXNUMINTERV =  1000;
CH_MAXSLOPE = str2double(get(handles.edit_maxslope_txt,'String'));
dataProperties.CH_MAXSLOPE = CH_MAXSLOPE;
dataProperties.FILTERPRM = handles.FILTERPRM;




waitbarHandle=mywaitbar(0,[],numPlots,'preparing preview');

try
    if ~isempty(handles.previewData.filteredMovieName) & get(handles.edit_check_filter,'Value')==0
        for i = 1:numPlots
            filteredImg=readmat(handles.previewData.filteredMovieName,selFrames(i));
            maxProj=max(filteredImg(:,:,:),[],3);
            maxProj=maxProj/max(maxProj(:));
            %open plot window
            figure(previewFig);
            moviePlot(i) = subplot(plotRows,plotCols,movieSequence(i));
            %show movie
            imshow(maxProj);
            title(['t = ',num2str(selFrames(i))]);
            %get histogram: generate subplot first to pass down to spotfind
            histoPlot = subplot(plotRows,plotCols,histoSequence(i));
            %calc and show histogram
            dataProperties.previewHandle = histoPlot;
            [cord(i), mnp] = spotfind(filteredImg,dataProperties);
            %blabla
            axes(moviePlot(i));
            hold on;
            for s=1:length(cord(i).sp)
                plot(cord(i).sp(s).cord(1),cord(i).sp(s).cord(2),'r+')
            end;
            %update waitbar handle
            mywaitbar(i/numPlots,waitbarHandle,numPlots);
        end
    else %there is no filtered movie yet. load unfiltered movie and do background subtraction if necessary
        [rawMovie,handles.previewData.movieName,handles.oldBGState] = correctBackground(handles.previewData.movieName,handles.currentBGState,handles.oldBGState);
        rawMovie = rawMovie(:,:,:,:,selFrames);
        filteredMovie = filtermovie(rawMovie,[handles.FILTERPRM]);
        for i = 1:numPlots
            %read filtered image
            filteredImg = squeeze(filteredMovie(:,:,:,:,i));
            maxProj=max(filteredImg(:,:,:),[],3);
            maxProj=maxProj/max(maxProj(:));
            %open plot window
            figure(previewFig);
            moviePlot(i) = subplot(plotRows,plotCols,movieSequence(i));
            %show movie
            imshow(maxProj);
            title(['t=',num2str(selFrames(i))]);
            %get histogram: generate subplot first to pass down to spotfind
            histoPlot = subplot(plotRows,plotCols,histoSequence(i));
            %calc and show histogram
            dataProperties.previewHandle = histoPlot;
            [cord(i), mnp] = spotfind(filteredImg,dataProperties);
            %blabla
            axes(moviePlot(i));
            hold on;
            for s=1:length(cord(i).sp)
                plot(cord(i).sp(s).cord(1),cord(i).sp(s).cord(2),'r+')
            end;            
            %update waitbar handle
            mywaitbar(i/numPlots,waitbarHandle,numPlots);
        end
    end
catch
    if findstr(lasterr,['Error using ==> get',char(10),'Invalid handle'])
        %do nothing special
        return
    else
        rethrow(lasterror) 
    end
end
close(waitbarHandle); 

%add button for fitting psf
handles.fitButtonH = uicontrol('Style','togglebutton','Units','pixels','Position',[1 1 70 35], 'String','Set PSF',...
    'TooltipString','Click here to start fitting the PSF','Tag','preview_fitButton_TB',...
    'Callback','editPropertiesGUI_preview_psFit(gcbo,[],guidata(findall(0,''Tag'',''EPGUI'')))');

handles.previewPSF.cord = cord;
handles.previewPSF.selFrames = selFrames;
handles.previewPSF.axesList = moviePlot;
handles.previewPSF.figureH = previewFig;
handles.previewPSF.imgSize = size(filteredImg);

%remember coordinates for psFit


guidata(hObject,handles);

%set colordef back to white (does not affect previewFigure)
colordef white;
