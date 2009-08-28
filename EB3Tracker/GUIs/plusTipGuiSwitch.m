function [handles]=plusTipGuiSwitch(hObject,eventdata,handles,callbackName)

% GUIs that use this function:
% plusTipTrackViz - for track plotting and movie making

switch callbackName

    case 'getProjPush'
        if ~isfield(handles,'loadProjList')
            handles.loadProjList=0;
        end

        if handles.loadProjList==0

            if handles.getStr==1
                handles.strList=inputGUI;
                [projList,projPathParsed]=getProj(handles.strList);
            else
                handles.strList=[];
                [projList,projPathParsed]=getProj;
            end
            handles.projList=projList;

        else
            % select one or more projList files
            [handles.projList]=combineProjListFiles(0);
            
            % if "select subset" is checked, ask user for search string(s)
            % and find matches from projList
            if handles.getStr==1
                handles.strList=inputGUI;
                roiDirList=struct2cell(handles.projList);
                roiDirList=roiDirList(2,:)';
                nStr=length(handles.strList);

                projCount=0;
                if ~isempty(roiDirList)
                    tempROI=ones(length(roiDirList),1);
                    for i=1:nStr
                        testStr = handles.strList{i};
                        tempROI=tempROI & cellfun(@(y) ~isempty(strfind(y,lower(testStr))),lower(roiDirList));
                    end

                    matches=find(tempROI);

                    handles.projList=handles.projList(matches);
                end
            end

        end

    case 'selectSavedRoiPushbutton'
        [FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat or roiMask.tif');
        if ~isequal(FileName,0)
            if ~isempty(strfind(FileName,'tif'))
                handles.roi=imread([PathName FileName]);
            elseif ~isempty(strfind(FileName,'mat'))
                p=load([PathName FileName]);
                handles.roi=p.roiYX;
            else
                errordlg('File chosen was not a roiYX.mat or roiMask.tif file. Please try again.','File Error');
            end
        end

    case 'startFrame'
        sFVal=get(hObject,'String');
        if isequal(lower(sFVal),'min')
            handles.timeRange(1)=1;
        else
            handles.timeRange(1)=str2double(sFVal);
        end
        set(hObject,'String',num2str(handles.timeRange(1)));

    case 'endFrame'
        eFVal=get(hObject,'String');
        if isequal(lower(eFVal),'max')
            handles.timeRange(2)=handles.projData.numFrames;
        else
            handles.timeRange(2)=str2double(eFVal);
        end
        set(hObject,'String',num2str(handles.timeRange(2)));

    case 'selectOutputDir'
        outDir=0;
        if isfield(handles,'dataDir')
            dirStart=handles.dataDir;
        else
            dirStart=pwd;
        end
        while isequal(outDir,0)
            outDir=uigetdir(dirStart,'Please select OUTPUT directory for movies');
        end
        handles.projData.movDir=outDir;
        cd(outDir);

    case 'selectTracksCheck'
        handles.ask4select=get(hObject,'Value');

    case 'dualPanelCheck'
        handles.rawToo=get(hObject,'Value');

    case 'showTracksCheck'
        handles.showTracks=get(hObject,'Value');

    case 'speedLimitEdit'
        velLimVal=get(hObject,'String');
        if isequal(lower(velLimVal),'max')
            handles.velLimit=inf;
        else
            handles.velLimit=str2double(velLimVal);
        end

    case 'indivTrackNumbersEdit'
        userInput=get(hObject,'String');
        handles.indivTrack=str2num(userInput)';

    case 'plotTracksPush'
        [handles.selectedTracks] = plusTipPlotTracks(handles.projData,[],...
            handles.timeRange,handles.img,handles.ask4select,...
            handles.plotCurrentOnly,handles.roi,handles.movieInfo);

        if ~isempty(handles.selectedTracks)

            temp=vertcat(handles.selectedTracks{:});
            handles.selectedTracks=unique(temp(:,1));
            [l w]=size(handles.selectedTracks);

            if l>w
                handles.selectedTracks=handles.selectedTracks';
            end

        end

    case 'aviCheckTrackMov'
        handles.doAvi=get(hObject,'Value');

    case 'selectedTracksDisplay'
        hObject=handles.selectedTracksDisplay;
        if isempty(handles.selectedTracks)
            set(hObject,'Enable','Off');
        else
            set(hObject,'Enable','On');
        end
        set(hObject,'String',num2str(handles.selectedTracks));

    case 'speedMovieButton'
        plusTipSpeedMovie(handles.projData,handles.timeRange,handles.velLimit,handles.roi,handles.doAvi);

    case 'trackMovieButton'
        plusTipTrackMovie(handles.projData,handles.indivTrack,handles.timeRange,...
            handles.roi,handles.magCoef,handles.showTracks,handles.showDetect,handles.doAvi,handles.rawToo);

    case 'aviCheckSpeedMov'
        handles.doAvi=get(hObject,'Value');
        
    case 'resetButton'
        closeGUI = handles.figure1; %handles.figure1 is the GUI figure
        guiPosition = get(handles.figure1,'Position'); %get the position of the GUI
        guiName = get(handles.figure1,'Name'); %get the name of the GUI
        eval(guiName) %call the GUI again
        close(closeGUI); %close the old GUI
        set(gcf,'Position',guiPosition); %set the position for the new GUI

    otherwise
        error('???')

end
