function [handles]=plusTipGuiSwitch(hObject,eventdata,handles,callbackName)


switch callbackName

    case 'getProjPush'
        if ~get(handles.getProjListFile_check,'Value')
            % If load projList is not checked
            if get(handles.getQueryStr_Check,'Value')
                handles.strList=inputGUI;
                projList=getProj(handles.strList);
            else
                handles.strList=[];
                projList=getProj;
            end
            handles.projList=projList;
            handles.nProjLists = 0;
        else
            % select one or more projList files
            [handles.projList handles.nProjLists]=combineProjListFiles(0);
            
            % if "select subset" is checked, ask user for search string(s)
            % and find matches from projList
            if get(handles.getQueryStr_Check,'Value')
                handles.strList=inputGUI;
                roiDirList=projList2Cell(handles.projList);
                nStr=length(handles.strList);

                if ~isempty(roiDirList)
                    tempROI=ones(length(roiDirList),1);
                    for i=1:nStr
                        testStr = handles.strList{i};
                        tempROI=tempROI & cellfun(@(y) ~isempty(strfind(y,lower(testStr))),lower(roiDirList(:,2)));
                    end
                    
                    matches=find(tempROI);
                    handles.projList=handles.projList(matches);
                end
            end
        end

    case 'selectedTracksDisplay'
        hObject=handles.selectedTracksDisplay;
        if isempty(handles.selectedTracks)
            set(hObject,'Enable','Off');
        else
            set(hObject,'Enable','On');
        end
        set(hObject,'String',num2str(handles.selectedTracks));


    case 'resetButton'
        closeGUI = handles.figure1; %handles.figure1 is the GUI figure
        guiPosition = get(handles.figure1,'Position'); %get the position of the GUI
        guiHandle = str2func(get(handles.figure1,'Name')); %get the name of the GUI
        close(closeGUI); %close the old GUI before calling it (singleton)
        guiHandle('Position',guiPosition) %call the GUI again        
        
    otherwise
        error('???')

end
