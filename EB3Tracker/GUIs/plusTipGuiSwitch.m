function [h]=plusTipGuiSwitch(hObject,eventdata,h,callbackName)

switch callbackName

    case 'chooseProjData'
        % load projData from meta
        [FileName,PathName] = uigetfile('*.mat','Select projData from project meta folder');
        if ~isequal(FileName,0)
            cd(PathName)
            cd ..
            h.dataDir=pwd;

            p=load([PathName filesep FileName]);
            h.projData=p.projData;

            % load tracksFinal from track
            anDir=formatPath(h.projData.anDir);
            trackDir=[anDir filesep 'track'];
            [listOfFiles]=searchFiles('.mat',[],trackDir,0);
            if ~isempty(listOfFiles)
                load([listOfFiles{1,2} filesep listOfFiles{1,1}])
                if ~exist('tracksFinal','var')
                    error('--plusTipTrackViz: tracksFinal missing...');
                end
            else
                error('--plusTipTrackViz: tracksFinal missing...');
            end
            h.tracksFinal=tracksFinal;
        end

    case 'selectSavedRoiPushbutton'
        [FileName,PathName] = uigetfile({'*.*'},'Select roiYX.mat or roiMask.tif');
        if ~isequal(FileName,0)
            if ~isempty(strfind(FileName,'tif'))
                h.roi=imread([PathName FileName]);
            elseif ~isempty(strfind(FileName,'mat'))
                p=load([PathName FileName]);
                h.roi=p.roiYX;
            else
                errordlg('File chosen was not a roiYX.mat or roiMask.tif file. Please try again.','File Error');
            end
        end
    case 'startFrame'
        sFVal=get(hObject,'String');
        if isequal(lower(sFVal),'min')
            h.timeRange(1)=1;
        else
            h.timeRange(1)=str2double(sFVal);
        end
    case 'endFrame'
        eFVal=get(hObject,'String');
        if isequal(lower(eFVal),'max')
            h.timeRange(2)=h.projData.numFrames;
        else
            h.timeRange(2)=str2double(eFVal);
        end
    case 'selectTracksCheck'
        val=get(hObject,'Value');
        if val==0
            h.ask4select=0;
        else
            h.ask4select=1;
        end
    case 'showTracksCheck'
        val=get(hObject,'Value');
        if val==0
            h.showTracks=0;
        else
            h.showTracks=1;
        end
    case 'speedLimitEdit'
        velLimVal=get(hObject,'String');
        if isequal(lower(velLimVal),'max')
            h.velLimit=inf;
        else
            h.velLimit=str2double(velLimVal);
        end
    case 'indivTrackNumbersEdit'
        userInput=get(hObject,'String');
        h.indivTrack=str2num(userInput)';
    case 'plotTracksPush'
        [h.selectedTracks] = plusTipPlotTracks(h.projData,[],...
            h.timeRange,h.img,h.ask4select,...
            h.plotCurrentOnly,h.roi,h.movieInfo);

        if ~isempty(h.selectedTracks)

            temp=vertcat(h.selectedTracks{:});
            h.selectedTracks=unique(temp(:,1));
            [l w]=size(h.selectedTracks);

            if l>w
                h.selectedTracks=h.selectedTracks';
            end

        end
    case 'aviCheckTrackMov'
        val=get(hObject,'Value');
        if val==0
            h.doAvi=0;
        else
            h.doAvi=1;
        end
    case 'selectedTracksDisplay'
        hObject=h.selectedTracksDisplay;
        if isempty(h.selectedTracks)
            set(hObject,'Visible','Off');
        else
            set(hObject,'Visible','On');
        end
        set(hObject,'String',num2str(h.selectedTracks));
    case 'speedMovieButton'
        plusTipSpeedMovie(h.projData,h.timeRange,h.velLimit,h.roi,h.doAvi);
    case 'trackMovieButton'
        plusTipTrackMovie(h.projData,h.indivTrack,h.timeRange,h.roi,h.magCoef,h.showTracks,h.showDetect,h.doAvi);
    case 'aviCheckSpeedMov'
        val=get(hObject,'Value');
if val==0
    h.doAvi=0;
else
    h.doAvi=1;
end
    case 'resetButton'
        closeGUI = h.figure1; %h.figure1 is the GUI figure
 
guiPosition = get(h.figure1,'Position'); %get the position of the GUI
guiName = get(h.figure1,'Name'); %get the name of the GUI
eval(guiName) %call the GUI again
 
close(closeGUI); %close the old GUI
set(gcf,'Position',guiPosition); %set the position for the new GUI

    otherwise
        error('???')

end
