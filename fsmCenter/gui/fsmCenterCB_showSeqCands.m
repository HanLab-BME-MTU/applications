function fsmCenterCB_showSeqCands
% fsmCenterCB_showSeqCands
%
% This function is a callback for the "More Tools" menu in the figure
% opened by the pushImShow_Callback function of fsmCenter
%
% This function allows the user to overlay the result of a speckle selection
% step on the currently displayed image.
%
% INPUT   None
%
% OUTPUT  None
%
% Andre Kerstens, 12/15/2004

% Current directory
oldDir=cd;

% Get handle to the invert menu
hMenu = findobj('Label','Show speckle selection');

if strcmp(get(hMenu, 'Checked'),'on')
    set(hMenu, 'Checked', 'off');
    
    % Delete the old cands dots
    currentH=findall(gca,'Tag','cands');
    if ~isempty(currentH)
        delete(currentH);
    end
else 
    set(hMenu, 'Checked', 'on');

    % Read current project path from fsmCenter
    hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
    if ~isempty(hFsm)
        handles=guidata(hFsm);
        
        % Put something here so that we start in the cands dir if there is one
        
        
        settings=get(handles.fsmCenter,'UserData');
        if ~isempty(settings)
            pos=find(strncmp(settings.subProjects,'tack',4));
            if ~isempty(pos)
                tackPath=[settings.projDir,filesep,char(settings.subProjects{pos}),filesep,'cands'];
            end
            if ~isempty(tackPath)
                if isdir(tackPath)
                    cd(tackPath);
                end
            end
        end
    end

    % Get slider and counter info
    hFsmSlider = findall (0, 'Tag', 'pictureSlide');
    hFsmCounter = findall (0, 'Style', 'text', 'Tag', 'pictureCount');
    
    % Fetch the jobvalues and image directory
    imageDirectory = handles.imageSeq.imageDirectory;
    firstImage     = handles.imageSeq.firstImage;
    lastImage      = handles.imageSeq.lastImage;
    imageRange     = handles.imageSeq.imageRange;
    imageNameList  = handles.imageSeq.imageNameList;
    bodyName       = handles.imageSeq.bodyName;
    imageRange     = handles.imageSeq.imageRange;
    
    % Get the current value of the slider, so that we know which frame the user
    % wants to see
    sliderValueTmp = get(hFsmSlider, 'Value');
    sliderValue = round(sliderValueTmp * imageRange);

    % Calculate the frame number to show
    imageNumber = (sliderValue - 1) + firstImage;
    
    % The user should select the cands file corresponding to the current image
    [fileName,dirName] = uigetfile(...
        {'*.mat;','MATLAB .mat files (*.mat)'
        '*.*','All Files (*.*)'},...
        'Select cands###.mat file');
    if(isa(fileName,'char') & isa(dirName,'char'))
        
        % Check whather a subpixel speckle file is selected
        spa = strfind(fileName, '_spa');
        if ~isempty(spa)
            % Remove the _spa extension
            fileName(spa:spa+3) = [];
        end
        
        [path,body,no,ext]=getFilenameBody([dirName,filesep,fileName]);
        
        formatStr = sprintf ('%%.%dd', length(no));
        imageNr = sprintf (formatStr, imageNumber);
        
        if ~isempty(spa)
            % A sub pixel speckle file has been selected
            candsName = [path body imageNr '_spa' ext];
        else
            candsName = [path body imageNr ext];
        end
        
        if ~exist(candsName)
            uiwait(errordlg(['Cands file ' candsName ' does not exist.'],'Error','modal'));
            set(hMenu, 'Checked', 'off');
            return
        end
        
        % Load the correct cands file
        try
            s=load(candsName);
            if ~isempty(spa)
                cands=s.candsSP;
            else
                cands=s.cands;
            end
        catch
            uiwait(errordlg('Invalid cands file.','Error','modal'));
            set(hMenu, 'Checked', 'off');
            return
        end
    else
        set(hMenu, 'Checked', 'off');
        return
    end
        
    % Extract speckle classes
    primary=find([cands.status]==1 & [cands.speckleType]==1);
    secondary=find([cands.status]==1 & [cands.speckleType]==2);
    tertiary=find([cands.status]==1 & [cands.speckleType]==3);
    higher=find([cands.status]==1 & [cands.speckleType]>3);

    % Extract speckle positions
    pPos=reshape([cands(primary).Lmax],2,length([cands(primary).Lmax])/2)';
    sPos=reshape([cands(secondary).Lmax],2,length([cands(secondary).Lmax])/2)';
    tPos=reshape([cands(tertiary).Lmax],2,length([cands(tertiary).Lmax])/2)';
    hPos=reshape([cands(higher).Lmax],2,length([cands(higher).Lmax])/2)';

    % Create a list of colors if the user wants to plot several cands on top of each other
    colors={'r','y','g','c','m'};          

    % Check whether there are already cands plotted
    currentH=findall(gca,'Tag','cands');

    % Select color for the plot depending on how many frames are
    %    already on the figure
    mx=0;
    n=get(currentH,'UserData');
    if ~isempty(n)
        for i=1:length(n)
            if n{i}>mx; mx=n{i}; end
        end
    end
    mx=mx+1;
    indx=mod(mx,5);
    if indx==0; indx=5; end

    % Find the 
    % Plot speckles 
    %    All speckles of a certain type are in one plot -   
    %    this allows to easily change their properties 
    hold on;
    h1=plot(pPos(:,2),pPos(:,1),'.','Color',colors{indx},'MarkerSize',6); % Primary speckles
    set(h1,'Tag','cands'); set(h1,'UserData',str2num(no));
    h2=plot(sPos(:,2),sPos(:,1),'+','Color',colors{indx},'MarkerSize',4); % Secondary peckles
    set(h2,'Tag','cands'); set(h2,'UserData',str2num(no));
    h3=plot(tPos(:,2),tPos(:,1),'^','Color',colors{indx},'MarkerSize',4); % Tertiary speckles
    set(h3,'Tag','cands'); set(h3,'UserData',str2num(no));
    h4=plot(hPos(:,2),hPos(:,1),'*','Color',colors{indx},'MarkerSize',4); % Higher-order speckles
    set(h4,'Tag','cands'); set(h4,'UserData',str2num(no));

    % Title
    hTitle=title('Speckles: . (1st order), + (2nd), ^ (3rd), * (4th and above)'); 
    set(hTitle,'Interpreter','none')

    hold off;
    
    % Assign cands data to handles struct
    handles.imageSeq.candsPath = path;
    handles.imageSeq.candsBody = body;
    handles.imageSeq.candsExt = ext;
    handles.imageSeq.candsNo = no;
    handles.imageSeq.candsSpa = spa;
    
    % Update the handles struct
    guidata(hFsm, handles);
    
    % Return loaded cands to MATLAB base workspace
    assignin('base','cands',cands);
end

% Go back to old directory
cd(oldDir);
