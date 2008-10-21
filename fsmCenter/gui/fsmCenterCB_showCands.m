function fsmCenterCB_showCands
% fsmCenterCB_showCands
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
% Aaron Ponti, 11/18/2003

% Current directory
oldDir=cd;

% Read current project path from fsmCenter
hFsm=findall(0,'Tag','fsmCenter','Name','fsmCenter');
if ~isempty(hFsm)
    handles=guidata(hFsm);
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

% The user should select the cands file corresponding to the current image
[fileName,dirName] = uigetfile(...
    {'*.mat;','MATLAB .mat files (*.mat)'
    '*.*','All Files (*.*)'},...
    'Select cands###.mat file');
if(isa(fileName,'char') && isa(dirName,'char'))
    try
        s=load([dirName,fileName]);
        %To cope with subpixel cands.
        if isfield(s,'candsSP')
            cands = s.candsSP;
        else
            cands=s.cands;
        end
    catch
        uiwait(errordlg('Invalid cands file.','Error','modal'));
        return
    end
else
    return
end

% Index
%To cope with subpixel cands filename. e.g. cands001_spa.
ind = strfind(fileName,'_spa');
if ~isempty(ind)
    fileName = fileName(1:ind(end)-1);
end
[path,body,no,ext]=getFilenameBody([dirName,filesep,fileName]);

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
set(h1,'Tag','cands'); set(h1,'UserData',str2double(no));
h2=plot(sPos(:,2),sPos(:,1),'+','Color',colors{indx},'MarkerSize',4); % Secondary peckles
set(h2,'Tag','cands'); set(h2,'UserData',str2double(no));
h3=plot(tPos(:,2),tPos(:,1),'^','Color',colors{indx},'MarkerSize',4); % Tertiary speckles
set(h3,'Tag','cands'); set(h3,'UserData',str2double(no));
h4=plot(hPos(:,2),hPos(:,1),'*','Color',colors{indx},'MarkerSize',4); % Higher-order speckles
set(h4,'Tag','cands'); set(h4,'UserData',str2double(no));

% Title
hTitle=title('Speckles: . (1st order), + (2nd), ^ (3rd), * (4th and above)'); 
set(hTitle,'Interpreter','none')

% Return loaded cands to MATLAB base workspace
assignin('base','cands',cands);

% Go back to old directory
cd(oldDir);
