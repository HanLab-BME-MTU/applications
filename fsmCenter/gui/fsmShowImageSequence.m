function fsmShowImageSequence
% fsmShowSlidingFrames finds the correct image and shows it in the
% figure opened by fsmCentercb_loadNewSequence
%
% SYNOPSIS       fsmShowImageSequence
%
% INPUT          none (it gets values from the slider created in fsmCentercb_loadNewSequence)
%
% OUTPUT         none (it updates the handles object directly)
%
% DEPENDENCIES   fsmShowImageSequence uses {nothing}
%                                  
%                fsmShowImageSequence is used by { fsmCentercb_loadNewSequence }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Dec 04          Initial release

% Current directory
oldDir=cd;

% Look for objects needed for information
hFsm = findall (0, 'Tag', 'fsmCenter');
hFsmSlider = findall (0, 'Tag', 'pictureSlide');
hFsmCounter = findall (0, 'Style', 'text', 'Tag', 'pictureCount');

% Use the hObject just found to get to the handles structure
handles = guidata (hFsm);

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

% Write the current frame number in the little window above the slider
set (hFsmCounter, 'String', num2str (imageNumber));

% Read the image frame from disk
fileName = char(imageNameList(imageNumber));
image = double(imread([imageDirectory,fileName]));

% Get a handle to the figure
figHandle = gcf;

% Get handle to the filter menu and see whether we have to do some
% filtering
hFiltMenu = findobj('Label','Filter image');
if strcmp(get(hFiltMenu, 'Checked'),'on')
   % Filter
   sigmaSeq = handles.imageSeq.sigma;
   fImg = Gauss2D(image,str2num(sigmaSeq));
   
   % Determine min and max without the borders
   bImg = fImg(sigmaSeq:end-sigmaSeq,sigmaSeq:end-sigmaSeq);
   bMin=min(bImg(:));
   bMax=max(bImg(:));
else
   fImg = image;
end

% Get handle to the invert menu and see whether we have to do some
% inverting
hInvMenu = findobj('Label','Invert image');
if strcmp(get(hInvMenu, 'Checked'),'on')
   % Invert
   iImg = (fImg-max(fImg(:)))/(min(fImg(:))-max(fImg(:)));
   
   if strcmp(get(hFiltMenu, 'Checked'),'on')
      % Determine min and max without the borders
      bImg = iImg(sigmaSeq:end-sigmaSeq,sigmaSeq:end-sigmaSeq);
      bMin=min(bImg(:));
      bMax=max(bImg(:));
   end
else
   iImg = fImg;
end

% Show the frame on the screen in the current figure
% hold on;
% if strcmp(get(hFiltMenu, 'Checked'),'on') %| strcmp(get(hInvMenu, 'Checked'),'on')
%    imshow(iImg,[bMin bMax]);
% else
%    imshow(iImg,[]);
% end
% hold off;

% Delete the axes currently shown in the figure on the screen
imgPtr=findall(gca,'Type','Image');
%delete (gca);
%delete (imgPtr); 

% Update figure
set(imgPtr,'CData',iImg);
set(imgPtr,'CDataMapping','scaled')
set(gca,'CLimMode','auto');
%set(gca,'CLim',[0 1]);
refresh;

% Make sure the axis are correct
axis([1 size(image,2) 1 size(image,1)]);

% Get handle to speckle plot menu and see whether these have to be plotted
% as well
hSpeckleMenu = findobj('Label','Show speckle selection');
if strcmp(get(hSpeckleMenu, 'Checked'),'on')
    % Plot speckles
   
    % Get cands file name info
    candsPath = handles.imageSeq.candsPath;
    candsBody = handles.imageSeq.candsBody;
    candsExt = handles.imageSeq.candsExt;
    
    formatStr = sprintf ('%%.%dd', 3);
    imageNr = sprintf (formatStr, imageNumber);

    candsName = [candsPath candsBody imageNr candsExt];

    if ~exist(candsName)
        uiwait(errordlg(['Cands file ' candsName ' does not exist.'],'Error','modal'));
        %set(hSpeckleMenu, 'Checked', 'off');
        
        % Delete the old cands dots
        currentH=findall(gca,'Tag','cands');
        if ~isempty(currentH)
            delete(currentH);
        end

        return
    end

    
    % Load the correct cands file
    try
        s=load(candsName);
        cands=s.cands;
    catch
        uiwait(errordlg('Invalid cands file.','Error','modal'));
        set(hSpeckleMenu, 'Checked', 'off');
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
    indx = 1;

    % Check whether there are already cands plotted and if yes delete them
    currentH=findall(gca,'Tag','cands');
    if ~isempty(currentH)
        delete(currentH);
    end
% 
%     % Select color for the plot depending on how many frames are
%     %    already on the figure
%     mx=0;
%     n=get(currentH,'UserData');
%     if ~isempty(n) & iscell(n)
%         for i=1:length(n)
%             if n{i}>mx; mx=n{i}; end
%         end
%     end
%     mx=mx+1;
%     indx=mod(mx,5);
%     if indx==0; indx=5; end
% 
    % Find the 
    % Plot speckles 
    %    All speckles of a certain type are in one plot -   
    %    this allows to easily change their properties 
    hold on;
    
    h1=plot(pPos(:,2),pPos(:,1),'.','Color',colors{indx},'MarkerSize',6); % Primary speckles
    set(h1,'Tag','cands'); set(h1,'UserData',str2num(imageNr));
    h2=plot(sPos(:,2),sPos(:,1),'+','Color',colors{indx},'MarkerSize',4); % Secondary peckles
    set(h2,'Tag','cands'); set(h2,'UserData',str2num(imageNr));
    h3=plot(tPos(:,2),tPos(:,1),'^','Color',colors{indx},'MarkerSize',4); % Tertiary speckles
    set(h3,'Tag','cands'); set(h3,'UserData',str2num(imageNr));
    h4=plot(hPos(:,2),hPos(:,1),'*','Color',colors{indx},'MarkerSize',4); % Higher-order speckles
    set(h4,'Tag','cands'); set(h4,'UserData',str2num(imageNr));

    % Title
    hTitle=title('Speckles: . (1st order), + (2nd), ^ (3rd), * (4th and above)'); 
    set(hTitle,'Interpreter','none')

    hold off;
end
    
% Add title
titleStr=[imageDirectory,fileName,' (',num2str(size(image,1)),'x',num2str(size(image,2)),')'];
set(figHandle,'Name',titleStr,'NumberTitle','off');

% Go back to old directory
cd(oldDir);
