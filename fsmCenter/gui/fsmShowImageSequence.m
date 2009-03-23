function fsmShowImageSequence
% fsmShowSlidingFrames finds the correct image and shows it in the
% figure opened by fsmCenterCB_loadNewSequence
%
% SYNOPSIS       fsmShowImageSequence
%
% INPUT          none (it gets values from the slider created in fsmCenterCB_loadNewSequence)
%
% OUTPUT         none (it updates the handles object directly)
%
% DEPENDENCIES   fsmShowImageSequence uses {nothing}
%                                  
%                fsmShowImageSequence is used by { fsmCenterCB_loadNewSequence }
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
sliderValue = get(hFsmSlider, 'Value');
sliderValue = round(sliderValue * imageRange);

% Calculate the frame number to show
imageNumber = (sliderValue - 1) + firstImage;

% Write the current frame number in the little window above the slider
set (hFsmCounter, 'String', num2str (imageNumber));

% Read the image frame from disk
imgListNr = sliderValue;
fileName = char(imageNameList(imgListNr));
image = double(imread([imageDirectory,fileName]));

% Get a handle to the figure
figHandle = gcf;

% Get handle to the filter menu and see whether we have to do some
% filtering
hFiltMenu = findobj('Label','Filter image');
if strcmp(get(hFiltMenu, 'Checked'),'on')
   % Filter
   sigmaSeq = handles.imageSeq.sigma;
   fImg = Gauss2D(image,str2double(sigmaSeq));
   
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

% Delete the axes currently shown in the figure on the screen
imgPtr=findall(gca,'Type','Image');

% Update figure
set(imgPtr,'CData',iImg);
set(imgPtr,'CDataMapping','scaled')
set(gca,'CLimMode','auto');
refresh;

% Make sure the axis are correct
axis([1 size(image,2) 1 size(image,1)]);

% See if the edge finder has to be run
hEdgeMenu = findobj('Label','Find edges');
if strcmp(get(hEdgeMenu, 'Checked'),'on')
    % Retrieve image from figure
    hImg=findall(gca,'Type','Image');
    img=get(hImg,'CData');

    % Get bitdepth from handles struct
    bitDepth = handles.imageSeq.bitDepth;

    % Check whether there are already edges plotted and if yes delete them
    currentH=findall(gca,'Tag','edge');
    if ~isempty(currentH)
        delete(currentH);
    end
    
    % Try to extract edges
    try
        img=double(img);

        % Normalize the image
        img=img/(2^bitDepth-1);

        [ans,img_edge,img_bg,edge_pixel,length_edge,frame_pos]=imFindCellEdge(double(img),'',0,'filter_image',1,'img_sigma',1,'bit_depth',2^bitDepth-1);
        figure(gcf);
        hold on
        e1=plot(edge_pixel(:,1),edge_pixel(:,2),'y-'); % Coordinates in edge_pixel are returned a [x y]n
        set(e1,'Tag','edge');
        hold off
        
        % Swap columns in img_edge from [x y] to [y x] (for compatibility)
        edge_pixel=edge_pixel(:,[2 1]);
        assignin('base','edgePixels',edge_pixel);
        assignin('base','imgEdge',img_edge);
        assignin('base','bwMask',img_bg);
    catch
        uiwait(errordlg('The function failed to retrieve edges.','Error'));
    end
end

% Get handle to speckle plot menu and see whether these have to be plotted
% as well
hSpeckleMenu = findobj('Label','Show speckle selection');
if strcmp(get(hSpeckleMenu, 'Checked'),'on')
    % Plot speckles
   
    % Get cands file name info
    candsPath = handles.imageSeq.candsPath;
    candsBody = handles.imageSeq.candsBody;
    candsExt = handles.imageSeq.candsExt;
    candsNo = handles.imageSeq.candsNo;
    candsSpa = handles.imageSeq.candsSpa;
    
    formatStr = sprintf ('%%.%dd', length(candsNo));
    imageNr = sprintf (formatStr, imageNumber);

    if ~isempty(candsSpa)
        % A sub pixel speckle file has been selected
        candsName = [candsPath candsBody imageNr '_spa' candsExt];
    else
        candsName = [candsPath candsBody imageNr candsExt];
    end
    
    if ~exist(candsName, 'file')
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
        if ~isempty(candsSpa)
            cands=s.candsSP;
        else
            cands=s.cands;
        end
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

    % Find the 
    % Plot speckles 
    %    All speckles of a certain type are in one plot -   
    %    this allows to easily change their properties 
    hold on;
    
    h1=plot(pPos(:,2),pPos(:,1),'.','Color',colors{indx},'MarkerSize',6); % Primary speckles
    set(h1,'Tag','cands'); set(h1,'UserData',str2double(imageNr));
    h2=plot(sPos(:,2),sPos(:,1),'+','Color',colors{indx},'MarkerSize',4); % Secondary peckles
    set(h2,'Tag','cands'); set(h2,'UserData',str2double(imageNr));
    h3=plot(tPos(:,2),tPos(:,1),'^','Color',colors{indx},'MarkerSize',4); % Tertiary speckles
    set(h3,'Tag','cands'); set(h3,'UserData',str2double(imageNr));
    h4=plot(hPos(:,2),hPos(:,1),'*','Color',colors{indx},'MarkerSize',4); % Higher-order speckles
    set(h4,'Tag','cands'); set(h4,'UserData',str2double(imageNr));

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
