function make3DspotMovie(inputList,inputProperties,ioOptions)
%MAKE3DSPOTMOVIE creates a movie from 3D spot coordinates and microscope movie data
%
%SYNOPSIS make3DSpotMovie(inputList,inputProperties,ioOptions)
%
%INPUT    inputList: (opt) structure containing the input which is to be plotted
%                 .spots:  nsp-by-5 matrix (nsp: number of spots). Each row
%                          defines [x y z color# size], where color# is a
%                          number used for reading the colormap.
%                          Coordinates are in pixels
%
%                 .image:  3D image to be displayed as 'background'
%
%                 one of both is optional!
%
%
%         inputProperties: (opt) structure with optional description of the input
%                 .colormap:         RGB-colormap {hsv}
%                 .maxColor:         highest value for given for color#
%                 .viewAngles:       [Az, El] in degrees (can be retrieved with 'view')
%                 .imgScalingFactor: {[1 1 1]} scaling factor if pixels in some
%                                    dimensions are 'longer' than the others
%
%         ioOptions: (opt) structure for saving figures or the movie
%                 .save2file:       [{0}/1] whether to save the images to file
%                 .save2filePath :  where to save the file, name
%                 .save2movie:      [0/{1}] whether to save the images
%                                   QT-movie
%                 .save2moviePath : where to save the file, name
%                 .verbose:         [0/{1}] currently without effect,
%                                   because there is no point in not
%                                   looking at the individual frames
%
%c: 1/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%defaults
azimuth = -37.5;%matlab defaults for view-angles
elevation = 30;%matlab defaults for view-angles
cMapLength = 64;
cMap = hsv(cMapLength);
maxColor = [];
imgScalingFactor = [1 1 1];
save2file = 0;
save2filePath = '';
save2fileName = '';
save2movie = 1;
save2moviePath = '';
save2movieName = '';
frameRate = 4;
verbose = 1;
doSpots = 1;
doImage = 1;
ask4view = 1;
figureSize = [360   514   560   420]; %default figure position

oldPath = pwd;

%-----TEST INPUT
if nargin == 0
    loadData = 1;
else
    loadData = 0;
end

if nargin > 0 & ~isempty(inputList)
    %test inputList: one can either omit spots or image but not both
    doSpots = isfield(inputList,'spots');
    doImage = isfield(inputList,'image');
    
    if ~(doSpots | doImage)
        error('neither field spots nor field image in inputList')
    end
    
    %if doSpots: when later finding the good ones, we will check for
    %correct matrix size
    
    %if doImage: we expect an image for all timePoints -> test now
    if doImage
        if isempty(inputList(1).image) | length(size(inputList(1).image))~=3
            error('empty image or not 3-dim')
        end
    end
    
end %test inputList

if nargin > 1 & ~isempty(inputProperties)
    %test inputProperties
    
    %colormap
    if isfield(inputProperties,'colorMap')
        if size(inputProperties.colorMap,2)==3
            cMap = colorMap;
            cMapLength = size(cMap,1);
        else
            error('bad size of colormap (n-by-3)')
        end
    end
    
    %maxColor
    if isfield(inputProperties,'maxColor')
        maxColor = inputProperties.maxColor;
    end
    
    %viewAngle
    if isfield(inputProperties,'viewAngle')
        if length(inputProperties.viewAngle) == 2
            azimuth = inputProperties.viewAngle(1);
            elevation = inputProperties.viewAngle(2);
            ask4view = 0;
        else
            error('viewAngle has to be [az,el] in degrees')
        end
    end
    
    %imgScalingFactor
    if isfield(inputProperties,'imgScalingFactor')
        imgScalingFactor = inputProperties.imgScaligFactor;
        minIsf = min(abs(imgScalingFactor(:)));
        if minIsf < 1
            imgScalingFactor = imgScalingFactor/minIsf;
        end
    else
        ask4imgScalingFactor = 1;
    end
    
end %test inputProperties

if nargin > 2 & ~isempty(ioOptions)
    %save2file
    if isfield(ioOptions,'save2file')
        if ioOptions.save2file == 0
            save2file = 0;
        elseif ioOptions.save2file == 1
            save2file = 1;
            %only now check for path
            if isfield(ioOptions,'save2filePath')
                %check whether we have been given a name or a path
                pathNameLength = length(ioOptions.save2filePath);
                fileSepList = strfind(ioOptions.save2filePath);
                if isempty(fileSepList)
                    %assume we're just being given the name
                    save2filePath = pwd;
                    save2fileName = ioOptions.save2filePath;
                elseif fileSepList(end)==pathNameLength
                    %path ends with filesep, so no name
                    save2filePath = ioOptions.save2filePath;
                    save2fileName = '';
                else
                    %1:lastFileSep-1 = path, rest = filename
                    save2filePath = ioOptions.save2filePath(1:fileSepList(end)-1);
                    save2fileName = ioOptions.save2filePath(fileSepList(end)+1:end);
                end
            end
            
        end
    end
    
    %save2movie
    if isfield(ioOptions,'save2movie')
        if ioOptions.save2movie == 0
            save2movie = 0;
        elseif ioOptions.save2movie == 1
            save2movie = 1;
            %only now check for path
            if isfield(ioOptions,'save2moviePath')
                %check whether we have been given a name or a path
                pathNameLength = length(ioOptions.save2moviePath);
                fileSepList = strfind(ioOptions.save2moviePath);
                if isempty(fileSepList)
                    %assume we're just being given the name
                    save2moviePath = pwd;
                    save2movieName = ioOptions.save2moviePath;
                elseif fileSepList(end)==pathNameLength
                    %path ends with filesep, so no name
                    save2moviePath = ioOptions.save2moviePath;
                    save2movieName = '';
                else
                    %1:lastFileSep-1 = path, rest = filename
                    save2moviePath = ioOptions.save2moviePath(1:fileSepList(end)-1);
                    save2movieName = ioOptions.save2moviePath(fileSepList(end)+1:end);
                end
            end
            
        end
    end
    
    %verbose
    if isfield(ioOptions,'verbose')
        if ioOptions.verbose == 0
            verbose = 0;
        elseif ioOptions.verbose == 1
            verbose = 1;
        end
    end
    
end %test ioOptions


%--------END TEST INPUT---------

%--------GET NAMES
if save2movie & isempty(save2movieName);
    helpTxt = 'Please select filename';
    if isempty(save2moviePath)
        helpTxt = [helpTxt, ' and pathname'];
    else
        cd(save2moviePath)
    end
    helpTxt = [helpTxt, ' to save the QT-movie! If you press ''cancel'', the movie will not be saved'];
    
    %tell the user what's going on
    h = helpdlg(helpTxt,'');
    uiwait(h);
    
    [save2movieName,save2moviePath] = uiputfile('*.mov','save results as QT-movie');
    
    %if user cancelled, nothing will be save                               
    if save2movieName == 0
        save2movie = 0;
    else
        %add extension
        if findstr(save2movieName,'.mov')
            %ok
        else
            save2movieName = [save2movieName,'.mov'];
        end
    end
    
end

if save2file & isempty(save2fileName);
    helpTxt = 'Please select filename';
    if isempty(save2filePath)
        helpTxt = [helpTxt, ' and pathname'];
    else
        cd(save2filePath)
    end
    helpTxt = [helpTxt, ' to save the movie frames! If you press ''cancel'', the frames will not be saved. Naming will be fname_1:fname_n'];
    
    %tell the user what's going on
    h = helpdlg(helpTxt,'');
    uiwait(h);
    
    [save2fileName,save2filePath] = uiputfile('*.jpg','save results as jpegs!');
    
    %if user cancelled, nothing will be save                               
    if save2fileName == 0
        save2file = 0;
        
    end
    
end
%----END GET NAMES


%---------LOAD FILES--------
if loadData
    
    
    %load data and calculate inputList
    [inputList,newImgScalingFactor,newMaxColor,doImage,doSpots] = make3DspotMovieBuildInputList;
    
    if ~isempty(newImgScalingFactor)
        imgScalingFactor = newImgScalingFactor;
    end
    if isempty(maxColor)
        maxColor = newMaxColor;
    end
    
end %if loadFiles
%-----END LOAD FILES--------


%start drawing



%calculate how to put background images, calculate blowUpVectors

if doImage
    %for speed reasons, I just calculate the vectors and everything for the
    %special case of 3D image stacks
    
    if imgScalingFactor(1) ~= imgScalingFactor(2)
        error('option not yet implemented: the x and y direction need to have the same scaling')
    end
    
    %calculate blowUpVector
    
    %calc ratio of z vs. xy. Has to be an integer
    pixelRatio = imgScalingFactor(3)/imgScalingFactor(1);
    pixelRatio = round(pixelRatio);
    if pixelRatio < 1
        pixelRatio = 1;
    end
    %movieSi
    imageSize = size(inputList(1).image);
    
    %create vector with which to blow up z-size
    [dummy, blowUpVector] = ndgrid(1:pixelRatio,1:imageSize(3));
    %...
    blowUpVector=blowUpVector(:); 
    
    %updateImageSize
    imageSize(3) = length(blowUpVector);
        
end %if doImage

%calculate axes limits
if doImage
    %limits are imageSize
    axesXLim = [0,imageSize(1)+1];
    axesYLim = [0,imageSize(2)+1];
    axesZLim = [0,imageSize(3)+1];
    %calculate aspect ratio
    boxSize = imageSize(1:3)+2;
    boxAspectRatio = boxSize/min(boxSize);
    
    %set this empty because we do not calculate it here
    allSpotCoords = [];
else
    %limits are extreme coords with 5% air
    
    %get extreme coords
    %collect vector of all spot coordinates and calc extreme coords
    allSpotCoords = cat(1,inputList.spots);
    allSpotCoords = allSpotCoords(:,1:3);
    xyzSpotMinMax = [min(allSpotCoords,[],1);max(allSpotCoords,[],1)];
    
    xyzMean = mean(xyzSpotMinMax,1);
    xyzMeanMean = [xyzMean;xyzMean];
    
    %create wider extremes
    xyzSpotMinMax = ((xyzSpotMinMax - xyzMeanMean) * 1.1) + xyzMeanMean;
    
    %write axes limits
    %limits are imageSize
    axesXLim = [xyzSpotMinMax(1,1),xyzSpotMinMax(2,1)];
    axesYLim = [xyzSpotMinMax(1,2),xyzSpotMinMax(2,2)];
    axesZLim = [xyzSpotMinMax(1,3),xyzSpotMinMax(2,3)];
    
    %calculate aspect ratio
    boxSize = diff(xyzSpotMinMax,1,1);
    boxAspectRatio = boxSize/min(boxSize);
    
end %if doImage axesLimits

%try to set best possible view
if ask4view & doSpots
    %collect vector of all spot coordinates
    if isempty(allSpotCoords)
        allSpotCoords = cat(1,inputList.spots);
        allSpotCoords = allSpotCoords(:,1:3);
    end
    if doImage %take into account blowUpVector
        allSpotCoords(:,3) = allSpotCoords(:,3)*pixelRatio;
    end
    
    %take only every fifth frame
    allSpotCoords = allSpotCoords([1:5:end],:);
    
    viewPrepFigH = figure('NumberTitle','off','Name','select view');
    plot3(allSpotCoords(:,1),allSpotCoords(:,2),allSpotCoords(:,3),'.','MarkerSize',16);
    line(allSpotCoords(:,1),allSpotCoords(:,2),allSpotCoords(:,3),'Color','k');
    set(gca,'XLim',axesXLim,'YLim',axesYLim,'ZLim',axesZLim,'Box','on','PlotBoxAspectRatio', boxAspectRatio);
    view(azimuth,elevation);
    
    h = helpdlg('Please choose a view angle and figure size (don''t forget this!), THEN close this window');
    %place dialogbox below figure
    %...somehow I don't get it with the positioning of the figure
    %     figPos = get(viewPrepFigH,'Position');
    %     dlgPos = get(h,'Position');
    %     dlgPos(2) = figPos(2) - dlgPos(4);
    set(h,'Position',[320.2500  272.2500  297.7500   79.5000]);
    uiwait(h);
    
    [azimuth,elevation] = view;
    if azimuth < 0
        azimuth = azimuth +360;
    end
    disp(['az and el in case you have to retry: ',num2str(azimuth),' ',num2str(elevation)])
    
    figureSize = get(viewPrepFigH,'Position');
    
    close(viewPrepFigH);
    
end

if doImage
    %create the imgPlotMatrices. The projection itself will always be the
    %same, but the matrix with all equal values will be different depending
    %on the view.
    
    %xy-projection
    [xyPlotMatrixX,xyPlotMatrixY] = ndgrid(1:imageSize(1),1:imageSize(2));
    %check Z-angle
    if elevation < 0
        xyPlotMatrixZ = repmat(imageSize(3)+1,[imageSize(1:2)]);
    else
        xyPlotMatrixZ = zeros(imageSize(1:2));
    end
    
    %xz-projection
    [xzPlotMatrixX,xzPlotMatrixZ] = ndgrid(1:imageSize(1),1:imageSize(3));
    %check azimuth
    if azimuth > 90 & azimuth < 270
        xzPlotMatrixY = zeros(imageSize([1,3]));
    else
        xzPlotMatrixY = repmat(imageSize(2)+1,[imageSize([1,3])]);
    end
    
    %yz-projection
    [yzPlotMatrixY,yzPlotMatrixZ] = ndgrid(1:imageSize(2),1:imageSize(3));
    %check azimuth
    if azimuth < 180
        yzPlotMatrixX = zeros(imageSize([2,3]));
    else
        yzPlotMatrixX = repmat(imageSize(1)+1,[imageSize([2,3])]);
    end
end %if doImage calculate plot-matrices

%calculate sphere data
if doSpots
    sizeSphere = 31; %smoothnes of sphere
    [xSphere,ySphere,zSphere] = sphere(sizeSphere-1);
    cData = ones(sizeSphere,sizeSphere,3);
end

%if we save to movie, launch makeQTMovie
if save2movie
    cd(save2moviePath);
    makeQTMovie('start',save2movieName);
end


%now, finally, loop through inputList, show images and ask for accepting
h = helpdlg(['now the software will generate the images.',...
        ' Please DO NOT change anything about them - just accept or discard individual images.',...
        ' Please DO NOT click on anything else while the movie is being stored - otherwise you will get the wrong image!'],...
    'Please read carefully!');
uiwait(h);

%initialize loop parameters
definitiveAnswer = [];
imageNumber = 1;

%loop
for t = 1:length(inputList)
    
    if ~strcmp(definitiveAnswer,'no')
        
        %show figure and axes
        figH = figure('NumberTitle','off','Name',['Figure ',num2str(imageNumber)],...
            'Position',figureSize,'Color',[0,0,0]);
        axH  = axes('XLim',axesXLim,'YLim',axesYLim,'ZLim',axesZLim,...
            'XTick',axesXLim,'YTick',axesYLim,'ZTick',axesZLim,...
            'XTickLabel',' ', 'YTickLabel',' ', 'ZTickLabel',' ',...
            'TickDir','in','Projection','perspective','Box','on','PlotBoxAspectRatio', boxAspectRatio,...
            'XColor',[1,1,1],'YColor',[1,1,1],'ZColor',[1,1,1],'Color',[0,0,0]);
        view(azimuth,elevation);
        colormap gray;
        hold on;
        
        if doSpots %plot nice spots
            if ~isempty(inputList(t).spots)
                for nsp = 1:size(inputList(t).spots)
                    %get sphere data
                    xs = xSphere*inputList(t).spots(nsp,5) + inputList(t).spots(nsp,2);
                    ys = ySphere*inputList(t).spots(nsp,5) + inputList(t).spots(nsp,1);
                    zs = zSphere*inputList(t).spots(nsp,5) + inputList(t).spots(nsp,3)*pixelRatio;
                    %get color data
                    sColor = cMap(round(inputList(t).spots(nsp,4)/maxColor*cMapLength),:);
                    cData(:,:,1) = 0.5*rand(sizeSphere)+0.5 * sColor(1);
                    cData(:,:,2) = 0.5*rand(sizeSphere)+0.5 * sColor(2);
                    cData(:,:,3) = 0.5*rand(sizeSphere)+0.5 * sColor(3);
                    %plot sphere
                    surf(xs,ys,zs,'CData',cData,'EdgeColor','none','FaceLighting','phong','FaceColor','interp');
                end %for nsp = 1:size(inputList(t).spots)
            end %if ~isempty(inputList(t).spots)
        end %if doSpots %plot nice spots
        
        if doImage
            %read image, adjust intensities
            image = inputList(t).image;
            image = image - min(image(:));
            %image = image - median(image(:)); %does help a little for
            %non-corrected images
            image = image/max(image(:));
            
            %create xy-image and plot
            xyImage = max(image,[],3);
            surf(xyPlotMatrixX,xyPlotMatrixY,xyPlotMatrixZ,xyImage,...
                'FaceColor','texturemap','EdgeColor','none','FaceLighting','none');
            
            %blow up image for other projections
            image = image(:,:,blowUpVector);
            
            %create xz-image and plot
            xzImage = squeeze(max(image,[],2));
            surf(xzPlotMatrixX,xzPlotMatrixY,xzPlotMatrixZ,xzImage,...
                'FaceColor','texturemap','EdgeColor','none','FaceLighting','none');
            
            %create xz-image and plot
            yzImage = squeeze(max(image,[],1));
            surf(yzPlotMatrixX,yzPlotMatrixY,yzPlotMatrixZ,yzImage,...
                'FaceColor','texturemap','EdgeColor','none','FaceLighting','none');
            
        end %if doImage
        
        %and there be light...
        camlight('left'); %(maybe adjust)
        
        if isempty(definitiveAnswer)
            %ask if accept
            opt.Position = [320.2500  272.2500  297.7500   79.5000]; %hard-coded pos
            answer = myQuestdlg('Do you want to include/save this frame?',...
                'Choose wisely','yes','no','yes to all','no to all','yes',opt);
        else
            answer = definitiveAnswer;
        end
        
        switch answer
            case 'yes'
                %call makeQTmovie and/or save frame
                if save2movie
                    %make figure current figure
                    figure(figH);
                    cd(save2moviePath);
                    makeQTMovie('addfigure');
                end
                
                if save2file
                    %make figure current figure
                    figure(figH);
                    %export figure as jpeg
                    cd(save2filePath);
                    fileName = [save2fileName,'_',num2str(imageNumber),'.jpg'];
                    frame = getframe(figH);
                    [I,map] = frame2im(frame);
                    if isempty(map)
                        % RGB image
                        imwrite(I,fileName, 'jpg', 'Quality', ...
                            100);
                    else
                        % Indexed image
                        writejpg_map(fileName, I, map);
                    end
                end
                
                imageNumber = imageNumber+1;
            case 'no'
                %don't keep
            case 'no to all'
                %don't keep and remember
                definitiveAnswer = 'no';
            case 'yes to all'
                %keep and remember
                definitiveAnswer = 'yes';
                
                %call makeQTmovie and/or save frame
                if save2movie
                    %make figure current figure
                    figure(figH);
                    cd(save2moviePath);
                    makeQTMovie('addfigure');
                end
                
                if save2file
                    %make figure current figure
                    figure(figH);
                    %export figure as jpeg
                    cd(save2filePath);
                    fileName = [save2fileName,'_',num2str(imageNumber),'.jpg'];
                    frame = getframe(figH);
                    [I,map] = frame2im(frame);
                    if isempty(map)
                        % RGB image
                        imwrite(I,fileName, 'jpg', 'Quality', ...
                            100);
                    else
                        % Indexed image
                        writejpg_map(fileName, I, map);
                    end
                end
                
                imageNumber = imageNumber+1;
            otherwise %user closed -> no to all
                definitiveAnswer = 'no';
        end
        
        close(figH)
        
    end
    
end %for t = 1:length(inputList)

%if we did a movie (and actually saved frames), finish it now
if save2movie & imageNumber > 1
    makeQTMovie('framerate',frameRate);
    makeQTMovie('finish');
end

cd(oldPath);

%%%%%%%%%%%%%%%  writejpg_map %%%%%%%%%%%%%%%%%
% Like the imwrite routine, but first pass the image data through the indicated
% RGB map. from makeQTMovie
function writejpg_map(name,I,map)

[y,x] = size(I);

% Force values to be valid indexes.  This fixes a bug that occasionally 
% occurs in frame2im in Matlab 5.2 which incorrectly produces values of I 
% equal to zero.
I = max(1,min(I,size(map,1)));

rgb = zeros(y, x, 3);
t = zeros(y,x);
t(:) = map(I(:),1)*255; rgb(:,:,1) = t;
t(:) = map(I(:),2)*255; rgb(:,:,2) = t;
t(:) = map(I(:),3)*255; rgb(:,:,3) = t;

imwrite(uint8(rgb),name,'jpeg','Quality',100);