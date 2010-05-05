function [fsmParam,status]=fsmPrepMain(fsmParam)
% fsmPrepMain is the main function of the fsmPrep module
%
% SYNOPSIS   [fsmParam,status]=fsmPrepMain(fsmParam)
%
% INPUT      fsmParam:   general parameter structure
%
% OUTPUT     fsmParam:   modified (when needed) parameter structure
%            status  :   it decides whether the program should continue after this module
%                        or should stop because of errors;
%                        status is set to 0 (error) in the beginning of a module and will
%                        be set to 1 at the end if the module completed successfully.
%                        status = 1 - if the module completed successfully
%                        status = 0 - if the module did not complete successfully
%
% DEPENDENCES   fsmPrepMain uses {}
%               fsmPrepMain is used by { fsmMain }
%
% Aaron Ponti, October 3rd, 2002


% Set initial module status
status=0;

% Check input parameter
if nargin~=1
    error('Input parameter fsmParam expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

userPath     = fsmParam.main.path;          % Working path (defined at the project level)
imagePath    = fsmParam.main.imagePath;     % Image path (defined at the project level)
imgNumber    = fsmParam.main.imgN;          % Number of image to be processed from the stack
xmin         = fsmParam.main.normMin;       % Lower intensity bound for intensity normalization
xmax         = fsmParam.main.normMax;       % Upper intensity bound for intensity normalization
noiseParam   = fsmParam.main.noiseParam;    % Parameters for the noise model
paramSpeckles= fsmParam.prep.paramSpeckles; % High-order speckle parameters
autoPolygon  = fsmParam.prep.autoPolygon;   % Automatic analisys of the image to extract cell boundaries
drawROI      = fsmParam.prep.drawROI;       % The user draws or loads a ROI to restrict analysis
subpixel     = fsmParam.prep.subpixel;      % significant speckles are localized with subpixel accuracy
psfsigma     = fsmParam.prep.psfSigma;      % true physical sigma of the image point-spread function, caluclated by sigma=0.21*(lambda/NA)/pixelsize
filtersigma  = fsmParam.prep.filterSigma;   % sigma used for the low-pass filtering; except where specifically
% stated differently by the user, filtersigma should have the same value as psfsigma;
% for filtersigma>psfsigma, image information is lost during filtering!!                                            % same value as
projDir = fsmParam.project.path;
edgeDir = fsmParam.project.edge;

oldPath = pwd;

% Change to userPath
cd(userPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMAGE STACK - GET FILES NAMES AND NUMBERS
%
% This part is not run in the case of a batch job
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(fsmParam,'batchJob')
    
    cd(imagePath);
    % The user must select the first image of the stack
    [fName,dirName] = uigetfile(...
        {'*.tif;*.tiff;*.jpg;*.jpeg;*.TIF','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
        '*.tif','TIF files (*.tif)'
        '*.tiff','TIFF files (*.tiff)'
        '*.jpg;','JPG files (*.jpg)'
        '*.jpeg;','JPEG files (*.jpeg)'
        '*.*','All Files (*.*)'},...
        'Select first image');
    if(isa(fName,'char') && isa(dirName,'char'))
        
        % Check that the user did not change the image directory
        %    (remove fileseps)
        if samdir(dirName,imagePath)==0
            uiwait(warndlg('If you really want to analyzes these images, please go back to fsmCenter and change the project settings.','Warning','modal'));
            return % Returns an error (status=0)
        end
        
        imageOne = imread([dirName,fName]);
        % Store image size
        imageSize=size(imageOne);
    else
        return % Returns an error (status=0)
    end
    cd(userPath);
    
    % Check whether the image is 8 bit and the selected bitdepth is higher
    imInfo=imfinfo([dirName,fName]);
    selBitDepth=log2(xmax+1);
    if imInfo.BitDepth==8 && selBitDepth>8
        msg=['You selected ',num2str(selBitDepth),' bit depth for your camera, but the images to be analyzed are only 8 bit. Continue anyway?'];
        button = questdlg(msg,...
            'Warning!','Yes','No','No');
        if strcmp(button,'No')
            return
        end
    end
    
    % Recover all file names from the stack
    outFileList=getFileStackNames([dirName,fName]);
    
    % recover the number of image selected
    [path,body,firstIndex]=getFilenameBody(outFileList{1}); %#ok<ASGLU>
    firstIndex=str2double(firstIndex);
    
    % Recover the number of the last image in the folder
    [path,body,no]=getFilenameBody(outFileList{end}); %#ok<ASGLU>
    
    % Prepare string number format
    s=length(num2str(no));
    strg=sprintf('%%.%dd',s);
    
    % Number of images: n
    if imgNumber==0
        n=length(outFileList);
    else
        n=min(imgNumber,length(outFileList));
        % Crop the portion of interest from the filename list
        outFileList=outFileList(1:n);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CALCULATE FACTORS FOR IMAGE INTENSITY CORRECTIONS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % factors=fsmPrepIntCorrFactors(outFileList,n,[xmin xmax]);
    factors=ones(1,n);   % The intensity correction has been removed; yet, the structure has been maintained
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % FIRST CHECK WHETHER WE HAVE TO LOAD A ROI
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if drawROI==2
        
        % The user is asked to load a ROI
        [userROIbwFileName,userROIbwPath]=uigetfile( ...
            {'userROI.mat'}, ...
            'Please load ''userROI.mat''');
        if ~(isa(userROIbwFileName,'char') && isa(userROIbwPath,'char'))
            return
        end
       
        try
            
            % Try to extract all the info we need
            s=load([userROIbwPath,userROIbwFileName]);
            userROIbw=s.userROIbw;
            userROIpoly=s.userROIpoly; %#ok<NASGU>
            
            % Check dimensions
            if size(userROIbw) ~= [imInfo.Height, imInfo.Width]
                
                % Error - inform the user that he will have to draw the roi
                errorMsg='The selected userROI.mat contains a polygon incompatible with your image size. You will be now asked to draw a ROI.';
                uiwait(errordlg(errorMsg,'Error','modal'));
                
                % Set drawROI=1, the user will draw
                drawROI=1;
                
                % And also update fsmParam
                fsmParam.prep.drawROI=1;
            else
                
                % Save polygon to current project
                eval(['save ',userPath,filesep,'userROI.mat userROIbw userROIpoly;']);
                
            end
            
        catch errMsg
            
            % Error - inform the user that he will have to draw the roi
            uiwait(errordlg(errMsg.message,'Error','modal'));
            
            % Set drawROI=1, the user will draw
            drawROI=1;
            
            % And also update fsmParam
            fsmParam.prep.drawROI=1;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ALLOW THE USER TO DRAW A ROI
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if drawROI==1
        
        % Allow the user to draw a ROI
        img=imreadnd2(outFileList{1},xmin,xmax);
        figure;
        imshow(img,[]);
        set(gcf,'Name','Please draw region of interest');
        set(gcf,'NumberTitle','off');
        
        % Select region of interest
        try
            [userROIbw,x,y]=roipoly;
            userROIpoly=[y x];
        catch
            close(gcf);
            return
        end
        
        % Close current figure
        close(gcf);
        
        % Save polygon to disk
        eval(['save ',userPath,filesep,'userROI.mat userROIbw userROIpoly;']);
        
    end
    
else
    
    % Remove the ''batchJob' field
    fsmParam=rmfield(fsmParam,'batchJob');
    
    % Read needed fields from fsmParam
    n=fsmParam.specific.imageNumber;
    imageSize=fsmParam.specific.imgSize;
    strg=fsmParam.specific.formString;
    outFileList=fsmParam.specific.fileList;
    factors=ones(1,n);
    firstIndex=fsmParam.specific.firstIndex;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% START PROCESSING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing progress bar
h = waitbar(0,'Preprocessing...');

% Added by Lin Ji on Jan. 24, 2005
% Search in the 'edge' directory first to see if edge has already
% been detected using edge tracking 'prPanel'.

bgMaskDir = [projDir filesep edgeDir filesep 'cell_mask'];
if autoPolygon == 1
    if isdir(bgMaskDir)
        dirList = dir(bgMaskDir);
        fileList = {dirList([dirList.isdir] == 0).name};
        bgMaskFileList = fileList(strmatch('mask_',fileList));
        
        %Get the index of the available mask files.
        bgMaskFileIndex = zeros(size(bgMaskFileList));
        for k = 1:length(bgMaskFileList)
            [path,body,no] = getFilenameBody(bgMaskFileList{k});
            bgMaskFileList{k}  = [bgMaskDir filesep bgMaskFileList{k}];
            bgMaskFileIndex(k) = str2double(no);
        end
        [bgMaskFileIndex,sortedI] = sort(bgMaskFileIndex);
        bgMaskFileList = bgMaskFileList(sortedI);
    else
        bgMaskFileList  = {};
        bgMaskFileIndex = [];
    end
    
    if isempty(bgMaskFileList)
        warndlg(['No cell edge mask files has been found in the selected ' ...
            '''edge'' directory: ' edgeDir '. Please run edge tracking first by ' ...
            'click ''Run edge tracker'' button in ''fsmCenter'' ' ...
            'or uncheck ''Automatic img segmentation''.'],'warning','modal');
        cd(oldPath);
        if ishandle(h)
            delete(h);
        end
        return;
    end
    
    %Check if edge detection has been done for the selected images.
    inRangeI = find(bgMaskFileIndex>=firstIndex & ...
        bgMaskFileIndex<=firstIndex+n-1);
    
    if length(inRangeI) < n
        answ = questdlg(['Cell edge mask files are missing for some images ' ...
            'in the selected ''edge'' directory: ' edgeDir ...
            '. Do you want to continue or go back and run edge tracking again by ' ...
            'click ''Run edge tracker'' button in ''fsmCenter''?'], ...
            'warning','Continue','Cancel','Continue');
        
        if strcmp(answ,'Cancel')
            cd(oldPath);
            if ishandle(h)
                delete(h);
            end
            return;
        end
    end
end

%%%%%% End of Lin Ji's change %%%%%%%%

% Go through the stack (n images)
for counter1=1:n
    
    % Index of the current image
    currentIndex=counter1+firstIndex-1;
    
    % Load and normalize the image
    img=imreadnd2(outFileList{counter1},xmin,xmax);
    
    % DON'T DO THAT UNTIL IT IS THOROUGHLY CHECKED
    % Load image
    %img = imread(outFileList{counter1});
    % Compute Ascombes transform
    %img = 2 * sqrt(double(img) + 3/8);
    % Normalize image
    %img = (img - xmin) / (xmax - xmin);

    if autoPolygon==1
        
        % Initialize successCE
        successCE=-1;
        
        % Extract cell outlines (b/w mask)
        %Added by Lin Ji on Jan 24, 2005. Check weather a mask file has
        % already existed.
        if ~isempty(bgMaskFileList)
            %Check if the current image has a 'bgMask' file already.
            rBgMaskInd = find(bgMaskFileIndex==currentIndex);
            
            if ~isempty(rBgMaskInd)
                bwMask = double(imread(bgMaskFileList{rBgMaskInd}));
                successCE = 1;
            end
        end
        
        % If imFindCellEdge returns successCE==-1 create a white mask too
        if successCE==-1
            bwMask=ones(size(img)); % imFindCellEdge failed to retrieve the edge
            %fprintf(1,'Edge extraction failed for frame %s.\n',num2str(currentIndex));
            fprintf(1,'Edge mask file for frame %s is missing.\n',num2str(currentIndex));
            fprintf(1,'Run edge tracker to generate this mask.\n');
        end
        
        % Multiply image with mask (to set background to 0)
        img=img.*bwMask;
    end
     
    
    % Prepare the image for the analysis
    img=fsmPrepPrepareImage(img,factors(counter1),[1 1 0 0; 0 0 imageSize(1) imageSize(2)],filtersigma);
    
    % Statistically test the local maxima to extract (significant) speckles
    fsmPrepMainSecondarySpeckles(img,strg,currentIndex,noiseParam,paramSpeckles,fsmParam);
    
    % Update wait bar
    waitbar(counter1/n,h);
    
    % Force matlab to update the waitbar
    drawnow;
end

% Close waitbar
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPDATING fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsmParam.specific.imgSize=imageSize;
fsmParam.specific.imageNumber=n;
fsmParam.specific.firstIndex=firstIndex;
fsmParam.specific.lastIndex=currentIndex;
fsmParam.specific.formString=strg;
fsmParam.specific.fileList=outFileList;
fsmParam.specific.intCorrFactors=factors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING MODULE STATUS TO 1 AND RETURNING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the status to 1 to mean that the module successfully finished
status=1;
