function varargout = create3DMovieROI(movieData,varargin)
%UNDER CONSTRUCTION, OBVIOUSLY!

dirName = 'ROI_';
filName = 'ROI info.mat';

%TEEEMMPT TEMP TEMP TEMP
showSteps = true(6,1);
clearProcesses = true;


%Load the movie into imaris. 
imApp = viewMovieImaris(movieData,[],showSteps);


%Get the original image coordinates. Round to avoid numerical error, and
%add one to first coordinate because imaris starts and ends its corodinates at
%outside edge of each pixel
origX =  round([imApp.mDataSet.mExtendMinX imApp.mDataSet.mExtendMaxX] ./ movieData.pixelSize_) + [1 0];
origY =  round([imApp.mDataSet.mExtendMinY imApp.mDataSet.mExtendMaxY] ./ movieData.pixelSize_) + [1 0];
origZ =  round([imApp.mDataSet.mExtendMinZ imApp.mDataSet.mExtendMaxZ] ./ movieData.zSpacing_) + [1 0];


%TEMP - THIS IS DUMB. Wait for the close event from imaris...For now Tell the user what the fuck is going on, and don't continue until they're
%done cropping
announcement = warndlg('WAIT! FIRST, crop the movie in imaris, THEN click OK, and DO NOT save changes when asked',...
                        '***Please read carefully***');
uiwait(announcement);

%Get the cropped image bounding coordinates, rounding to avoid numerical
%error
try
    cropX =  round([imApp.mDataSet.mExtendMinX imApp.mDataSet.mExtendMaxX] ./ movieData.pixelSize_) + [1 0];
    cropY =  round([imApp.mDataSet.mExtendMinY imApp.mDataSet.mExtendMaxY] ./ movieData.pixelSize_) + [1 0];
    cropZ =  round([imApp.mDataSet.mExtendMinZ imApp.mDataSet.mExtendMaxZ] ./ movieData.zSpacing_) + [1 0];
    imClosed = false;
catch    
    imClosed = true;
end
     
% imClosed = false;
% load([movieData.outputDirectory_ filesep 'ROI selection' filesep 'ROI_1.mat']);
% %Adjust for previous error in coordinates:
% cropX = cropX + [1 0];cropY = cropY + [1 0];cropZ = cropZ + [1 0];
% origX = origX + [1 0];origY = origY + [1 0];origZ = origZ + [1 0];

%% ------------ Output ------------ %%

if imClosed || all([cropX cropY cropZ] == [origX origY,origZ])
    disp('No cropping performed, saving nothing.')
else        

    %Get number of current ROI
    iROI = numel(movieData.rois_)+1;
    
    %Setup ROI directory
    outDir = [movieData.outputDirectory_ filesep dirName num2str(iROI)];
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    
    %Save the ROI specification to file    
    outFile = [outDir filesep filName];                
    save(outFile,'cropX','cropY','cropZ','origX','origY','origZ');        
    
    %Create the ROI moviedata and specify info location    
    movieData.addROI(outFile,outDir,~clearProcesses);
    %Save the ROI movieData to the ROI folder
    movieData.rois_(iROI).setPath(outDir);
    %movieData.rois_(iROI).setFilename('ROImovieData.mat');%Can't do this -
    %you're an idiot and hard-coded the filename into the constructor!
    
    %Can't do this either - the objects are handles so it deletes from the
    %parent!
%     %If requested, remove processes from the ROI
%     if clearProcesses
%         nProc = numel(movieData.rois_(iROI).processes_);
%         for j = nProc:-1:1
%             movieData.rois_(iROI).deleteProcess(j);
%         end
%     end
%     
    movieData.rois_(iROI).save;
    
    
end

if nargout > 0
    varargout{1} = imApp;
elseif exist('imApp','var')
    try
        imApp.Quit;
    catch
    end
end
    

