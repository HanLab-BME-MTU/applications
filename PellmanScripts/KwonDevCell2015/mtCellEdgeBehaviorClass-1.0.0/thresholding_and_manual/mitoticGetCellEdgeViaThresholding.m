function [ output_args ] = mitoticGetCellEdgeViaThresholding( ML,varargin)
% mitoticGetCellEdgeViaThresholding: segment via thresholding the cell edge. 
%
%
%% 
% INPUT: 
% movieList: (REQUIRED) : a movieData or movieList object (created in
%                         movieSelectorGUI.m 
%
% outFilename: (PARAM) : character array
%                        Default: 'masks'
%                        Name of file for saving masks 
%
% outPath: (PARAM) :     character array 
%                        Default: MD.outputDirectory_ for individual
%                        projects
%         
% threshMethod: (PARAM) : character array  ('Rosin','MinMax','Otsu',[])
%                        Default : [] (flags to make a troubleshoot overlay)
%                        where you can choose your segmentation of choice.
%                        Note only looks at first frame but this was
%                        typically sufficient. 
% 
% largestOnly : (PARAM) : logical  
%                        Default : true 
%                        Keep only the largest connected component of the
%                        thesholded mask 
%
% imopenMask: (PARAM) : logical 
%                        Default : false 
%                        Smooth mask via a imopen : Default is to use strel 
%                        'disk'
% 
% radiusDisk: (PARAM) : scalar 
%                       Default : 3 
%                       Radius for the strel 'disk'
%                       Only applicable if imopenMask set to true

%
% % OUTPUT: Creates a file 'masks' in MD.outputDirectory_ containing
%        threshold masks for the entire movie. 
%        Note no smoothing of the edge is
%        performed by default. However, the largest conneected component is chosen by default.
%        By default the inside the 'masks' folder a folder 'overlays'is likewise created
%        where the mask outline is plotted over the original image for each
%        frame. These options, along with the thresholding method used, can be modified : see above                    
%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('outFilename','masks',@(x) ischar(x));

ip.addParameter('outPath',[]); 

ip.addParameter('threshMethod',[]);

ip.addParameter('TSOverlays',true); 

ip.addParameter('largestOnly',true); %

% Smoothing via opening parameters
ip.addParameter('imopenMask',false); 
ip.addParameter('radiusDisk',3); 

ip.parse(varargin{:});
%% Set up 
if isa(ML,'MovieData')
    projList{1,1} = [ML.movieDataPath_ filesep ML.movieDataFileName_];
else
    projList = ML.movieDataFile_;
end

if ip.Results.imopenMask 
    forFilename = ['SmoothViaOpen_DiskRad' num2str(ip.Results.radiusDisk)]; 
else 
    forFilename = []; 
end 
%% Make interactive overlay of the three different thresholding methods if required
% NOTE: currently will not loop only looks at the first frame of the 
% first project- typically this was completley sufficient to find the
% correct method for the dataset. 
if isempty(ip.Results.threshMethod)    
    load(projList{1}); 
    % ask user which type to use
    names = {'Rosin','MinMax','Otsu'};
    stopFlag = 0;
    while stopFlag == 0
        
        [s,v] = listdlg('PromptString',{'Please Select'; 'a Threshold Method';'To Compare Methods' ; 'Select All'},...
            'ListString',names);
        
        if v == 1
            if length(s) == 1
                threshMethod = names{s} ;
                close gcf
                stopFlag = 1;
            else
                % get first image
                % do calc for both overlay
                %
                img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{1}])) ;
   
                img = img./((2^MD.camBitdepth_)-1);
              
                % filter the image 
                imgFiltered = filterGauss2D(img,1);
                maxSignal = max(imgFiltered(:));
                minSignal = min(imgFiltered(:));
                %normalize image between 0 and 1
                imgFilteredNorm = (imgFiltered - minSignal)/(maxSignal - minSignal);
                
                [levelRosin]  = thresholdRosin(imgFilteredNorm);
                
                [levelMaxMin] = thresholdFluorescenceImage(imgFilteredNorm);
                
                [levelOtsu] = graythresh(imgFilteredNorm);
                
                maskR = imgFilteredNorm>=levelRosin; % get mask
                maskM = imgFilteredNorm>=levelMaxMin;
                maskO = imgFilteredNorm>=levelOtsu;
               
                
                if ip.Results.imOpenMask
                    maskM = imopen(maskM,strel('disk',ip.Results.radiusDisk));
                    maskR = imopen(maskR,strel('disk',ip.Results.radiusDisk));
                    maskO = imopen(maskO,strel('disk',ip.Results.radiusDisk));
                end
                
                
                % clean it up
                maskR = imfill(maskR,'holes');
                maskM = imfill(maskM,'holes');
                maskO= imfill(maskO,'holes');
                
                if ip.Results.largestOnly
                    % get largest CC
                    maskR = logical(getLargestCC(maskR));
                    maskM = logical(getLargestCC(maskM));
                    % maskDir = [projList{iProj} filesep 'masks'];
                    maskO = logical(getLargestCC(maskO));
                end
                roiYXR = bwboundaries(maskR);
                roiYXM = bwboundaries(maskM);
                roiYXO = bwboundaries(maskO);
                figure;
                subplot(1,2,1)
                imshow(-img,[]);
                hold on
                cellfun(@(x) plot(x(:,2),x(:,1),'b'),roiYXR);
                cellfun(@(x) plot(x(:,2),x(:,1),'m'),roiYXM);
                cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYXO);
                text(5,10,'Rosin','color','b');
                text(5,25,'MinMax','color','m');
                text(5,40,'Otsu','color','g');
                subplot(1,2,2)
                hist(imgFilteredNorm(:),100)
                xlabel('Normalized Img Intensity')
                ylabel('Number of Pixels');
                
                hold on
                [x,~] = hist(imgFilteredNorm(:),100);
                
                line([levelRosin,levelRosin],[0,max(x)],'color','b');
                line([levelMaxMin,levelMaxMin],[0 max(x)],'color','m');
                line([levelOtsu,levelOtsu],[0 max(x)],'color','g');
                ms =  msgbox('Examine Thresholding: Zoom if Needed, Click OK When Finished');
                uiwait(ms)
            end
            
        else
            x = questdlg('Try Again?');
            if (strcmpi(x,'no') || strcmpi(x,'cancel'))
                stopFlag =1;
                error('User Exited: No Analysis Performed');
                
            end
        end
    end
else 
    threshMethod = ip.Results.threshMethod; 
end % isempty ThreshMethod


%% Thresholding

% Start movieList loop 
for iProj  =1 : numel(projList)
    load(projList{iProj}); 
    display(['Performing ' threshMethod 'Thresholding for ' projList{iProj}])
    
    % Set up mask Directory
        if isempty(ip.Results.outPath)
            maskDir = [MD.outputDirectory_ filesep ip.Results.outFilename];
        else
            maskDir = [ip.Results.outPath filesep ip.Results.outFilename];
        end

        if ~isdir(maskDir)
            mkdir(maskDir); 
        end 
        
        % Make overlay directories if necessary
        if ip.Results.TSOverlays   
            if ~isdir([maskDir filesep 'overlays'])    
                mkdir([maskDir filesep 'overlays']);
            end
        end 

  % Start thresholding individual movie  
    for iMask = 1:MD.nFrames_;
        
        % Load image 
        fileNameIm = [MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iMask}];
        img = double(imread(fileNameIm));
        img = img./((2^MD.camBitdepth_)-1);
        
        % filter image 
        imgFiltered = filterGauss2D(img,1);
        maxSignal = max(imgFiltered(:));
        minSignal = min(imgFiltered(:));
        
        %normalize image between 0 and 1
        imgFilteredNorm = (imgFiltered - minSignal)/(maxSignal - minSignal);
        if strcmpi(threshMethod,'Rosin')
            
            [level]  = thresholdRosin(imgFilteredNorm);
        elseif strcmpi(threshMethod,'MinMax')
            [level] = thresholdFluorescenceImage(imgFilteredNorm);
        else
            [level] = graythresh(imgFilteredNorm);
        end
        
        % create the mask 
        mask = imgFilteredNorm>=level; % get mask
       
        if ip.Results.imopenMask
            mask = imopen(mask,strel('disk',ip.Results.radiusDisk));
        end
        % clean it up
        mask = imfill(mask,'holes');
        
        if ip.Results.largestOnly
            % get largest connected component
            mask = logical(getLargestCC(mask));
        end
        
        if ip.Results.TSOverlays
            roiYX = bwboundaries(mask);
            [ny,nx] = size(img); 
            setFigure(nx,ny,'off'); 
           
            imshow(-img,[]);
            hold on
            cellfun(@(x) plot(x(:,2),x(:,1),'r'),roiYX);
            saveas(gcf,[maskDir filesep 'overlays' filesep threshMethod forFilename 'Overlay' num2str(iMask, '%03d') '.tif' ]);
            close gcf
        end
               
        imwrite(mask,[ maskDir filesep 'mask' threshMethod forFilename num2str(iMask,'%03d') '.tif']);
        
        
       
    end % for iMask
    
    
end % for iProj

end