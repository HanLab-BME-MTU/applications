function [] = mitoticPerformManualFunctions(ML,varargin)
% mitoticPerformManualFunctions: simple function to make a spindleMask, 
% manually document the pole locations, and/or check segmentation masks. 

%% 
% INPUT: 
% movieList: (REQUIRED) : movieList object 
%
% inFilename: (PARAM) : character array
%                        Default: 'masks'
%                        By default looks in the 
%                        [MD.outputDirectory_ filesep 'masks'];
%                        for masks to load 
%
% spindleMask: (PARAM) : logical 
%                       Default : true 
% 
% documentPoles : (PARAM) : logical  
%                        Default : true 
%
%
% checkCellEdgeSeg : (PARAM) : logical 
%                     Default : true 
%% Check Input 
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('inFilename','masks',@(x) ischar(x));

ip.addParameter('spindleMask',true);

ip.addParameter('documentPoles',true);

ip.addParameter('checkCellEdgeSeg',true); %

ip.parse(varargin{:});
%% Set Up 
if isa(ML,'MovieData')
    projList{1,1} = [ML.movieDataPath_ filesep ML.movieDataFileName_];
else
    projList = ML.movieDataFile_;
end
%% Start Project Wrapper 
for iProj = 1:numel(projList)
    
    % Set UP
    load(projList{iProj});
    
    maskDirIn = [MD.outputDirectory_ filesep ip.Results.inFilename];
    maskDirOut = [MD.outputDirectory_];
    maskInternal = zeros(MD.imSize_);
    mask = zeros(MD.imSize_);
    
    img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{1}])) ;
    
    img = img./((2^MD.camBitdepth_)-1);
    normImg = (img-min(img(:)))./(max(img(:)) - min(img(:)));
    %%  Make a spindle mask
    if ip.Results.spindleMask
        
        if exist([MD.outputDirectory_ filesep 'roiMaskSpindle.tif'],'file');
            maskInternal =  logical(imread(([MD.outputDirectory_ filesep 'roiMaskSpindle.tif'])));
            display(['Spindle Mask Found for ' projList{iProj} ': Loading']); 
        else
            
            h = msgbox('please make a spindle roiMask');
            uiwait(h);
            figure;
            
            [maskInternal, polyXcoord,polyYcoord] = roipoly(normImg);
            
            imwrite(maskInternal,[maskDirOut filesep 'roiMaskSpindle.tif'],'tif');
            close(gcf)
        end
    end % spindle mask
    %% Test for bad masks
    if ip.Results.checkCellEdgeSeg
        
        if exist([MD.outputDirectory_ filesep 'badMaskFlag.mat'],'file');
            display([projList{iProj} 'cell edge mask has already been checked : Skipping']);
        else 
            [listOfMasks] = mitoticSearchFiles('.tif',[],maskDirIn,0,'all',1);
            maskExternal = logical(imread(listOfMasks{1}));
            mask = maskExternal - maskInternal;
            
            figure;
            imshow(img,[]);
            hold on
            roiYX = bwboundaries(mask);
            cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
            reply1 = questdlg('Is this a good mask?') ;
            if strcmpi(reply1,'no')
                
                badMaskFlag = 0; %
            elseif strcmpi(reply1,'yes')
                badMaskFlag = 1;
            end
            
            save([MD.outputDirectory_ filesep 'badMaskFlag.mat'],'badMaskFlag'); %
        
        end
    end % checkCellEdgeSeg
%% start pole documentation 

if ip.Results.documentPoles
    
    if exist([MD.outputDirectory_ filesep 'poleInfo.mat'],'file'); 
           display(['PoleInfo for ' projList{iProj} ' has been found  : Skipping']);
    else 
    poleMaskTotal = zeros(size(img));
    clickPole = 1;
    ptCount =1;
    while clickPole == 1
        
        imshow(img,[]);
        hold on
        roiYX = bwboundaries(mask);
        cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
        spy(poleMaskTotal,'r',50); % plot the poles thus far ;
        reply2 = questdlg('Document New Pole?');
        
        
        if strcmpi(reply2,'yes')
            %         imshow(img,[])
            %         hold on
            %         cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX);
            %         spy(poleMaskTotal,'r',50); % plot the poles thus far ;
            h=impoint;
            position = wait(h);
            idx = sub2ind(size(img),round(position(2)),round(position(1)));
            coords(ptCount,2) = position(2);
            coords(ptCount,1) = position(1);
            temp = zeros(size(img));
            temp(idx) = 1;
            poleMasksInd(:,:,ptCount) = temp ;
            poleMaskTotal(idx) = 1; % add points to mask
            ptCount = ptCount+1;
        elseif (strcmpi(reply2,'no') && ptCount ==1)
            coords(ptCount,2) = nan;
            coords(ptCount,1) = nan;
            poleInfo.coords = coords;
            poleInfo.numPoles = 0 ;
            save([projList(iProj).anDir filesep 'poleInfo.mat'],'poleInfo');
            clickPole = 0;
        else
            % finish up and save
            if ~isempty(coords)
                poleInfo.coords = coords; % save coordinations of  poles
                poleInfo.poleMasksInd = poleMasksInd;
                poleInfo.poleMaskTotal = poleMaskTotal;
                numPoles = length(coords(:,1));
                poleInfo.numPoles = numPoles; % save info regarding the number of poles
                if numPoles > 2
                    poleInfo.mult = 1;
                else poleInfo.mult =0;
                end
                
                save([MD.outputDirectory_ filesep 'poleInfo.mat'],'poleInfo');
                
            end % isemtpy coords
            clear coords poleMasksInd poleMaskTotal poleInfo
            
            clickPole = 0; % stop iteration
        end % if strcmpi
        
    end % while
    
    close(gcf)
    end  % if exist 
end  % if documentPoles
end  % for iProj
         
         
         
         
end

     
         
        

 


    



