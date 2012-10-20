classdef Channel3D < Channel
    
    %Class definition for a channel for a 3d movie
    
    properties(SetAccess = protected, GetAccess = public)
        
        %Z-Parameters (# of slices, slice spacing) are general to
        %MovieData, not channel, so are not stored here.
        
        
    end
    
    methods (Access = public)
        
        %Constructor  - same as Channel
        function obj = Channel3D(varargin)
                        
            obj = obj@Channel(varargin{:});                        
                        
        end
        
        % ------ Set / Get Methods ------- %
        
        
        %Overloads the default to make sure the owner is 3D
        function setOwner(obj,owner)
            %TEMP - check that owner is either the object owner or its
            %parent (for ROIs)
            if isa(owner,'MovieData3D')
                    obj.owner_ = owner;
            else
                error('The Channel can only be owned by a MovieData3D object!')
            end            
        end
        
        %---- Sanity Check ----%
        %Verifies that the channel specification is valid, and returns
        %properties of the channel
                
        %Overloads Channel.m sanity check to deal with 3d images
        function [width height nFrames nSlices ...
                  timeInterval binning zSpacing] = sanityCheck(obj)
        % Check the validity of each channel and return image info
            
            % Exception: channel path does not exist
            assert(logical(exist(obj.channelPath_, 'dir')), ...
                    'Channel paths is not a valid directory! Please double check the channel paths.')
            
            % Check the number of file extensions
            [fileNames nofExt] = imDir(obj.channelPath_,true);
            switch nofExt
                case 0
                    % Exception: No proper image files are detected
                    error('No proper image files are detected in:\n\n%s\n\nValid image file extension: tif, TIF, STK, bmp, BMP, jpg, JPG.',obj.channelPath_);
                                                  
                case 1
                    nFrames = length(fileNames); 
                    
                otherwise
                    % Exception: More than one type of image 
                    % files are in the current specific channel
                    error('More than one type of image files are found in:\n\n%s\n\nPlease make sure all images are of same type.', obj.channelPath_); 
            end
            
            % Check the consistency of image size in current channel
            imSize = zeros(nFrames,3);
            timeStamp = nan(nFrames,1);
            binning = zeros(nFrames,2);
            zSpacing = zeros(nFrames,1);
            for iFrame = 1:nFrames
                
                %Try to read the info from the header first, since it's much faster
                try
                    imInfo = stkinfo([obj.channelPath_ filesep fileNames(iFrame).name]);
                    %Extract the relevant parameters from the header.
                    imSize(iFrame,1) = imInfo.Height;
                    imSize(iFrame,2) = imInfo.Width;
                    imSize(iFrame,3) = imInfo.NumZPlanes;
                    zSpacing(iFrame) = imInfo.ZSpacing;
                    
                    %Get time stamp, converting to seconds and rounding -
                    %the data is stored only to 1-second precision and any
                    %other values are due to numerical error.
                    timeStamp(iFrame,:) = round(datenum(imInfo.DateTime,...
                                    'yyyy:mm:dd  HH:MM:SS') * 86400);
                    binning(iFrame,:) = imInfo.Binning;
                                        
                    
                catch em
                    disp(['Couldn''t read image header for file ' ...
                         [obj.channelPath_ filesep fileNames(iFrame).name] ...
                         ' : ' em.message ' loading entire image...']);
                    
                    %If we can't read the header, just go ahead and load
                    %the whole image
                    try
                        %Try stackRead first since it's slightly faster.
                        currIm = stackRead([obj.channelPath_ filesep fileNames(iFrame).name]);                        
                    catch %#ok<CTCH>
                        %If the compression format is unsupported,
                        %tif3dread may be able to open it.
                        currIm = tif3Dread([obj.channelPath_ filesep fileNames(iFrame).name]);
                    end
                    
                    if ndims(currIm) ~=3
                        error(['The image for frame ' num2str(iFrame) ' is not 3D!']);
                    end
    
                    imSize(iFrame,:) = size(currIm);
                    
                end
                
                
                
            end
            
            if size(unique(imSize,'rows'),1)>1
                error('All images must have the same width,heigth and number of z-slices! Check images!');            
            end
            
            if imSize(1,3) <= 1
                error('The images in the specified directory are not 3D!!')
            end
            
            if size(unique(binning,'rows'),1)>1
                error('All images must have the same pixel binning! Check images!');            
            end            
            
            if numel(unique(zSpacing))>1
                error('All images must have the same z-spacing! Check images!')
            end
            
            %Convert to desired output format
            width = imSize(1,2);
            height = imSize(1,1);
            nSlices = imSize(1,3);
            binning = binning(1,:);
            zSpacing = zSpacing(1) * 1e3;%Metamorph stores z-spacing in microns, so convert to nm
            
            %There are often small variations in the actual acquisition
            %times, so we process these carefully and check for odd frame
            %times
            if nFrames > 1
                dT = diff(timeStamp);
                uniqueDts = unique(dT);%Get the different intervals that were found
                nDts = numel(uniqueDts);
                nOfEach = arrayfun(@(x)(nnz(dT == x)),uniqueDts);%See how many frames had each interval
                [~,iMostCommon] = max(nOfEach);%Find the most common interval
                timeInterval = uniqueDts(iMostCommon);

                if nDts > 0
                    %If there were some frames with different time intervals,
                    %warn the user.                
                    warning('Channel3D:OutlierTimeStamp',...
                        'There were %d outlier time stamp(s) out of %d images! \n For the images in %s \n the most common time interval was %s seconds, \n while the outlier time intervals were different by as much as %s second(s)!',...
                        nDts,nFrames,obj.channelPath_,num2str(timeInterval),...
                        num2str(max(abs(timeInterval-dT))));%We use num2str because the number display used within warning is shitty                    
                end 
                    
            else
                timeInterval = NaN;
            end            
        end
                        
        
    end    
    
    
end

            