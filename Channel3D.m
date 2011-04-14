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
            if isempty(obj.owner_)
                if isa(owner,'MovieData3D')
                    obj.owner_ = owner;
                else
                    error('The Channel can only be owned by a MovieData3D object!')
                end

            else
                error('This channel already has an owner, and this property cannot be changed!');
            end
        end
        
        %---- Sanity Check ----%
        %Verifies that the channel specification is valid, and returns
        %properties of the channel
                
        %Overloads Channel.m sanity check to deal with 3d images
        function [width height nFrames nSlices] = sanityCheck(obj)
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
            for iFrame = 1:nFrames
                %TEMP - Go through each frame, and load the image to check
                %size (imfinfo.m does not return the number of
                %z-slices for metamorph STK files). Hopefully at some point
                %I'll figure out how to do this without loading the
                %images.-HLE
                
                %Try stackread first since it's faster
                try
                    currIm = stackRead([obj.channelPath_ filesep fileNames(iFrame).name]);
                catch %#ok<CTCH>
                    %It it's an un-supported compression format, try
                    %tif3dread.m also
                    currIm = tif3Dread([obj.channelPath_ filesep fileNames(iFrame).name]);
                end
                
                if ndims(currIm) ~=3
                    error(['The image for frame ' num2str(iFrame) ' is not 3D!']);
                end
                    
                imSize(iFrame,:) = size(currIm);
                
            end
            
            if size(unique(imSize,'rows'),1)>1
                error('All images must have the same width,heigth and number of z-slices! Check images!');            
            end
            
            %Convert to desired output format
            width = imSize(1,2);
            height = imSize(1,1);
            nSlices = imSize(1,3);
                                            
        end
                        
        
    end    
    
    
end

            