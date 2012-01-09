classdef MovieData3D < MovieData
    %Movie management class for 3D movies where each timepoint is stored as
    %a image stack
    %
    %Hunter Elliott,
    %6/2010
    %
    properties (SetAccess = protected, GetAccess = public)
        
        %Image Parameters
                
        nSlices_                % Number of Z-Slices in a stack.
        zSpacing_               % Spacing between z-slices, nm
        
    
    end
    
                      
    methods (Access = public)
        
        
        function obj = MovieData3D(channels,outputDirectory,...
                                    preBinPixSizeXY)
            %MOVIEDATA3D constructor for a MovieData object describing a 3D movie
            %
            % movieData3D = MovieData3D
            % movieData3D = MovieData3D(channels,outputDir,preBinPixSizeXY)
            %
            % This function creates a MovieData3D object describing a 3D movie
            % with a single .STK file per time-point, with each channel stored
            % in a separate directory. 
            %
            % **NOTE*** After the MovieData3D object is created, the
            % sanityCheck method will be called. This method will perform two
            % actions: First it will determine the binning, number of z planes,
            % time interval and spacing between z-planes by reading the image
            % headers, and second it will verify that all the information in
            % the movieData3D is valid and in accordance across channels.
            %
            % Input:
            %
            %   channels - a Channel3D object or an array of channel3D objects
            %   describing the image channels of the movie.
            %   Optional. If not input, the user will be asked to select
            %   directories containing images.
            %
            %   outputDir - The directory wher the output of any processing
            %   performed on the movie will be stored, and the location where
            %   the movieData will be saved to disk.
            %   Optional. If not input, the user will be asked to select a
            %   directory.
            %
            %   preBinPixSizeXY - The image pixel size PRIOR TO BINNING. The
            %   binning will be read from the image header and the final pixel
            %   size will take the binning into account.
            %   Optional. If not input, the user will be asked to enter it.
            %
            % Output:
            %
            % movieData3D - The MovieData3D object describing the movie. The
            % object will also be saved to disk in the movie's analysis
            % directory.
            %
            % Hunter Elliott
            % 4/2011
            %
            
            mdFileName = 'movieData.mat'; %Name for saving MovieData3d object to disk
                                
            if nargin < 1 || isempty(channels)
                channels = Channel3D;
            elseif ~isa(channels,'Channel3D')
                error('The first input must be a valid Channel3D object or array of Channel3D objects!');                
            end

            nChan = numel(channels);

            for j = 1:nChan

                if isempty(channels(j).channelPath_);
                    tmp = uigetdir(pwd,['Select the directory with stacks for channel ' num2str(j) ':']);
                    if tmp == 0
                        error('You must specify a directory to continue!')
                    else
                        channels(j).setChannelPath(tmp);                        
                    end                    
                end
                if ~exist(channels(j).channelPath_,'dir') || isempty(imDir(channels(j).channelPath_))
                    error(['The directory specified for channel ' num2str(j) ' is not a valid directory containing image stacks!'])
                end

            end

            superArgs{1} = channels;

            if nargin < 2 || isempty(outputDirectory);
                outputDirectory = uigetdir(pwd,'Select a directory to store the output:');

                if outputDirectory == 0
                    error('You must specify an output directory!')
                end

            end
            if ~exist(outputDirectory,'dir')
                error('Invalid output directory!')
            end
            
            superArgs{2} = outputDirectory;
            superArgs(3:4) = {'movieDataPath_',outputDirectory};%I just force the moviedata to be in the output directory for simplicity...
            superArgs(5:6) = {'movieDataFileName_',mdFileName};                        
            
            if nargin < 3 || isempty(preBinPixSizeXY)

                preBinPixSizeXY = str2double(inputdlg('Enter the PRE-BINNING, XY image pixel size in nm:'));

                if isnan(preBinPixSizeXY) || isempty(preBinPixSizeXY)
                    error('Invalid XY pixel size!')
                end                

            end
            %We don't pass the pixel size to the MovieData constructor
            %because it is pre-binning and we need to convert it during the
            %sanity check.                                    
            
            obj = obj@MovieData(superArgs{:});
                                    
            %Set the MovieData object as the owner of all the channels
            for j = 1:nChan
                obj.channels_(j).setOwner(obj);
            end
            
            %Call the sanity check to set parameters and check images etc.
            obj.sanityCheck([],[],0,preBinPixSizeXY)
            
            
            %Save the new movieData to file
            obj.save;

        end
        
       
        function sanityCheck(obj, movieDataPath, movieDataFileName,askUser,...
                             preBinPixSizeXY)
        % Sanity check - verifies that all the movie information is
        % correct, corrects the specified paths if the movie has been
        % moved, and stores imaging parameters if not already present.     
                       
            % Ask user by default for relocation
            if nargin < 4, askUser = true; end
            if nargin < 5, preBinPixSizeXY = []; end
            
            % Check if the path and filename stored in the movieData are the same
            % as the ones provided in argument. They can differ if the movieData
            % MAT file has been renamed, move or copy to another location.
            if nargin > 1 && ~isempty(movieDataPath) && ~isempty(movieDataFileName)
                
                %Remove ending file separators if any
                endingFilesepToken = [regexptranslate('escape',filesep) '$'];
                path1 = regexprep(obj.movieDataPath_,endingFilesepToken,'');
                path2 = regexprep(movieDataPath,endingFilesepToken,'');
                if  ~strcmp(path1, path2)
                    
                    if askUser
                        relocateMsg=sprintf(['The movie data located in \n%s\n has been relocated to \n%s.\n'...
                            'Should I try to relocate the components of the movie data as well?'],path1,path2);
                        confirmRelocate = questdlg(relocateMsg,'Movie Data','Yes','No','Yes');
                    else
                        confirmRelocate = 'Yes';
                    end
                    
                    if strcmp(confirmRelocate,'Yes')
                        obj.relocate(movieDataPath); 
                    else
                        obj.setMovieDataPath(newMovieDataPath);
                    end
                end
            
                if  ~strcmp(obj.movieDataFileName_, movieDataFileName)
                    obj.movieDataFileName_ = movieDataFileName; 
                end
            
            end
            nChan = numel(obj.channels_);
            width = zeros(1,nChan);
            height = zeros(1,nChan);
            nFrames = zeros(1,nChan);
            nSlices = zeros(1,nChan);
            timeInts = zeros(1,nChan);
            binning = zeros(2,nChan);
            zSpacing = zeros(1,nChan);
            
            %Call the individual channel sanity checks, and get the image
            %parameters for each channel
            for i = 1:nChan
                [width(i) height(i) nFrames(i) nSlices(i) ...
                    timeInts(i) binning(:,i) zSpacing(i)] = ...
                                            obj.channels_(i).sanityCheck;
            end            
            
            %Go through all the parameters and make sure they are the same
            %for each channel, and that they agree with what is stored in
            %the MovieData3D object.
            
            assert(max(width)==min(width) && max(height)==min(height), ...
                'Image sizes are inconsistent in different channels.\n\n')
            
            if ~isempty(obj.imSize_)
                assert(obj.imSize_(2) == width(1) && obj.imSize_(1) ==height(1), 'Record shows image size has changed in this movie.')
            else
                obj.imSize_ = [height(1) width(1)];
            end            

            assert(max(nFrames) == min(nFrames), ...
                'Different number of frames are detected in different channels. Please make sure all channels have same number of frames.')            
                        
            if ~isempty(obj.nFrames_)
                assert(obj.nFrames_ == nFrames(1), 'Record shows the number of frames has changed in this movie.')
            else
                obj.nFrames_ = nFrames(1);
            end
            
            if numel(unique(nSlices)) > 1 
                error('Image z-slice numbers are inconsistent between channels!')
            end            
            if ~isempty(obj.nSlices_)
                assert(obj.nSlices_ == nSlices(1), 'Record shows image number of z-slices has changed in this movie.')
            else
                obj.nSlices_ = nSlices(1);
            end
            
            if obj.nFrames_ > 1 && numel(unique(timeInts)) > 1
                error('The time intervals differ between the channels!')
            end
            if ~isempty(obj.timeInterval_)
                if obj.timeInterval_ ~= timeInts(1) && obj.nFrames_ > 1
                    error('The time interval specified for this movie does not agree with the image headers!')
                end
            elseif obj.nFrames_ > 1
                obj.timeInterval_ = timeInts(1);
            end                                        
            if size(unique(binning),2) > 1
                error('The binning differs between channels!')
            end
            if ~isempty(obj.binning_)
                assert(all(obj.binning_ == binning(:,1)'),'Record shows that the binning has changed on the images in this movie!')
            else
                obj.binning_ = binning(:,1)';
            end
            if numel(unique(zSpacing)) > 1
                error('The z-spacing differs between channels of the movie!')
            end
            if ~isempty(obj.zSpacing_)
                assert(obj.zSpacing_ == zSpacing(1),'Record shows that the z-spacing has changed in this movie!')
            else
                obj.zSpacing_ = zSpacing(1);
            end
            
            %Finally, set or check the pixel size taking the binning into
            %account            
            if ~isempty(preBinPixSizeXY)                
                if ~isempty(obj.pixelSize_)
                    assert(obj.pixelSize_ == (preBinPixSizeXY*obj.binning_(1)),...
                        'Stored post-binning pixel size disagrees with input pre-binned pixel size and image binning!')
                else
                    obj.pixelSize_ = preBinPixSizeXY*obj.binning_(1);
                end                                
            end

        end
        
        
        % ------ Set / Get Methods ----- %
        
        %At some point I should convert these to method used in MovieData,
        %where calling a specific set method is not required.
        function setzSpacing(obj,zSpacing)            
            if ~isempty(obj.zSpacing_)
                error('Z Spacing has already been set and cannot be changed!')
            end            
            if nargin ~= 2 || isempty(zSpacing) || ...
                    zSpacing <= 0 || numel(zSpacing) ~= 1 || ~isreal(zSpacing) %you never know what they may try ;)
                error('Invalid z-spacing specification. The value must be a positive scalar specifying the z-spacing in nm!')
            end            
            obj.zSpacing_ = zSpacing;
        end
        
    end
end
                