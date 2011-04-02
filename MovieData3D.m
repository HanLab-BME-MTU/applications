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
        
        %Constructor
        function obj = MovieData3D(channels,outputDirectory,...
                                    pixelSize,timeInterval,...
                                    zSpacing)


            if nargin < 1 || isempty(channels)
                channels = Channel3D;
            elseif ~isa(channels,'Channel3D')
                error('The first input must be a valid Channel3D object!');                
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
            superArgs{3} = outputDirectory;%I just force the moviedata to be in the output directory for simplicity...
            superArgs{4} = 'movieData.mat';            
            superArgs{5} = ''; %Notes not used yet...
            
            if nargin < 3 || isempty(pixelSize)

                pixelSize = str2double(inputdlg('Enter the XY pixel size in nm:'));

                if isnan(pixelSize) || isempty(pixelSize)
                    error('Invalid XY pixel size!')
                end                

            end

            superArgs{6} = pixelSize;

            if nargin < 4 || isempty(timeInterval)
                timeInterval = str2double(inputdlg('Enter the time interval in seconds:'));

                if isnan(timeInterval) || isempty(timeInterval)
                    error('Invalid time interval!')
                end
            end

            superArgs{7} = timeInterval;
            
            obj = obj@MovieData(superArgs{:});
                        

            if nargin < 5 || isempty(zSpacing)
                
                zSpacing = str2double(inputdlg('Enter the Z stack spacing in nm:'));

                if isnan(zSpacing) || isempty(zSpacing)
                    error('Invalid Z spacing!')
                end                

            end
            
            obj.zSpacing_ = zSpacing;
            
            %Set up the image size and number fields
            chanPaths = obj.getChannelPaths;
            nIm = cellfun(@(x)(numel(imDir(x))),chanPaths);
            if numel(unique(nIm)) > 1
                error('All channel directories must contain the same number of images!')
            end
            obj.nFrames_ = nIm(1);                

            %Check the image size and store it
            imNames = obj.getImageFileNames;
            imSize = zeros(nChan,3);
            for j = 1:nChan
                tmp = size(stackRead([chanPaths{j} filesep imNames{j}{1}])); %Just check the first image...
                if numel(tmp) ~= 3 || min(tmp) < 2
                    error(['The images in channel ' num2str(j) ' are not 3D stacks!'])
                end
                imSize(j,:) = tmp;                    
            end
            if size(unique(imSize,'rows'),1) > 1
                error('All channels must have images of the same dimension!')
            end

            obj.imSize_ = imSize(1,1:2);
            obj.nSlices_ = imSize(1,3);
            
            %Set the MovieData object as the owner of all the channels
            for j = 1:nChan
                obj.channels_(j).setOwner(obj);
            end
            %Save the new movieData to file
            obj.saveMovieData;

        end
        
        % Sanity check - verifies that all the movie information is
        % correct, and corrects the specified paths if the movie has been
        % moved
        function sanityCheck(obj, movieDataPath, movieDataFileName,askUser)
                        
           
            % Ask user by default for relocation
            if nargin < 4, askUser = true; end
            
            % Check if the path and filename stored in the movieData are the same
            % as the ones provided in argument. They can differ if the movieData
            % MAT file has been renamed, move or copy to another location.
            if nargin > 1
                
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
                        obj.relocateMovieData(movieDataPath); 
                    else
                        obj.setMovieDataPath(newMovieDataPath);
                    end
                end
            
                if  ~strcmp(obj.movieDataFileName_, movieDataFileName)
                    obj.movieDataFileName_ = movieDataFileName; 
                end
            
            end
            width = zeros(1, length(obj.channels_));
            height = zeros(1, length(obj.channels_));
            nFrames = zeros(1, length(obj.channels_));
            nSlices = zeros(1, length(obj.channels_));
            
            for i = 1: length(obj.channels_)
                [width(i) height(i) nFrames(i) nSlices(i)] = obj.channels_(i).sanityCheck;
            end
            
            assert(max(nFrames) == min(nFrames), ...
                'Different number of frames are detected in different channels. Please make sure all channels have same number of frames.')
            
            assert(max(width)==min(width) && max(height)==min(height), ...
                'Image sizes are inconsistent in different channels.\n\n')

            if numel(unique(nSlices)) > 1 
                error('Image z-slice numbers are inconsistent between channels!')
            end
		
            % Define imSize_ and nFrames_;     
            if ~isempty(obj.nFrames_)
                assert(obj.nFrames_ == nFrames(1), 'Record shows the number of frames has changed in this movie.')
            else
                obj.nFrames_ = nFrames(1);
            end
            
            if ~isempty(obj.imSize_)
                assert(obj.imSize_(2) == width(1) && obj.imSize_(1) ==height(1), 'Record shows image size has changed in this movie.')
            else
                obj.imSize_ = [height(1) width(1)];
            end
            
            if ~isempty(obj.nSlices_)
                assert(obj.nSlices_ == nSlices(1), 'Record shows image number of z-slices has changed in this movie.')
            else
                obj.nSlices_ = nSlices(1);
            end                
            

        end
        
        
        % ------ Set / Get Methods ----- %
        
        
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
                