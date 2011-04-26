function setup3DProjectMovieData(projectFolder,pixSizeXY,forceReplace)
%SETUPPROJECTMOVIEDATA sets up the MovieData3D for a collection of 3D movies
%
% setupProjectMovieData
% setupProjectMovieData(projectFolder)
% setupProjectMovieData(projectFolder,pixSizeXY)
% setupProjectMovieData(projectFolder,pixSizeXY,forceReplace)
%
% This goes through every movie folder in the specified project directory
% and creates and saves a MovieData3D object for each movie. The movies in
% the project directory must be set up so that each movie is in its own
% directory, and each movie directory has a sub-directory named "images"
% which contains the movie images. Movie folders can be set up in this way
% using setupProjectImageFolders.m and separateNumberedFiles.m
%   It is assumed that the images are stored as .TIF or .STK files with one
% file per timepoint. Any imaging parameters which can be read from the
% image headers will be stored in the movieData object.
%
% The analysis directory for each movie will be set to the movie directory,
% and the movieData object will be stored there.
%
% Input:
%   projectFolder - The parent directory containing all the movie folders.
%   This is the same folder you specified when using setupProjectImageFolders.m
%   Optional. If not input, the user will be asked to select a directory.
%
%   pixSizeXY - Physical size of a single image pixel - PRIOR TO BINNING -
%   in nanometers. This value will be used for all movies, but the binning
%   for each movie will be determined from the image header, and the pixel
%   size stored in the MovieData3D object will take the binning into
%   account.
%   Optional. If not input, no pixel size will be stored in the MovieData3D
%   objects.
%   
%   forceReplace - If true, a new movieData will be created even if there
%   is an existing one in the same directory. This will erase all
%   processing that has been logged in the movieData.
%   Optional. Default is false.
%
% 
% Output: 
% 
%   The newly created movieData structures will be saved in each movie's
%   analysis directory as a file named movieData.mat.
% 
% 
% Hunter Elliott 
% 4/2011
%

if nargin < 1 || isempty(projectFolder)
    projectFolder = uigetdir(pwd,'Select a project directory containing movie directories:');
    if projectFolder == 0 %If the user clicked cancel
        return
    end
end

if nargin < 2
    pixSizeXY = [];
end

if nargin < 3 || isempty(forceReplace)
    forceReplace = false;
end

%Get the folders for each movie
movieFolders = dir([projectFolder filesep]);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders)); 

   
nMovies = length(movieFolders);


%Go through each and try to set up the movieData
for j = 1:nMovies
            
    disp(['Setting up movie ' num2str(j) ' of ' num2str(nMovies)])           
    
    %Get current folder path for readability
    currDir = [projectFolder filesep movieFolders(j).name];
    
    clear chans;
    
    %Check for an existing movieData
    if forceReplace || ~exist([currDir filesep 'movieData.mat'],'file')
    
        %Look for sub-folder named "images"        
        if exist([currDir filesep 'images'],'dir')
            
                        
            %Check for sub directories containing the images for each
            %channel
            chanDir = dir([currDir filesep 'images']);
            chanDir = chanDir(arrayfun(@(x)(x.isdir && ... %Retain only the directories.
                            ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),chanDir));
            
            %Check if it is a single- or multiple-channel movie and set up the channels            
            if ~isempty(chanDir)
            
                %Determine the number of images and image size in each channel
                %directory.
                nChanDir = numel(chanDir);
                for i = 1:nChanDir
                    %Find images in this directory
                    chans(i) = Channel3D([currDir filesep 'images' filesep chanDir(i).name]); %#ok<AGROW> Can't initialize due to private fields                    
                end                                                
                                
            else %If it's a single-channel movie...
                chans = Channel3D([currDir filesep 'images']);                
            end                        
            
            try
                %try creating the moviedata. This will also call the sanity
                %check, so we catch any potential errors
                MD = MovieData3D(chans,currDir,pixSizeXY); %#ok<NASGU>
                disp('MovieData setup!')                                    
            catch errMess
                disp(['Problem setting up movie ' num2str(j) ' : ' errMess.message]);                    
            end
                            
        else
           disp(['Movie folder ' currDir ' has no sub-directory named "images" - cannot setup movieData!']) 
        end
    else
        disp(['Movie ' num2str(j) ' of ' num2str(nMovies) ' already had a movieData.mat! Doing nothing...'])
        
    end
end
