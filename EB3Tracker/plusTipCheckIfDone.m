function [allProjects,notDone]=plusTipCheckIfDone(fileType,chckFrNum)
% plusTipCheckIfDone checks whether detection or post-processing has been finished for a group of movies
%
% SYNOPSIS: [allProjects,notDone]=plusTipCheckIfDone(fileType)
%
% INPUT: fileType  : 1 to check if detection is done 
%                    2 to check if tracking is done
%                    3 to check if post-processing is done (default)
%          user will be asked to select one or more projList.mat files
%                    containing the list(s) of movies to check  (see getProj
%                    for how to create projList files)
%        chckFrNum:  1 to add how many frames there are in each movie (time
%                    consuming)
%
% OUTPUT: allProjects: nMovies x 4 matrix where column 1 contains the
%                      directory name and columns 2:4 contain the
%                      timestamps of creation of the movieInfo.mat file, the
%                      trackResults.mat file, and the projData.mat file, if
%                      they exist. if chckFrNum=1, there will be a 5th
%                      column with the number of frames for each movie.
%         notDone    : vector of movie numbers corresponding to projList for
%                      which the querried file does not exist
%             


allProjects=[];
notDone=[];

if nargin < 1
    fileType=3;
end
if nargin < 2
    chckFrNum=0;
end

% allow user to concatenate multiple project lists
[projList]=combineProjListFiles;
if isempty(projList)
    error('No projects selected.')
end

% get list of all projects
[allProjList]=projList2Mat(projList);

% list of all the files in each FEAT directory
dirContents=cellfun(@(i) dir([i filesep 'feat']),allProjList,'uniformoutput',0);
% convert each structure to a cell array
dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
% get index for movieInfo file, if it exists
movieInfoLoc=cellfun(@(i) find(strcmpi(i(1,:),'movieInfo.mat')),dirContentsCell,'uniformoutput',0);
% get timestamp for each movieInfo file
dates1=cellfun(@(i,j) i(2,j),dirContentsCell,movieInfoLoc,'uniformoutput',0);


% list of all the files in each TRACK directory
dirContents=cellfun(@(i) dir([i filesep 'track']),allProjList,'uniformoutput',0);
% convert each structure to a cell array
dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
% get index for projData file, if it exists
trackResultsLoc=cellfun(@(i) find(strcmpi(i(1,:),'trackResults.mat')),dirContentsCell,'uniformoutput',0);
% get timestamp for each projList file
dates2=cellfun(@(i,j) i(2,j),dirContentsCell,trackResultsLoc,'uniformoutput',0);

% list of all the files in each META directory
dirContents=cellfun(@(i) dir([i filesep 'meta']),allProjList,'uniformoutput',0);
% convert each structure to a cell array
dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
% get index for projData file, if it exists
projListLoc=cellfun(@(i) find(strcmpi(i(1,:),'projData.mat')),dirContentsCell,'uniformoutput',0);
% get timestamp for each projList file
dates3=cellfun(@(i,j) i(2,j),dirContentsCell,projListLoc,'uniformoutput',0);

dates=[dates1 dates2 dates3];

% find out if any of them don't have dates
notDone=find(cell2mat(cellfun(@(i) isempty(i),eval(['dates' num2str(fileType)]),'uniformoutput',0)));

if chckFrNum==1
    % get how many images are in each project
    homeDir=pwd;
    nImages=cell(length(allProjList),1);
    for iProj=1:length(allProjList)
        cd(allProjList{iProj,1});
        cd ..
        cd('images')
        nImages{iProj,1}=size(searchFiles('.tif',[],pwd,0),1);
    end
    cd(homeDir)
    allProjects=[allProjList dates nImages];
else
    allProjects=[allProjList dates];
end

