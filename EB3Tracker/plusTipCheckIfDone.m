function [allProjects,notDone]=plusTipCheckIfDone(fileType)
% plusTipCheckIfDone checks whether detection or post-processing has been finished for a group of movies
%
% SYNOPSIS: [allProjects,notDone]=plusTipCheckIfDone(fileType)
%
% INPUT: fileType: 1 to check if detection is done 
%                  2 to check if tracking is done
%                  3 to check if post-processing is done (default)
%        user will be asked to select one or more projList.mat files
%                  containing the list(s) of movies to check  (see getProj
%                  for how to create projList files)
%
% OUTPUT: allProjects: nMovies x 2 matrix where column 1 contains the
%                      directory name and column 2 contains the timestamp of
%                      creation for either the movieInfo.mat file (if fileType
%                      is 1) or the projData.mat file (if fileType is 2)
%         notDone    : vector of movie numbers corresponding to projList for
%                      which the querried file does not exist
%             


allProjects=[];
notDone=[];

if nargin < 1
    fileType=3;
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

allProjects=[allProjList dates];