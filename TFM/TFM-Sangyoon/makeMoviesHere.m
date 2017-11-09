function []=makeMoviesHere(analysisFolder)
%makeMoviesHere(analysisFolder) finds all folders starting with 'movie' and make movies of
% them in avi format (should be readable in mac too) with the same cropping
% information
% eg:
% analysisFolder='/project/bioinformatics/Danuser_lab/P01adhesion/analysis/Sangyoon/NA_RecruitmentProject/Alexia/2015-07-17/Vinculin5/Colocalization/analysis1'
% makeMoviesHere(analysisFolder)
% output movies will be stored in each movie folder

% Find folders containing 'movie'
if nargin<1
    analysisFolder=pwd;
end
folderNames=dir([analysisFolder filesep 'movies*']);
folderNames={folderNames.name}';

for kk=1:numel(folderNames)
    rootFolder = [analysisFolder filesep folderNames{kk}];
    subfolderNames=dir([rootFolder filesep 'tif*']);
    subfolderNames=subfolderNames(arrayfun(@(x) x.isdir,subfolderNames));
    subfolderNames={subfolderNames.name}';
    if kk==1
        rect=cell(numel(subfolderNames),1);
    end
    for jj=1:numel(subfolderNames)
        curFolder = [rootFolder filesep subfolderNames{jj}];
        if kk==1
            rect{jj,1} = makeMovieFromImageSequence(curFolder,30,true);
        else
            makeMovieFromImageSequence(curFolder,30,true,rect{jj,1});
        end
    end
end

% filnally same the same sized movie 

