function [mData]=generateMovieDatabase(projList)
% generateMoveieDatabase creates a dataset array for movies in projList

% INPUT : projList: list of projects created by getProj. If no input given,
%                   the user can select more than one projList file.
% OUTPUT: mData   : dataset array (see statistics toolbox help) containing
%                   the following columns for each movie: date, target,
%                   oligo number, movie number, roi number, sub-roi number, and
%                   project path. all these are extracted from the project
%                   path name, so the column names only make sense if the
%                   directory structure follows the pattern:
%                   ...\081120\Allstar\1027280\01\roi_1


% allow user to concatenate multiple project lists if no input given
if nargin<1
    [projList]=combineProjListFiles;
end

projPathParsed{size(projList,1),7}=[];
for iProj=1:size(projList,1)

    currentROI=formatPath(projList(iProj,1).anDir);
    % parse the path to get "words" used to identify target, oligo,
    % movie, and roi
    nChar=length(currentROI);
    if ispc
        filesepLoc=regexp(currentROI,'\\');
    else
        filesepLoc=regexp(currentROI,'\/');
    end
    wordStart=[1 filesepLoc+1]; wordEnd=[filesepLoc-1 nChar];
    words=cell(length(wordStart),1);
    for iWord=1:length(wordStart)
        words{iWord,1}=currentROI(wordStart(iWord):wordEnd(iWord));
    end

    % assign values to cell
    if strfind(words{end},'roi')
        projPathParsed{iProj,1}=words{end-4,1}; % date
        projPathParsed{iProj,2}=words{end-3,1}; % target
        projPathParsed{iProj,3}=words{end-2,1}; % oligo
        projPathParsed{iProj,4}=words{end-1,1}; % movie number
        projPathParsed{iProj,5}=words{end-0,1}; % roi number
    else
        projPathParsed{iProj,1}=words{end-6,1}; % date
        projPathParsed{iProj,2}=words{end-5,1}; % target
        projPathParsed{iProj,3}=words{end-4,1}; % oligo
        projPathParsed{iProj,4}=words{end-3,1}; % movie number
        projPathParsed{iProj,5}=words{end-2,1}; % roi number
        % end-1 is for subRoi header directory
        projPathParsed{iProj,6}=words{end-0,1}; % sub-roi number
    end
    projPathParsed{iProj,7}=currentROI;         % project directory path
end


nMovs=size(projPathParsed,1);
NameObs = strcat({'Movie'},num2str((1:nMovs)','%d'));

mData=dataset({projPathParsed(:,1:7),'date','target','oligo','movNum','roiNum','subNum','path'},'ObsNames',NameObs);

% make targets nominal and order labels so that Allstar is first and Negctr is second
mData.target=nominal(mData.target);
targetLabels=getlabels(mData.target);
nTargets=length(targetLabels);
% get indices of  Allstar and Negctr in the set of labels
allstarIdx=find(vertcat(cellfun(@(y) strcmp(y,'Allstar'),targetLabels)));
negctrIdx=find(vertcat(cellfun(@(y) strcmp(y,'Negctr'),targetLabels)));
% determine new order
if allstarIdx~=1 & negctrIdx~=2 & allstarIdx<negctrIdx
    newOrder=[allstarIdx negctrIdx 1:allstarIdx-1 allstarIdx+1:negctrIdx-1 negctrIdx+1:nTargets]';
else
    newOrder=[1:nTargets]';
end
% apply new order
mData.target=reorderlevels(mData.target,targetLabels(newOrder));
