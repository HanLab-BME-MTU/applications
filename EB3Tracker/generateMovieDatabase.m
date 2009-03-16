function [mData]=generateMovieDatabase(projList)

% allow user to concatenate multiple project lists if no input given
if nargin<1
    temp=[];
    userEntry='y';
    while strcmp(userEntry,'y')
        [fileName,pathName] = uigetfile('*.mat','Select projList.mat file');
        load([pathName filesep fileName]);

        temp=[temp; projList];

        disp(['Selected: ' pathName filesep fileName])
        userEntry = lower(input('Select another projList file? y/n ','s'));

    end
    clear projList;
    projList=temp;
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
        projPathParsed{iProj,6}=words{end-0,1}; % sub-roi number
    end
    projPathParsed{iProj,7}=currentROI;         % project directory path
end




nMovs=size(projPathParsed,1);
NameObs = strcat({'Movie'},num2str((1:nMovs)','%d'));

mData=dataset({projPathParsed(:,1:7),'date','target','oligo','movNum','roiNum','subNum','path'},'ObsNames',NameObs);

% make targets nominal and order labels so that Allstar is first and Negctr
% is second
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


%mData.date=ordinal(mData.date);
% dateLabels=getlabels(mData.date);
% 
%mData.oligo=nominal(mData.oligo);
% oligoLabels=getlabels(mData.oligo);
% 
%mData.movNum=nominal(mData.movNum);
% movieLabels=getlabels(mData.movNum);
% 
%mData.roiNum=nominal(mData.roiNum);







