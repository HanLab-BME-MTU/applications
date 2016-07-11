function [ groupList ] = mitoticCreateGroupListFromMovieLists(varargin)
% createGroupListFromMovieLists : small helper to create a groupList
% formats from a set of movieLists. 
%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;
%ip.KeepUnmatched = true;

ip.addOptional('movieListCell',[]); 
ip.addOptional('groupNames',[]);
ip.addParameter('outputDirectory',[]); 
ip.addParameter('groupListName','groupList'); 


ip.addParameter('makeSubRoiGroupList',true); 
% 
ip.addParameter('subRoiFilename',[], @(x) ischar(x) || isempty(x));
ip.addParameter('bipolar',true,@(x) islogical(x)); % 
ip.addParameter('sortGroupList',true, @(x) islogical(x)); 
ip.addParameter('subRoiGroupListFilename',[]); 
ip.addParameter('replaceStr',[]); 



ip.parse(varargin{:});
p = ip.Results; 
%%
%% Check input 
if isempty(ip.Results.movieListCell)
 selectMore = 'yes';
 count = 1;
 while strcmpi(selectMore,'yes'); 
    [filename,pathname]   = uigetfile(pwd,'Choose a MovieList File');
    % collect MD paths. 
    load([pathname filesep filename]); 
    list = ML.movieDataFile_; 
    projList(:,2) = cellfun(@(x) [mitoticUpDirectory(x,1) filesep 'TrackingPackage'],list,'uniformoutput',0); 
    grpName{1} = strrep(ML.movieListFileName_,'.mat',''); 
    MDs = load(cellfun(@(x) x, ML.movieDataFile_,'uniformoutput',0));
    projList(:,3) = cellfun(@(x) x.getChannelPaths{1},MDs,'uniformoutput',0); 
    nMovies = numel(list); 
    projList(:,1) = repmat(grpName,[nMovies,1]); 
    groupList{count} = projList; 
    clear projList
    count = count+1; 
    selectMore = questdlg('Get More Movie Lists?') ; 
    
 end 
else 
    for i= 1:numel(ip.Results.movieListCell)
        ML = ip.Results.movieListCell{i};
        list = ML.movieDataFile_;
        projList(:,2) = cellfun(@(x) [mitoticUpDirectory(x,1) filesep 'TrackingPackage'],list,'uniformoutput',0);
        grpName{1} = strrep(ML.movieListFileName_,'.mat','');
        MDs = cellfun(@(x) load(x),ML.movieDataFile_,'uniformoutput',0);
        projList(:,3) = cellfun(@(x) x.MD.getChannelPaths{1},MDs,'uniformoutput',0); 
        nMovies = numel(list);
        projList(:,1) = repmat(grpName,[nMovies,1]);
        groupList{i} = projList;
        clear projList
       
    end
    
end 

groupList = vertcat(groupList{:});
if ~isempty(ip.Results.outputDirectory)
    save([ip.Results.outputDirectory filesep ip.Results.groupListName],'groupList');
end

if ip.Results.makeSubRoiGroupList
   
   groupList = mitoticMakeSubRoiGroupList(groupList,p);
         
end 


end

