function plusTipBatchSubRoi(excludeRegion,micFromEdge,midPoint,minFrames)
% plusTipBatchSubRoi allows projects to be picked for sub-roi selection
%
% plusTipBatchSubRoi(excludeRegion,micFromEdge,midPoint,minFrames)
%
% this function allows the user to pick multiple projects from projList
% file(s) and divide each one into two regions: the central region and the
% peripheral region of the cell. you have to draw a region around the whole
% cell and can optionally select regions to exclude
%
% INPUT:
% excludeRegion: 1 to select one or more circular regions in which tracks
%                beginning within them should be excluded from analysis
% micFromEdge  : microns from the cell edge (user-chosen) for the boundary
% midPoint     : 1 if within-subRoi tracks must have half or more of their
%                lives within the region, 0 if they must have only some
%                number of frames within the region
% minFrames    : if midPoint=0, the number of frames to use as criterion
%                for whether a track is in the region



[projList]=combineProjListFiles;
if ~isempty(projList)
    % here we filter out any sub-directories
    a=struct2cell(projList);
    if isempty(strfind(a{1,1},'roi_'))
        a=a(2,:)';
    else
        a=a(1,:)';
    end
    a=sort(a);
    b=cellfun(@isempty, strfind(a,'sub'));
    a=a(b);
       
    % allow multiple projects to be selected
    [selection,selectionList]=listSelectGUI(a,[],'move',1);

    % if a project was selected, save projData info and get data
    if ~isempty(selection)
        projList=projList(selection,1);
    else
        projList=[];
    end
else
    msgbox('No projects selected.')
    projList=[];
end


fractionFromEdge=[];
savedROI=[];

for i=1:length(projList)
    projPath=formatPath(projList(i).anDir);
    p=load([projPath filesep 'meta' filesep 'projData.mat']);
    sourceProjData=p.projData;
    plusTipSubdivideRoi(sourceProjData,fractionFromEdge,savedROI,excludeRegion,micFromEdge,midPoint,minFrames);    
end