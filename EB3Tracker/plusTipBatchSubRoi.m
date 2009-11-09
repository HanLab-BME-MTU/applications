function plusTipBatchSubRoi(excludeRegion,micFromEdge,midPoint,minFrames)


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