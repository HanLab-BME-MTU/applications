function [outputData3S,outputDataMP] = threeVsTwoD(idlist, dataProperties)
%function to compare 3D motion analysis with 2D data

%----Test input

if nargin < 1 || isempty(idlist)
    inputData = loadIdlistsFromFileList;
else
    inputData(1).idlist = idlist;
    inputData(1).dataProperties = dataProperties;
end

%----End test input

%---- calc new idlist and anaDat

for i = 1:length(inputData)%loop through all data
    
    %get idlist & dataProperties
    idlist = inputData(i).idlist;
    dataProperties = inputData(i).dataProperties;
    
    %calculate depth-of-field
    pixel_Z = dataProperties.PIXELSIZE_Z;
    wvl     = dataProperties.WVL;
    nRef    = 1.33;
    NA      = dataProperties.NA;
    nZSlices = dataProperties.movieSize(3);
    
    dofInterval = 1.1*wvl*nRef/(NA^2); %we accept a signal if it is at 50% intensity in the focal plane (ca. 1.1 dof)
    
    
    %loop through idlist to
    % 3S) find in which interval the tags are
    % MP) set the z-coord to zero
    
    minMaxSliceNum = zeros(length(idlist),2);
    idlistMP = idlist;
    
    for t = 1:length(idlist)
        if ~isempty(idlist(t).linklist)
            
            %do the easy one first: set idlistMP.linklist(z-coord) to zero
            idlistMP(t).linklist(:,11) = 0;
            
            %find min and max z-distance of the tags
            tagPos= (idlist(t).linklist(:,11));
            %store
            minMaxTagPos(t,:) = [min(tagPos),max(tagPos)];
        end
    end
    
    %now check in which frames the tags could be in focus
    tagInterval = diff(minMaxTagPos,1,2);
    badIntervalIdx = find(tagInterval>dofInterval | tagInterval == 0);
    goodIntervalIdx = find(tagInterval<=dofInterval & tagInterval ~= 0);
    
    minMaxTagPos(badIntervalIdx,:) = 0;
    
    % according to Chad, it is possible to adjust the focal plane in 0.5
    % seconds. So we just assume there is a perfect microscopist, and keep
    % all inFocus
    outOfFocusFrames = badIntervalIdx;
    if ~isempty(goodIntervalIdx) | length(goodIntervalIdx)==1
        do3S = 1;
    else
        do3S = 0;
    end
    
%     %for all nonzero minmax: do a continuous histogram to find the best
%     %focal plane
%     %unsing non-normalized plot with interval-length leads at least in one
%     %case to a loner anaDat
%     if ~isempty(goodIntervalIdx) | length(goodIntervalIdx)==1
%         [yHist,xHist] = contHisto([mean(minMaxTagPos(goodIntervalIdx,:),2),diff(minMaxTagPos(goodIntervalIdx,:),1,2)],'norm',0,0);
%         %[yHist,xHist] = contHisto([mean(minMaxTagPos(goodIntervalIdx,:),2),repmat(dofInterval/2,size(goodIntervalIdx))],'norm',1,1);
%         
%         [dummy,yMaxIdx] = max(yHist);
%         focalPlane = xHist(yMaxIdx);
%         
%         %if the focus region would lead us outside of the movie, adjust
%         if focalPlane-dofInterval/2 < 0
%             focalPlane = dofInterval/2;
%         end
%         if focalPlane+dofInterval/2 > nZSlices*pixel_Z
%             focalPlane = nZSlices*pixel_Z-dofInterval/2;
%         end
%         
%         
%         %find all the not ok frames (= find ~okframes)
%         goodFramesList = (minMaxTagPos(:,1) >= focalPlane-dofInterval/2 & minMaxTagPos(:,2) <= focalPlane+dofInterval/2);
%         inFocusFrames = find(goodFramesList);
%         outOfFocusFrames = find(~goodFramesList);
%         do3S = 1;
%         if isempty(inFocusFrames) | length(inFocusFrames) == 0
%             do3S = 0;
%         end
%     else
%         outOfFocusFrames = [1:length(idlist)];
%         do3S = 0;
%     end
    
    idlist3S = idlistMP;
    %remove all frames where tags are out of focus
    [idlist3S(outOfFocusFrames).linklist] = deal([]);
    
    %calc anaDat
    anaDatMP = adgui_calc_anaDat(idlistMP, dataProperties,'idlisttrack_L');
    if do3S
        anaDat3S = adgui_calc_anaDat(idlist3S, dataProperties,'idlisttrack_L');
        
        %assign outputData
        outputData3S(i).idlist = idlist3S;
        outputData3S(i).dataProperties = dataProperties;
        outputData3S(i).anaDat = anaDat3S;
        outputData3S(i).frameList = inFocusFrames;
        outputData3S(i).focalPlane = focalPlane;
    else
        outputData3S(i).idlist = [];
        outputData3S(i).dataProperties = [];
        outputData3S(i).anaDat = [];
        outputData3S(i).frameList = [];
        outputData3S(i).focalPlane = [];
    end
    
    outputDataMP(i).idlist = idlistMP;
    outputDataMP(i).dataProperties = dataProperties;
    outputDataMP(i).anaDat = anaDatMP;
end

%---- end calc new idlist and anaDat