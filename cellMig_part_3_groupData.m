function [groupData]=cellMig_part_3_groupData(groupData,folderToAdd,yModu_kPa,cc,cond,justPlot)
clear persistent;

if nargin<1 || isempty(groupData)
    try
        load('xGroupedData.mat');
        numDataSets = length(groupData.dataSetList);
        numClasses  = length(groupData.kPaClass);
    catch exception
        display('Couldnt find grouped data assume that it is the first data set');
        numDataSets=0;
        numClasses=0;
        groupData=[];
    end
else
    numDataSets = length(groupData.dataSetList);
    numClasses  = length(groupData.kPaClass);
end

try
    numNotUsed = length(groupData.notUsed);
catch exception
    numNotUsed = 0;
end

try
    load([folderToAdd,filesep,'xResults.mat']);
catch exception
    display(['Couldnt add results for folder: ',folderToAdd,'! xResult.mat not found!']);
    return;
end
    

imageFileListNuclei=[folderToAdd,filesep,'Nuclei'];
try
    imageFileListNuclei=getFileListFromFolder(imageFileListNuclei);
catch exception
    [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
        'Select First nuclei Image');
    
    if ~ischar(filename) || ~ischar(pathname)
        return;
    end
    
    imageFileListNuclei = getFileStackNames([pathname filesep filename]);
end
[~,fnameFirstNucleiImg]=getFilenameBody(imageFileListNuclei{1});

if nargin<6 || isempty(justPlot)
    justPlot=0;
end

% 'xResults.mat','signFac','cell2edgeDist','tracksMatxCord','tracksMatyCord
% 'velMatxCord','velMatyCord','velMatMag','densityMeasurement','coveredDist
% 'roughness','freqSpec'
% minTrackLength','timeWindow','dt','pixSize_um','toDoList','badFlag'

dFramesCorr=24;      % The correlation analysis is done only every dFramesCorr-frames
dR=100;             % The disc width for calculating I(r+dr);
minlength=12; %12;       % The minimal movie length in hours, shorter movies are disregarded!           
timeWindowCorr=6    % Time window for averaging velocities for correlation analysis.
bandWidth=250;      % Band width of strips with iso-distance to the wound edge for correlation analysis. binPix=100 used for the density seems to be too small!
maxDist=max(cell2edgeDist(:,end)); % Bands with a larger distance are to the wound edge are disregarded in the correlation analysis.
%maxDist=1000;

binPix =densityMeasurement.binPix; % Band width of strips with iso-distance to the wound edge for density analysis.
numFrames  = toDoList(end);  % Number for good frames in the movie.

if ~justPlot
    %**********************************************************************
    % Check if data set is already grouped
    %**********************************************************************
    if numDataSets>0 && sum(strcmp(fnameFirstNucleiImg,groupData.dataSetList))>0
        display('The data set is already in the group list, nothing to do!')
        return;
    end
    
    %**********************************************************************
    % Check if data set is useful (it might be too short)
    %**********************************************************************
    if toDoList(end)*dt<minlength
        % If the sheet has been tracked less than 12h:
        % Disregard it:
        display(['Not using data set (too short): ',folderToAdd])
        
        if numNotUsed==0
            groupData.notUsed{1}=fnameFirstNucleiImg;
            groupData.minlength(1)=minlength;
        else
            groupData.notUsed{end+1}=fnameFirstNucleiImg;
            groupData.minlength(end+1)=minlength;
        end
        % save('xGroupedData.mat','groupData','-v7.3');
        return;
    end
    
    %**********************************************************************
    % Perform the local neigh measurement
    %**********************************************************************
    %myMean=@(x,y) mean([repmat(x,length(y),1) y],2);
    %retX=@(x,y) x;
    %retY=@(x,y) y;
    
    % Don't store all the values, that would be crazy:
    display('Measure local neighborhood:...')
    for frame=1:dFramesCorr:numFrames
%         text=['Frame: ',num2str(frame,'%03.f'),' of ',num2str(numFrames)];
%         progressText(frame/numFrames,text);
        
        % find the tracks which start at this frame
        xPos = tracksMatxCord(:,frame:numFrames); % Must be the same as tracksMatxCord because of cellDistClean (keep minTrackLength same)!!!
        yPos = tracksMatyCord(:,frame:numFrames); % Must be the same as tracksMatyCord because of cellDistClean (keep minTrackLength same)!!!
        c2eD =  cell2edgeDist(:,frame:numFrames);
        
        % Take only cells that are there at the start frame:
        cv1=~isnan(xPos(:,1));
        xPos=xPos(cv1,:);
        yPos=yPos(cv1,:);
        c2eD=c2eD(cv1,:);
        
        numCells=size(xPos,1);
        
        % This is to find the pairs of cell indices which correspond to
        % each row of pdist:
        % p is simply the ind in the vector:
        p=(1:(numCells-1)*numCells/2)';
        % vi is the index of "Spaltenumbruch" in the matrix:
        vi=1:(numCells-1);
        vi=vi*numCells-vi.*(vi+1)/2;        
%         try
%             % find the corresponding indices with mex function:
%             c2=comp2Vecs(p,vi');
%         catch
            % The corresponding Matlab code which is five times slower than
            % the mex function:
            c2=zeros((numCells-1)*numCells/2,1);
            for i=1:length(p)
              c2(i)=sum(vi<p(i),2)+1;
            end
%         end
        c1=p+c2.*(c2+1)/2-numCells*(c2-1);
        % The construction of c1 and c2 has been checked thoroughly!
        
        c2cDvec    =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
        c2edCOM    =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
        dirDisMean =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
        trvDisMean =zeros(numCells*(numCells-1)/2,numFrames-frame+1);
        % travDisMean is the mean of the distances traveled by the pair of
        % cells.            
            
        if numCells>1 % in case that there is at least 1 pair:
            for subframe=frame:numFrames
                % Calculate the cell to cell distance, the distance of the
                % center of mass to the wound edge and the direct distance and the traveled distance:
                c2cDvec(:,subframe-frame+1)= pdist([xPos(:,subframe-frame+1) yPos(:,subframe-frame+1)]);
                % c2edCOM(:,subframe-frame+1)= pdist( c2eD(:,subframe-frame+1),myMean);

                dirDist = sqrt((xPos(:,subframe-frame+1)-xPos(:,1)).^2+(yPos(:,subframe-frame+1)-yPos(:,1)).^2);

                trvDist = sum(sqrt((xPos(:,2:(subframe-frame+1))-xPos(:,1:(subframe-frame))).^2+(yPos(:,2:(subframe-frame+1))-yPos(:,1:(subframe-frame))).^2),2);

                %idVec=(1:numCells)';
                %c2 = pdist(idVec,retX)';
                %c1 = pdist(idVec,retY)';

                % [c1 c2]
                dirDisC1=dirDist(c1);
                dirDisC2=dirDist(c2);
                dirDisMean(:,subframe-frame+1)=mean([dirDisC1 dirDisC2],2);

                trvDisC1=trvDist(c1);
                trvDisC2=trvDist(c2);
                trvDisMean(:,subframe-frame+1)=mean([trvDisC1 trvDisC2],2);


                c2edC1=c2eD(c1,subframe-frame+1);
                c2edC2=c2eD(c2,subframe-frame+1);            
                c2edCOM(:,subframe-frame+1)=mean([c2edC1 c2edC2],2);
            end       
        else % there is not a single cell pair, might happen in single cell data:
            c2cDvec=NaN;
            c2edCOM=NaN;
            dirDisMean=NaN;
            trvDisMean=NaN;
        end
        % Take only the pairs which are there at the beginning and which
        % are close together:
        maxPairDist=50;
        cv1=~isnan(c2cDvec(:,1));
        cv2=c2cDvec(:,1)<maxPairDist;
        cvtot=cv1&cv2;
        
        c2cDvecOI    = c2cDvec(cvtot,:);
        c2edCOMOI    = c2edCOM(cvtot,:);
        dirDisMeanOI = dirDisMean(cvtot,:);
        trvDisMeanOI = trvDisMean(cvtot,:);
        
        % bin those according to the distance to the wound edge at the very
        % beginning:
        if maxDist==0
           bin=ones(size(c2edCOMOI(:,1)));
        else
           [~,bin]=histc(c2edCOMOI(:,1),0:bandWidth:maxDist+binPix);
        end
        
        if ~isempty(bin)
            for binID=1:max(bin)
                % find all cells that are in bin:
                ptsId=(bin==binID);
                % store the d(t) for each pair:
                neigh(binID,frame).c2cDis  =    c2cDvecOI(ptsId,:);
                neigh(binID,frame).c2cDD   =    c2cDvecOI(ptsId,:)-repmat(c2cDvecOI(ptsId,1),1,size(c2cDvecOI,2)); % change in relative distance over time        
                neigh(binID,frame).dirDis  = dirDisMeanOI(ptsId,:);  % This is the mean value of the cell pair
                neigh(binID,frame).trvDis  = trvDisMeanOI(ptsId,:);  % This is the mean value of the cell pair
                
                % Also store the average distance to the wound edge:
                neigh(binID,frame).c2edMean = mean(c2edCOMOI(ptsId));
                neigh(binID,frame).c2edSTD  =  std(c2edCOMOI(ptsId));
            end
        else  % this happens for single cell migration if there is not a single pair
            neigh(1,frame).c2cDis  = [];
            neigh(1,frame).c2cDD   = [];
            neigh(1,frame).dirDis  = [];  % This is the mean value of the cell pair
            neigh(1,frame).trvDis  = [];  % This is the mean value of the cell pair
            neigh(1,frame).c2edMean = [];
            neigh(1,frame).c2edSTD  = [];
        end
        clear c2cDvec c2edCOM dirDisMean c2edMean;
    end
    display('done');
    % plotNeighResults(neigh,maxDist,dFramesCorr)
    % close all
    
    %**********************************************************************
    % Perform the correlation analysis (takes a long time!)
    %**********************************************************************
    % prefilter the velocity matrix:
    display('Measure correlations:...')
    [posMatxCorr,posMatyCorr,vMatxCorr,vMatyCorr,~]=conv2CordMatVelMat(tracksMatxCord,tracksMatyCord,minTrackLength,timeWindowCorr); %,toDoList);???
    for frame=1:dFramesCorr:(numFrames-timeWindowCorr)
        
%         text=['Frame: ',num2str(frame,'%03.f'),' of ',num2str(numFrames-timeWindowCorr)];
%         progressText(frame/(numFrames-timeWindowCorr),text);

        %Remove the NaNs:
        checkVec=~isnan(vMatxCorr(:,frame));
        
        xPos    = posMatxCorr(checkVec,frame); % Must be the same as tracksMatxCord because of cellDistClean (keep minTrackLength same)!!!
        yPos    = posMatyCorr(checkVec,frame); % Must be the same as tracksMatyCord because of cellDistClean (keep minTrackLength same)!!!  
        xVel    = vMatxCorr(checkVec,frame);
        yVel    = vMatyCorr(checkVec,frame);                    
        cellDistClean = cell2edgeDist(checkVec,frame);
        
        % Bin the data according to the cell to edge distance:
        if maxDist==0
            bin=ones(size(cellDistClean(:,1)));
        else
            [~,bin]=histc(cellDistClean,0:bandWidth:maxDist+binPix);
        end
        
        for binID=1:max(bin)
            % find all cells that are in bin:
            ptsId=(bin==binID);
            % calculate the velocity correlation function:
            [currCorr]=velCorrFunc(xPos,yPos,xVel,yVel,dR,ptsId);
            
            % Also store the average distance to the wound edge:
            currCorr.c2edMean = mean(cellDistClean(ptsId));
            currCorr.c2edSTD  = std(cellDistClean(ptsId));
            
            % store the correlation functions for this frame and distance
            % to the wound edge. Also store R?!:           
            corr(binID,frame)=currCorr;
        end
    end
    % plotCorrResults(corr,1000,24);
    % close all;
    display('done');
    
    %**********************************************************************
    % Sort in the data set
    %**********************************************************************
    if numDataSets==0
        % Sort in the first value:
        classID = numClasses+1;
        setId   = numDataSets+1;

        % these values are used to check where to sort in the new data set:
        kPaClass(classID).yModu_kPa      = yModu_kPa;
        kPaClass(classID).cc             = cc;
        kPaClass(classID).cond           = cond;

        kPaClass(classID).minTrackLength = minTrackLength;
        kPaClass(classID).timeWindow     = timeWindow;
        kPaClass(classID).binPix         = densityMeasurement.binPix;
        kPaClass(classID).dt             = dt;
        kPaClass(classID).pixSize_um     = pixSize_um;

        numSetsInCl=1;
        kPaClass(classID).numSetsInCl               = numSetsInCl;
        kPaClass(classID).maxFrame                  = [toDoList(end) toDoList(end)];
        kPaClass(classID).setIdInList(numSetsInCl)  = setId;
        kPaClass(classID).signFac                   = repmat(signFac,size(cell2edgeDist,1),1);
        kPaClass(classID).tracksMatxCord            = tracksMatxCord;
        kPaClass(classID).tracksMatyCord            = tracksMatyCord;
        kPaClass(classID).cell2edgeDist             = cell2edgeDist;
        kPaClass(classID).velMatxCord               = signFac*velMatxCord;
        kPaClass(classID).velMatyCord               = velMatyCord;
        kPaClass(classID).velMatMag                 = velMatMag;
        kPaClass(classID).coveredDist               = coveredDist;

        % Quick fix for single cell data. Should be removed later on:
        [~,numFramesNew  ]=size(cell2edgeDist);
        if sum(isnan(roughness))==length(roughness) && numFramesNew<length(roughness)
            roughness(numFramesNew+1:end)=[];
        end
        kPaClass(classID).roughness                 = roughness;
        
        kPaClass(classID).binData(numSetsInCl).cells       = densityMeasurement.cells;
        kPaClass(classID).binData(numSetsInCl).area        = densityMeasurement.area;
        kPaClass(classID).binData(numSetsInCl).normDensity = densityMeasurement.density/nanmean(densityMeasurement.density(:,1));
        kPaClass(classID).binData(numSetsInCl).zeroDensity = nanmean(densityMeasurement.density(:,1));
        
        kPaClass(classID).freqSpec(numSetsInCl).ssA        = freqSpec.ssA;
        kPaClass(classID).freqSpec(numSetsInCl).f          = freqSpec.f;
        
        kPaClass(classID).corr         = corr;        
        kPaClass(classID).neigh = neigh;
            

        % these are global properties of the grouped data.
        groupData.dataSetList{setId} = fnameFirstNucleiImg;
        groupData.glbMaxFrame        = kPaClass(classID).maxFrame ;

        % put in the first kPaClass:
        groupData.kPaClass(classID)  = kPaClass(classID);
    else
        % We have already checked that the file isn't in there
        setId   = numDataSets+1;
        groupData.dataSetList{setId} = fnameFirstNucleiImg;
        groupData.glbMaxFrame        = [min(groupData.glbMaxFrame(1),toDoList(end))    max(groupData.glbMaxFrame(2),toDoList(end))];
        
        % Find the appropriate kPaClass:
        found=0;
        currClass=1;
        while found==0 && currClass<=numClasses
            if groupData.kPaClass(currClass).yModu_kPa==yModu_kPa && strcmp(groupData.kPaClass(currClass).cc,cc) && strcmp(groupData.kPaClass(currClass).cond,cond) && groupData.kPaClass(currClass).minTrackLength==minTrackLength && groupData.kPaClass(currClass).timeWindow==timeWindow && groupData.kPaClass(currClass).dt==dt && groupData.kPaClass(currClass).pixSize_um == pixSize_um
                found=1;
                classID = currClass;
            else
                currClass=currClass+1;
            end
        end
        
        if found==1
            % sort in the values:
            % first find the frame number of xGroupData and the new data set:
            [~         ,numFramesGroup]=size(groupData.kPaClass(classID).cell2edgeDist);
            [numColsNew,numFramesNew  ]=size(cell2edgeDist);
            maxFrames=max(numFramesGroup,numFramesNew);
            
            numSetsInCl   = groupData.kPaClass(classID).numSetsInCl + 1;
            groupData.kPaClass(classID).numSetsInCl             = numSetsInCl;
            groupData.kPaClass(classID).maxFrame                = [min(groupData.kPaClass(classID).maxFrame(1),toDoList(end))    max(groupData.kPaClass(classID).maxFrame(2),toDoList(end))];
            groupData.kPaClass(classID).setIdInList(numSetsInCl)= setId;
            
            groupData.kPaClass(classID).signFac                 = vertcat(groupData.kPaClass(classID).signFac                  , repmat(signFac,numColsNew,1));
            groupData.kPaClass(classID).tracksMatxCord          = vertcat(padarray(groupData.kPaClass(classID).tracksMatxCord  , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(tracksMatxCord     , [0,maxFrames-numFramesNew], NaN,'post'));
            groupData.kPaClass(classID).tracksMatyCord          = vertcat(padarray(groupData.kPaClass(classID).tracksMatyCord  , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(tracksMatyCord     , [0,maxFrames-numFramesNew], NaN,'post'));            
            groupData.kPaClass(classID).cell2edgeDist           = vertcat(padarray(groupData.kPaClass(classID).cell2edgeDist   , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(cell2edgeDist      , [0,maxFrames-numFramesNew], NaN,'post'));
            groupData.kPaClass(classID).velMatxCord             = vertcat(padarray(groupData.kPaClass(classID).velMatxCord     , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(signFac*velMatxCord, [0,maxFrames-numFramesNew], NaN,'post'));
            groupData.kPaClass(classID).velMatyCord             = vertcat(padarray(groupData.kPaClass(classID).velMatyCord     , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(        velMatyCord, [0,maxFrames-numFramesNew], NaN,'post'));
            groupData.kPaClass(classID).velMatMag               = vertcat(padarray(groupData.kPaClass(classID).velMatMag       , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(          velMatMag, [0,maxFrames-numFramesNew], NaN,'post'));
            groupData.kPaClass(classID).coveredDist             = vertcat(padarray(groupData.kPaClass(classID).coveredDist     , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(coveredDist        , [0,maxFrames-numFramesNew], NaN,'post'));
            
            % Quick fix for single cell data. Should be removed later on:
            if sum(isnan(roughness))==length(roughness) && numFramesNew<length(roughness)
                roughness(numFramesNew+1:end)=[];
            end
            groupData.kPaClass(classID).roughness               = vertcat(padarray(groupData.kPaClass(classID).roughness       , [0,maxFrames-numFramesGroup], NaN,'post'),padarray(roughness          , [0,maxFrames-numFramesNew], NaN,'post'));
            
            groupData.kPaClass(classID).binData(numSetsInCl).cells        = densityMeasurement.cells;
            groupData.kPaClass(classID).binData(numSetsInCl).area         = densityMeasurement.area;
            groupData.kPaClass(classID).binData(numSetsInCl).normDensity  = densityMeasurement.density/nanmean(densityMeasurement.density(:,1));
            groupData.kPaClass(classID).binData(numSetsInCl).zeroDensity  = nanmean(densityMeasurement.density(:,1)); 
            % padd all with NaNs:
            [grpRows,grpCols]=size(groupData.kPaClass(classID).binData(1).normDensity);
            [newRows,newCols]=size(groupData.kPaClass(classID).binData(numSetsInCl).normDensity);
            for msrm=1:numSetsInCl
                [curRows,curCols]=size(groupData.kPaClass(classID).binData(msrm).normDensity);
                groupData.kPaClass(classID).binData(msrm).normDensity=padarray(groupData.kPaClass(classID).binData(msrm).normDensity, [max(grpRows,newRows)-curRows,max(grpCols,newCols)-curCols], NaN,'post');
                groupData.kPaClass(classID).binData(msrm).cells      =padarray(groupData.kPaClass(classID).binData(msrm).cells      , [max(grpRows,newRows)-curRows,max(grpCols,newCols)-curCols], 0  ,'post');
                groupData.kPaClass(classID).binData(msrm).area       =padarray(groupData.kPaClass(classID).binData(msrm).area       , [max(grpRows,newRows)-curRows,max(grpCols,newCols)-curCols], 0  ,'post');
            end
            
            % Insert the roughness measures:
            groupData.kPaClass(classID).freqSpec(numSetsInCl).ssA = freqSpec.ssA;
            groupData.kPaClass(classID).freqSpec(numSetsInCl).f   = freqSpec.f;            
            [grpFreqs,grpFrames]=size(groupData.kPaClass(classID).freqSpec(1).ssA);
            [newFreqs,newFrames]=size(groupData.kPaClass(classID).freqSpec(numSetsInCl).ssA);
            for msrm=1:numSetsInCl
                [curFreqs,curFrames]=size(groupData.kPaClass(classID).freqSpec(msrm).ssA);
                groupData.kPaClass(classID).freqSpec(msrm).ssA=padarray(groupData.kPaClass(classID).freqSpec(msrm).ssA, [max(grpFreqs,newFreqs)-curFreqs,max(grpFrames,newFrames)-curFrames], NaN,'post');
                groupData.kPaClass(classID).freqSpec(msrm).f  =padarray(groupData.kPaClass(classID).freqSpec(msrm).f  , [max(grpFreqs,newFreqs)-curFreqs,0], NaN,'post');
            end
            
            groupData.kPaClass(classID).corr         = groupCorrResults( groupData.kPaClass(classID).corr        ,corr);
            groupData.kPaClass(classID).neigh = groupNeighResults(groupData.kPaClass(classID).neigh,neigh);

        else
            classID = numClasses+1;
            
            % these values are used to check where to sort in the new data set:
            kPaClass(classID).yModu_kPa      = yModu_kPa;
            kPaClass(classID).cc             = cc;
            kPaClass(classID).cond           = cond;
            kPaClass(classID).minTrackLength = minTrackLength;
            kPaClass(classID).timeWindow     = timeWindow;
            kPaClass(classID).binPix         = densityMeasurement.binPix;
            kPaClass(classID).dt             = dt;
            kPaClass(classID).pixSize_um     = pixSize_um;
            
            numSetsInCl=1;
            kPaClass(classID).numSetsInCl               = numSetsInCl;
            kPaClass(classID).maxFrame                  = [toDoList(end) toDoList(end)];
            kPaClass(classID).setIdInList(numSetsInCl)  = setId;
            kPaClass(classID).signFac                   = repmat(signFac,size(cell2edgeDist,1),1);
            kPaClass(classID).tracksMatxCord            = tracksMatxCord;
            kPaClass(classID).tracksMatyCord            = tracksMatyCord;
            kPaClass(classID).cell2edgeDist             = cell2edgeDist;
            kPaClass(classID).velMatxCord               = signFac*velMatxCord;
            kPaClass(classID).velMatyCord               = velMatyCord;
            kPaClass(classID).velMatMag                 = velMatMag;
            kPaClass(classID).coveredDist               = coveredDist;
            
            % Quick fix for single cell data. Should be removed later on:
            [~,numFramesNew  ]=size(cell2edgeDist);
            if sum(isnan(roughness))==length(roughness) && numFramesNew<length(roughness)
                roughness(numFramesNew+1:end)=[];
            end
            kPaClass(classID).roughness                 = roughness;
        
            kPaClass(classID).roughness                 = roughness;
            
            kPaClass(classID).binData(numSetsInCl).cells       = densityMeasurement.cells;
            kPaClass(classID).binData(numSetsInCl).area        = densityMeasurement.area;
            kPaClass(classID).binData(numSetsInCl).normDensity = densityMeasurement.density/nanmean(densityMeasurement.density(:,1));
            kPaClass(classID).binData(numSetsInCl).zeroDensity = nanmean(densityMeasurement.density(:,1));
            
            kPaClass(classID).freqSpec(numSetsInCl).ssA        = freqSpec.ssA;
            kPaClass(classID).freqSpec(numSetsInCl).f          = freqSpec.f;
            
            kPaClass(classID).corr         = corr;
            kPaClass(classID).neigh = neigh;
            
            groupData.kPaClass(classID)  = kPaClass(classID);
        end
    end
    % save('xGroupedData.mat','groupData','-v7.3');
end


% tracksMatxCord(~checkVec,:)=[];
if justPlot==1
    dframes=24;
    marker=['ro','bs','m*','c+','gd','yo','ks'];
    % plot the results for each kPaClass:
    for classID=1:length(groupData.kPaClass)
        figure()
        label=1;
        % here one could take either the global or class specific values
        % for the smallest/largest value for the last frame.
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        
        % doesn't work because of the sign-factor!!!
        % instead you could also put in your own values:        
        timeWindowPlot     = groupData.kPaClass(classID).timeWindow;
        minTrackLengthPlot = groupData.kPaClass(classID).minTrackLength;

        timeWindowPlot=6;
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        dt = groupData.kPaClass(classID).dt;
        fac_um_min=pixSize_um/(timeWindowPlot*dt*60);
        
        [xPlot,yPlot,vxPlot,vyPlot,vPlot,checkVecPlot]=conv2CordMatVelMat(groupData.kPaClass(classID).tracksMatxCord,groupData.kPaClass(classID).tracksMatyCord,minTrackLengthPlot,timeWindowPlot); %,toDoList);???
        c2edPlot=groupData.kPaClass(classID).cell2edgeDist(checkVecPlot,:);
        signFac =groupData.kPaClass(classID).signFac(checkVecPlot);
        vxPlot  =repmat(signFac,[1,size(vxPlot,2)]).*vxPlot;
        

        
        for frame=1:dframes:(numFrames-timeWindowPlot)
            %plot(tracksMatxCord(:,frame),
            plot(pixSize_um*c2edPlot(:,frame),fac_um_min*vxPlot(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
            % count the number of points in the plot:
            numPts(classID,frame)=sum(~isnan(c2edPlot(:,frame)) & ~isnan(vxPlot(:,frame)));
            hold on
            M{label}=[num2str((frame-1)*dt),'h',' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),'; n=',num2str(numPts(classID,frame))];
            label=label+1;
        end
        xlabel('Distance to the wound edge in [um]')
        ylabel('x-component of the velocity [um/min]')
        title(['Plot for stiffness: ',num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1}])
        legend(M);
        % saveas(gcf,['fig_x_vel_over_dist_to_edge_raw',num2str(groupData.kPaClass(classID).yModu_kPa),'.tiff'],'tiffn');
        clear M
        hold off
    end


    for classID=1:length(groupData.kPaClass)
        binPix =groupData.kPaClass(classID).binPix;
        maxDist=max(groupData.kPaClass(classID).cell2edgeDist(:,end));
        figure()
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        label=1;
        h=[];       
        
        % instead you could also put in your own values:        
        timeWindowPlot     = groupData.kPaClass(classID).timeWindow;
        minTrackLengthPlot = groupData.kPaClass(classID).minTrackLength;
        
        timeWindowPlot=6;
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        fac_um_min=pixSize_um/(timeWindowPlot*groupData.kPaClass(classID).dt*60);
        dt = groupData.kPaClass(classID).dt;
        
        [xPlot,yPlot,vxPlot,vyPlot,vPlot,checkVecPlot]=conv2CordMatVelMat(groupData.kPaClass(classID).tracksMatxCord,groupData.kPaClass(classID).tracksMatyCord,minTrackLengthPlot,timeWindowPlot); %,toDoList);???
        c2edPlot=groupData.kPaClass(classID).cell2edgeDist(checkVecPlot,:);
        signFac =groupData.kPaClass(classID).signFac(checkVecPlot);
        vxPlot  = repmat(signFac,[1,size(vxPlot,2)]).*vxPlot;

        for frame=1:dframes:(numFrames-timeWindowPlot)
            %Remove the NaNs:
            checkVec=~isnan(vxPlot(:,frame));
            velVecxClean =vxPlot(checkVec,frame);
            cellDistClean=c2edPlot(checkVec,frame);

            % bin the data:
            [n,bin]=histc(cellDistClean,0:binPix:maxDist+binPix);

            for binID=1:max(bin)
                velVecxMean(binID)=mean(velVecxClean(bin==binID));
                velVecxSEM(binID) =std(velVecxClean(bin==binID))/sqrt(n(binID));

                cellDistMean(binID)=mean(cellDistClean(bin==binID));
                cellDistSEM(binID) =std(cellDistClean(bin==binID))/sqrt(n(binID));
            end

            currh=errorbarxy(pixSize_um*cellDistMean,fac_um_min*velVecxMean,pixSize_um*cellDistSEM,fac_um_min*velVecxSEM,[],[],marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)),marker(2*(mod(frame,7)+1)-1));
            h=horzcat(h,currh(1));
            hold on
            M{label}=[num2str((frame-1)*dt),'h',' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),'; n=',num2str(numPts(classID,frame))];
            label=label+1;
            
            clear velVecxMean velVecxSEM cellDistMean cellDistSEM
        end
        xlabel('Distance to the wound edge in [um]')
        ylabel('x-component of the velocity [um/min')
        title(['Plot for stiffness: ',num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1}])
        legend(h,M);
        % saveas(gcf,['fig_x_vel_over_dist_to_edge_cum',num2str(groupData.kPaClass(classID).yModu_kPa),'.tiff'],'tiffn');
        clear M
        hold off
    end
    
    
    for classID=1:length(groupData.kPaClass)
        binPix =groupData.kPaClass(classID).binPix;
        maxDist=max(groupData.kPaClass(classID).cell2edgeDist(:,end));
        figure()
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        label=1;
        h=[];
        
        % instead you could also put in your own values:        
        timeWindowPlot     = groupData.kPaClass(classID).timeWindow;
        minTrackLengthPlot = groupData.kPaClass(classID).minTrackLength;
        
        timeWindowPlot=6;
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        fac_um_min=pixSize_um/(timeWindowPlot*groupData.kPaClass(classID).dt*60);
        dt = groupData.kPaClass(classID).dt;

        [xPlot,yPlot,vxPlot,vyPlot,vPlot,checkVecPlot]=conv2CordMatVelMat(groupData.kPaClass(classID).tracksMatxCord,groupData.kPaClass(classID).tracksMatyCord,minTrackLengthPlot,timeWindowPlot); %,toDoList);???
        c2edPlot=groupData.kPaClass(classID).cell2edgeDist(checkVecPlot,:);
        

        for frame=1:dframes:(numFrames-timeWindowPlot)
            %Remove the NaNs:
            checkVec      = ~isnan(vPlot(:,frame));
            speedClean    =        vPlot(checkVec,frame);
            cellDistClean =     c2edPlot(checkVec,frame);
            
            % bin the data:
            [n,bin]=histc(cellDistClean,0:binPix:maxDist+binPix);
            
            for binID=1:max(bin)
                speedMean(binID)=mean(speedClean(bin==binID));
                speedSEM(binID) = std(speedClean(bin==binID))/sqrt(n(binID));
                
                cellDistMean(binID)=mean(cellDistClean(bin==binID));
                cellDistSEM(binID) = std(cellDistClean(bin==binID))/sqrt(n(binID));
            end
            
            currh=errorbarxy(pixSize_um*cellDistMean,fac_um_min*speedMean,pixSize_um*cellDistSEM,fac_um_min*speedSEM,[],[],marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)),marker(2*(mod(frame,7)+1)-1));
            h=horzcat(h,currh(1));
            hold on
            M{label}=[num2str((frame-1)*dt),'h',' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),'; n=',num2str(numPts(classID,frame))];
            label=label+1;
            
            clear speedMean speedSEM cellDistMean cellDistSEM
        end
        xlabel('Distance to the wound edge in [um]')
        ylabel('Speed [um/min]')
        title(['Plot for stiffness: ',num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1}])
        legend(h,M);
        % saveas(gcf,['fig_x_vel_over_dist_to_edge_cum',num2str(groupData.kPaClass(classID).yModu_kPa),'.tiff'],'tiffn');
        clear M
        hold off
    end
    
    figure()
    label=1;
    h=[];
    for classID=1:length(groupData.kPaClass)
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        dt = groupData.kPaClass(classID).dt;
        
        for frame=1:numFrames
            %Remove the NaNs:
            checkVec=~isnan(groupData.kPaClass(classID).coveredDist(:,frame));
            numVal=length(checkVec);
            
            if numVal~=groupData.kPaClass(classID).numSetsInCl
                display(['Inconsistency, not all results are plotted in figure for covered distance at frame: ', um2str(frame)])
            end
            
            covDistMean(frame) =mean(groupData.kPaClass(classID).coveredDist(checkVec,frame));
            covDistSEM(frame)  =std(groupData.kPaClass(classID).coveredDist(checkVec,frame))/sqrt(numVal);
        end
        currh=errorbar(dt*((1:numFrames)-1),pixSize_um*covDistMean,pixSize_um*covDistSEM,marker(2*(mod(classID,7)+1)-1:2*(mod(classID,7)+1)));
        h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
        hold on
        M{label}=[num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1},' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl)];
        label=label+1;
        
        clear covDistMean covDistSEM
    end
    xlabel('Time [h]')
    ylabel('Covered distance [um]')
    legend(h,M);
    %saveas(gcf,['fig_trav_dist','.tiff'],'tiffn');
    clear M
    hold off
    
    for classID=1:length(groupData.kPaClass)
        binPix =groupData.kPaClass(classID).binPix;
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        dt = groupData.kPaClass(classID).dt;
        figure()
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        label=1;
        h=[];
        % calculate the average densities for all frames and all bins:
        % Here we keep only bins that are present in all datasets (others
        % will be NaN)!
        densityAll  = cat(3,groupData.kPaClass(classID).binData.normDensity);
        densityMean = mean(densityAll,3);
        densitySEM  = std(densityAll,[],3)/sqrt(groupData.kPaClass(classID).numSetsInCl);
        
        % also sum up the number of cells
        cellsAll=cat(3,groupData.kPaClass(classID).binData.cells);
        cellsSum=sum(cellsAll,3);
        
        for frame=1:dframes:numFrames
            [numBins,~]=size(densityMean);
            % disregard the last two bins:
            dBin=0;
            currh=errorbar(pixSize_um*binPix*(1:numBins-dBin),densityMean(1:numBins-dBin,frame),densitySEM(1:numBins-dBin,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
            h=horzcat(h,currh(1));
            hold on
            numCells(classID,frame)=sum(cellsSum(1:numBins-dBin,frame));
            M{label}=[num2str((frame-1)*dt),'h',' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),' ; n=',num2str(numCells(classID,frame))];
            label=label+1;
        end
        xlabel('Distance to the wound edge in [um]')
        ylabel('Normalized cell density [1]')
        title(['Plot for stiffness: ',num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1}])
        legend(h,M);
        % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn');
        clear M
        hold off
    end
    % close all
    
    %**********************************************************************
    % Perform measurement of local neigh:
    %**********************************************************************
    figure('Name','Local neighborhood');
    clear M
    label=1;
    h=[];
    %figure('Name',['Velocity correlation, for ',num2str(groupData.kPaClass(classID).yModu_kPa),' kPa-Class = ',groupData.kPaClass(classID).cc{1}],'NumberTitle','off')
    numFrames=groupData.glbMaxFrame(1);
    numBins=4;
    for binR=1:numBins
        numPlotFrames=length(1:dframes:numFrames);
        k=1;
        for frame=1:dframes:numFrames
            subplot(numPlotFrames,numBins,(k-1)*numBins+binR)
            
            for classID=1:length(groupData.kPaClass)
                pixSize_um=groupData.kPaClass(classID).pixSize_um;
                
                % Calculate the mean separation speed:
                neigh=groupData.kPaClass(classID).neigh(binR,frame);
                
                c2edClassMean    = nanmean(neigh.c2edMean);
                c2edSTDClassMean = nanmean(neigh.c2edSTD);
                
                % c2cDis c2cDD dirDis trvDis c2edMean c2edSTD
                               
                % count paris in each bin:
                Ntot=sum(~isnan(neigh.c2cDD),1);
                c2cDDClassMean   =              nanmean(neigh.c2cDD,1);                
                c2cDDClassSEM95  = facSEMtoSEM95*nanstd(neigh.c2cDD,[],1)./sqrt(Ntot);
                
                % Change in cell 2 cell distance, normalized by the average
                % distance traveled by the cell pair:
                normDDClassMean  =              nanmean(neigh.c2cDD./neigh.trvDis,1);
                normDDClassSEM95 = facSEMtoSEM95*nanstd(neigh.c2cDD./neigh.trvDis,[],1)./sqrt(Ntot);
                
                % c2cDDClassMean
                frameCount=length(c2cDDClassMean);
                % plot the results:
                %currh=errorbar(dt*((1:frameCount)-1),pixSize_um*c2cDDClassMean ,pixSize_um*c2cDDClassSEM95 ,['.',marker(2*(mod(classID,7)+1)-1)]);
                currh=errorbar(dt*((1:frameCount)-1),normDDClassMean,normDDClassSEM95,marker(2*(mod(classID,7)+1)-1:2*(mod(classID,7)+1)));
                h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
                hold on
                M{label}=[num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1},' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),' ; <n>=',num2str(mean(Ntot),'%.1f')];
                label=label+1;
                hold on;
            end
            title(['Cell to edge distance: ',num2str(pixSize_um*c2edClassMean,'%.1f'),'+-',num2str(pixSize_um*c2edSTDClassMean,'%.1f'),'[um]'])
            ylabel(['Change of cell to cell dist, dD(r,t=',num2str(dt*(frame-1)),'h)']);
            xlabel('dt [h]');
            ylim([ 0 0.4]);
            xlim([ 0 20]);
            if frame ==1 && binR==numBins
                legend(h,M);
            end
            clear M
            label=1;
            h=[];
            hold off           
            k=k+1;
        end
    end
    % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn')
    
    
    % return;
    
    
    %**********************************************************************
    % Perform persistence measurement:
    %**********************************************************************
    figure('Name','Persistence');
    clear M
    label=1;
    h=[];
    %figure('Name',['Velocity correlation, for
    %',num2str(groupData.kPaClass(classID).yModu_kPa),' kPa-Class = ',groupData.kPaClass(classID).cc{1}],'NumberTitle','off')
    numBins=4;
    trgDFrames=24; % this = 4h
    numFrames=min([1+trgDFrames*4,groupData.glbMaxFrame(2)]);
    for binc2e=1:numBins
        numPlotFrames=length(1:dframes:(numFrames-trgDFrames));
        k=1;
        for frame=1:dframes:(numFrames-trgDFrames)
            subplot(numPlotFrames,numBins,(k-1)*numBins+binc2e)
            
            for classID=1:length(groupData.kPaClass)
                pixSize_um=groupData.kPaClass(classID).pixSize_um;
                dt = groupData.kPaClass(classID).dt;
        
                % cut out the relevant part of the coordinate and cell to edge
                % distance field:
                
                % adapt the number of dT steps:
                maxDFrames=min(trgDFrames,groupData.kPaClass(classID).maxFrame(2)-frame);
                
                xPos=groupData.kPaClass(classID).tracksMatxCord(:,frame:frame+maxDFrames);
                yPos=groupData.kPaClass(classID).tracksMatyCord(:,frame:frame+maxDFrames);
                c2ed=groupData.kPaClass(classID).cell2edgeDist(:,frame);
                
                % Bin the data according to the cell to edge distance:
                [~,bin]=histc(c2ed,0:bandWidth:maxDist+binPix);
                
                % find all cells that are in bin:
                ptsId=(bin==binc2e);
                xPos=xPos(ptsId,:);
                yPos=yPos(ptsId,:);
                
                % mean distance to the wound edge:
                c2edClassMean=mean(c2ed(ptsId));
                c2edSTDClassMean=std(c2ed(ptsId));
                
                % calculate the persistence function:
                %[pMean,pSEM95,Ntot]=calcPersistence(xPos,yPos);
                [pMean,pSEM95,Ntot,pxMean,pxSEM95,pyMean,pySEM95]=calcPersistence(xPos,yPos);
                
                % The persistence
                %currh=errorbar(dt*((1:maxDFrames)-1),pMean,pSEM95,['.',marker(2*(mod(classID,7)+1)-1)]);
                %currh=errorbar(dt*((1:maxDFrames)-1),pxMean,pxSEM95,['.',marker(2*(mod(classID,7)+1)-1)]);
                currh=errorbar(dt*((1:maxDFrames)-1),pyMean,pySEM95,['.',marker(2*(mod(classID,7)+1)-1)]);
                
                h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
                hold on
                M{label}=[num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1},' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),' ; <n>=',num2str(mean(Ntot),'%.1f')];
                label=label+1;
                hold on;
            end
            title(['Cell to edge distance: ',num2str(pixSize_um*c2edClassMean,'%.1f'),'+-',num2str(pixSize_um*c2edSTDClassMean,'%.1f'),'[um]'])
            ylabel(['Persistence p(r,t=',num2str(dt*(frame-1)),'h)']);
            xlabel('dt [h]');
            ylim([0 1]);
            xlim([ 0 dt*(maxDFrames-1)]);
            if frame ==1 && binc2e==numBins
                legend(h,M);
            end
            clear M
            label=1;
            h=[];
            hold off
            
            k=k+1;
        end
    end
    % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn');
    
    
    
    %**********************************************************************
    % Perform measurement of frequencies:
    %**********************************************************************
    for classID=1:length(groupData.kPaClass)
        binPix =groupData.kPaClass(classID).binPix;
        pixSize_um=groupData.kPaClass(classID).pixSize_um;
        dt = groupData.kPaClass(classID).dt;
        
        figure()
        numFrames=groupData.kPaClass(classID).maxFrame(1);
        label=1;
        h=[];
        % calculate the average densities for all frames and all bins:
        % Here we keep only bins that are present in all datasets (others
        % will be NaN)!
        
        %     kPaClass(classID).freqSpec(numSetsInCl).ssA        = freqSpec.ssA;
        %     kPaClass(classID).freqSpec(numSetsInCl).f          = f;
        
        % group the spectra:
        specAll    = cat(3,groupData.kPaClass(classID).freqSpec(:).ssA);
        specMean   = mean(specAll,3);
        specSEM95  = facSEMtoSEM95*std(specAll,[],3)/sqrt(groupData.kPaClass(classID).numSetsInCl);
        
        % group the frequency
        fAll   = horzcat(groupData.kPaClass(classID).freqSpec.f);
        fMean  = nanmean(fAll,2);
        % This should be zero:
        fSEM95 = facSEMtoSEM95*std(fAll,[],2)/sqrt(groupData.kPaClass(classID).numSetsInCl);
        
        if sum(fSEM95)>eps
            error('Frequency grouping went wrong');
        end
        
        for frame=1:dframes:numFrames
            %currh=errorbar(fMean/pixSize_um,pixSize_um*specMean(:,frame),pixSize_um*specSEM95(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
            currh=plot(fMean/pixSize_um,pixSize_um*specMean(:,frame),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
            
            h=horzcat(h,currh(1));
            hold on
            M{label}=[num2str((frame-1)*dt),'h',' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl)];
            label=label+1;
        end
        xlabel('frequency [1/um]')
        ylabel('amplitude [um]')
        title(['Plot for stiffness: ',num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1}])
        legend(h,M);
        xlim([0,0.01])
        ylim([0,200])
        % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn');
        clear M
        hold off
    end
    
    
    % figure('Name','Velocity correlation');
    % for classID=1:length(groupData.kPaClass)
    %     %figure('Name',['Velocity correlation, for ',num2str(groupData.kPaClass(classID).yModu_kPa),' kPa-Class = ',groupData.kPaClass(classID).cc{1}],'NumberTitle','off')
    %     numFrames=groupData.kPaClass(classID).maxFrame(1);
    %     numBins=4;
    %     for binR=1:numBins
    %         numPlotFrames=length(1:dframes:(numFrames-timeWindowCorr));
    %         k=1;
    %         for frame=1:dframes:(numFrames-timeWindowCorr)
    %             subplot(numPlotFrames,numBins,(k-1)*numBins+binR)
    %
    %             % calculate the average correlation function for a certain r,t:
    %             % Thereby disregard all columns where one or more NaN's appear
    %             currCorr=groupData.kPaClass(classID).corr(binR,frame);
    %             [numRows,numCols]=size(currCorr.RMean);
    %
    %             RClassMean=sum(currCorr.RMean.*currCorr.Ntot)./sum(currCorr.Ntot);
    %             % RClassSTD =?
    %
    %             % First normalize then average:
    %             % average the correlations normalized by the autocorrelation:
    %             funcVelClassMean  = sum((currCorr.funcVel ./padarray(currCorr.funcVel(:,1) ,[0 numCols-1],'replicate','post')) .*currCorr.Ntot)./sum(currCorr.Ntot);
    %             funcVelxClassMean = sum((currCorr.funcVelx./padarray(currCorr.funcVelx(:,1),[0 numCols-1],'replicate','post')) .*currCorr.Ntot)./sum(currCorr.Ntot);
    %             funcVelyClassMean = sum((currCorr.funcVely./padarray(currCorr.funcVely(:,1),[0 numCols-1],'replicate','post')) .*currCorr.Ntot)./sum(currCorr.Ntot);
    %
    %             % average the correlations normalized by the autocorrelation:
    %             RMean4CosaClassMean = sum(currCorr.RMean4CosaMean.*currCorr.Ntot)./sum(currCorr.Ntot);
    %
    %             % average the correlations (normalized by the autocorrelation):
    %             %cosaMeanClassMean = sum((currCorr.cosaMean ./padarray(currCorr.cosaMean(:,1) ,[0 numCols-1],'replicate','post')) .*currCorr.Ntot)./sum(currCorr.Ntot);
    %             cosaMeanClassMean = sum(currCorr.cosaMean.*currCorr.Ntot)./sum(currCorr.Ntot);
    %
    %             % average Cell to edge distance for this bin:
    %             c2edClassMean    = mean(currCorr.c2edMean);
    %             c2edSTDClassMean = mean(currCorr.c2edSTD);
    %
    %             % First normalize then average:
    %             % plot(RClassMean,funcVelClassMean,'sb');
    %             hold on;
    %
    %             % Calculate the correlation from scratch:
    %             R=groupData.kPaClass(classID).corr(binR,frame).R;
    %             for c2cRBin=1:length(R)
    %                     NtotClass(c2cRBin)        = sum(R(c2cRBin).Ntot);
    %                     RClassMean(c2cRBin)       = nanmean(R(c2cRBin).pts2ptR);
    %                     RClassSEM95(c2cRBin)      = facSEMtoSEM95*nanstd(R(c2cRBin).pts2ptR)/sqrt(NtotClass(c2cRBin));
    %                     [I, Isem95]               = corrFunc(R(c2cRBin).xVelCtr,R(c2cRBin).xVelPts,R(c2cRBin).yVelCtr,R(c2cRBin).yVelPts);
    %                     funcVelClassMean(c2cRBin) = I;
    %                     funcVelClassSEM95(c2cRBin)= Isem95;
    %                     funcVelxClassMean(c2cRBin)= corrFunc(R(c2cRBin).xVelCtr,R(c2cRBin).xVelPts);
    %                     funcVelyClassMean(c2cRBin)= corrFunc(R(c2cRBin).yVelCtr,R(c2cRBin).yVelPts);
    %
    %                     % The cosa needs a special treatment, since here the caclulation could
    %                     % introduce additional NaNs (if vel=0). These positions have to be
    %                     % canceled out in the distance vector:
    %                     R4CosaClassMean(c2cRBin) = mean(R(c2cRBin).pts2ptR(~isnan(R(c2cRBin).cosda)));
    %                     R4CosaClassSTD(c2cRBin)  =  std(R(c2cRBin).pts2ptR(~isnan(R(c2cRBin).cosda)));
    %                     cosaClassMean(c2cRBin)   = nanmean(R(c2cRBin).cosda);
    %                     cosaClassSTD(c2cRBin)    =  nanstd(R(c2cRBin).cosda);
    %             end
    %             % First average then normalize:
    %             % Normalize the correlation by the autocorrelation (If
    %             % doNorm==1 in corrFunc nothing will be done in the folowing
    %             % lines since I(1)=1 anyways.):
    %             funcVelClassMeanNorm  = funcVelClassMean/funcVelClassMean(1);
    %             funcVelClassSEM95Norm = funcVelClassSEM95/funcVelClassMean(1);
    %
    %             errorbarxy(RClassMean,funcVelClassMeanNorm,RClassSEM95,funcVelClassSEM95Norm,[],[],'.r','r');
    %             %errorbar(RClassMean,funcVelClassMean,'sb');
    %             title(['Cell to edge distance: ',num2str(c2edClassMean,'%.1f'),'+-',num2str(c2edSTDClassMean,'%.1f'),'[Pix]'])
    %             ylabel(['Correlation I(r,t=',num2str(frame),')']);
    %             xlabel('R [pix]');
    %             ylim([-1 1]);
    %             xlim([ 0 1000]);
    %             hold on;
    %             k=k+1;
    %         end
    %     end
    %     % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn');
    % end
    
    
    
    figure('Name','Velocity correlation');
    clear M
    label=1;
    h=[];
    %figure('Name',['Velocity correlation, for ',num2str(groupData.kPaClass(classID).yModu_kPa),' kPa-Class = ',groupData.kPaClass(classID).cc{1}],'NumberTitle','off')
    numFrames=groupData.glbMaxFrame(1);
    numBins=4;
    for binR=1:numBins
        numPlotFrames=length(1:dframes:(numFrames-timeWindowCorr));
        k=1;
        for frame=1:dframes:(numFrames-timeWindowCorr)
            subplot(numPlotFrames,numBins,(k-1)*numBins+binR)
            
            for classID=1:length(groupData.kPaClass)
                pixSize_um=groupData.kPaClass(classID).pixSize_um;
                dt = groupData.kPaClass(classID).dt;
                
                c2edClassMean    = mean(groupData.kPaClass(classID).corr(binR,frame).c2edMean);
                c2edSTDClassMean = mean(groupData.kPaClass(classID).corr(binR,frame).c2edSTD);
                
                % Calculate the correlation from scratch:
                R=groupData.kPaClass(classID).corr(binR,frame).R;
                for c2cRBin=1:length(R)
                    NtotClass(c2cRBin)         = sum(R(c2cRBin).Ntot);
                    RClassMean(c2cRBin)        = nanmean(R(c2cRBin).pts2ptR);
                    RClassSEM95(c2cRBin)       = facSEMtoSEM95*nanstd(R(c2cRBin).pts2ptR)/sqrt(NtotClass(c2cRBin));
                    % The correlation of the full velocity:
                    [I, Isem95]                = corrFunc(R(c2cRBin).xVelCtr,R(c2cRBin).xVelPts,R(c2cRBin).yVelCtr,R(c2cRBin).yVelPts);
                    funcVelClassMean(c2cRBin)  = I;
                    funcVelClassSEM95(c2cRBin) = Isem95;
                    % The correlation of ux:
                    [I, Isem95]                = corrFunc(R(c2cRBin).xVelCtr,R(c2cRBin).xVelPts);
                    funcVelxClassMean(c2cRBin) = I;
                    funcVelxClassSEM95(c2cRBin)= Isem95;
                    % The correlation of uy:
                    [I, Isem95]                = corrFunc(R(c2cRBin).yVelCtr,R(c2cRBin).yVelPts);
                    funcVelyClassMean(c2cRBin) = I;
                    funcVelyClassSEM95(c2cRBin)= Isem95;
                    
                    
                    % The cosa needs a special treatment, since here the caclulation could
                    % introduce additional NaNs (if vel=0). These positions have to be
                    % canceled out in the distance vector:
                    R4CosaClassMean(c2cRBin)  = mean(R(c2cRBin).pts2ptR(~isnan(R(c2cRBin).cosda)));
                    R4CosaClassSEM95(c2cRBin) = facSEMtoSEM95*std(R(c2cRBin).pts2ptR(~isnan(R(c2cRBin).cosda)))/sqrt(sum(~isnan(R(c2cRBin).cosda)));
                    cosaClassMean(c2cRBin)    = nanmean(R(c2cRBin).cosda);
                    cosaClassSEM95(c2cRBin)   = facSEMtoSEM95*nanstd(R(c2cRBin).cosda)/sqrt(sum(~isnan(R(c2cRBin).cosda)));
                end
                % First average then normalize:
                % Normalize the correlation by the autocorrelation (If
                % doNorm==1 in corrFunc nothing will be done in the folowing
                % lines since I(1)=1 anyways.):
                funcVelClassMeanNorm  = funcVelClassMean /funcVelClassMean(1);
                funcVelClassSEM95Norm = funcVelClassSEM95/funcVelClassMean(1);
                
                funcVelxClassMeanNorm  = funcVelxClassMean /funcVelxClassMean(1);
                funcVelxClassSEM95Norm = funcVelxClassSEM95/funcVelxClassMean(1);
                
                funcVelyClassMeanNorm  = funcVelyClassMean /funcVelyClassMean(1);
                funcVelyClassSEM95Norm = funcVelyClassSEM95/funcVelyClassMean(1);
                
                % The correlation of the full velocity:
                currh=errorbarxy(pixSize_um*RClassMean,funcVelClassMeanNorm,pixSize_um*RClassSEM95,funcVelClassSEM95Norm,[],[],['.',marker(2*(mod(classID,7)+1)-1)],marker(2*(mod(classID,7)+1)-1));
                % currh=errorbar(pixSize_um*RClassMean,funcVelClassMeanNorm,funcVelClassSEM95Norm,['.',marker(2*(mod(classID,7)+1)-1)]);
                
                % The correlation by cosa:
                % currh=errorbar(pixSize_um*R4CosaClassMean,cosaClassMean,cosaClassSEM95,['+',marker(2*(mod(classID,7)+1)-1)]);
                
                % The correlation of ux:
                % currh=errorbarxy(pixSize_um*RClassMean,funcVelxClassMeanNorm,pixSize_um*RClassSEM95,funcVelxClassSEM95Norm,[],[],['.',marker(2*(mod(classID,7)+1)-1)],marker(2*(mod(classID,7)+1)-1));
                
                % The correlation of uy:
                % currh=errorbarxy(pixSize_um*RClassMean,funcVelyClassMeanNorm,pixSize_um*RClassSEM95,funcVelyClassSEM95Norm,[],[],['.',marker(2*(mod(classID,7)+1)-1)],marker(2*(mod(classID,7)+1)-1));
                
                h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
                hold on
                M{label}=[num2str(groupData.kPaClass(classID).yModu_kPa),'kPa',' = ',groupData.kPaClass(classID).cc{1},' of: ',groupData.kPaClass(classID).cond{1},' ; N=',num2str(groupData.kPaClass(classID).numSetsInCl),' ; <n>=',num2str(mean(NtotClass),'%.1f')];
                label=label+1;
                hold on;
            end
            title(['Cell to edge distance: ',num2str(pixSize_um*c2edClassMean,'%.1f'),'+-',num2str(pixSize_um*c2edSTDClassMean,'%.1f'),'[um]'])
            ylabel(['Correlation I(r,t=',num2str(dt*(frame-1)),'h)']);
            xlabel('R [um]');
            ylim([-0.1 1]);
            xlim([ 0 pixSize_um*1000]);
            if frame ==1 && binR==numBins
                legend(h,M);
            end
            clear M
            label=1;
            h=[];
            hold off
            
            clear NtotClass RClassMean RClassSEM95
            clear funcVelClassMean funcVelClassSEM95
            clear funcVelxClassMean funcVelxClassSEM95Norm
            clear funcVelyClassMean funcVelyClassSEM95Norm
            clear R4CosaClassMean R4CosaClassSEM95
            clear cosaClassMean cosaClassSEM95
            
            k=k+1;
        end
    end
    % saveas(gcf,['fig_x_vel_over_dist_to_edge','.tiff'],'tiffn')
end