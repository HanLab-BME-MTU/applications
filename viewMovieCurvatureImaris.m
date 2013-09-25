function imarisApp = viewMovieCurvatureImaris(movieData,varargin)

%Uses imaris to visualize the surface curvature, fluorescence, and cortical
%fluorescence samples simultaneously.


%Hunter Elliott
%2/2013




%% ---------- Parameters ---------- %%

iProcChan = 1;


nBins = 100;%Number of bins in histogram for correlation exploration

%% ------ Input ----- %%

nChan = numel(movieData.channels_);
chanNames = movieData.getChannelNames;

ip = inputParser;
ip.addParamValue('ChannelIndex',1:nChan,@(x)(numel(x) >= 1 && all(isposint(x))));%Raw fluorescence to show
ip.addParamValue('CurveTypeIndex',[],@(x)(all(isposint(x))));%Curvature types to show as channel
ip.addParamValue('SampleChannelIndex',[],@(x)(all(isposint(x))));%Cortical samples to show
ip.addParamValue('SampleTypeIndex',[],@(x)(numel(x) <= 1 && all(isposint(x))));%Intensity sample type to show as channel
ip.addParamValue('ExploreCorrelation',true,@(x)(numel(x) == 1 &&  islogical(x)));%Intensity sample type to show as channel
ip.parse(varargin{:});
p = ip.Results;

pixXY = movieData.pixelSize_;
pixZ = movieData.zSpacing_;           

%Possible curvature/ intensity measure fields
[curvTypes,curvNames,curvUnits,curvConv] = getCurveTypeFields(pixXY,true);%Use microns because imaris doesn't handle displaying very small numbers well
[intTypes,intNames] = getIntTypeFields;


if isempty(p.CurveTypeIndex)
    p.CurveTypeIndex = listdlg('ListString',curvNames,'PromptString','Curvature samples to display:');    
end

if isempty(p.SampleChannelIndex)
    p.SampleChannelIndex = listdlg('ListString',chanNames,'PromptString','Intensity samples to display:');    
end

if isempty(p.SampleTypeIndex)
    p.SampleTypeIndex = listdlg('ListString',intNames,'PromptString','Intensity samples type to display:','SelectionMode','single');    
end


%% ------------------------ Init ------------------------ %%

if p.ExploreCorrelation && (numel(p.SampleChannelIndex)~= 1  || numel(p.CurveTypeIndex) ~= 1)
    warning('Correlation exploration can only be run with one intensity sample and one curvature type! Disabling exploration')
    p.ExploreCorrelation = false;
end

if p.ExploreCorrelation
    expChans = 1;
else
    expChans = 0;
end



imNames = movieData.getImageFileNames(p.ChannelIndex);
imDirs = movieData.getChannelPaths(p.ChannelIndex);

iceConn = IceImarisConnector;
iceConn.startImaris;
imarisApp = iceConn.mImarisApplication;

%Create a blank scene
imarisScene = imarisApp.GetFactory.CreateDataContainer;
imarisApp.SetSurpassScene(imarisScene)

%Add lighting and frame objects to scene
imarisScene.AddChild(imarisApp.GetFactory.CreateLightSource,0); %add the light to the scene
imarisScene.AddChild(imarisApp.GetFactory.CreateFrame,0); %add the frame to the scene


iMI = movieData.getProcessIndex('MaskedIntensity3DProcess',1,1);


nFrames = movieData.nFrames_;

iMG = movieData.getProcessIndex('MaskGeometry3DProcess',1,1);

if isempty(iMI) || isempty(iMG)
    error('Masked intensity analysis and geometry analysis must both have been run!')
end

%imSize = [movieData.imSize_,movieData.nSlices_];

curvTypeToShow = curvTypes(p.CurveTypeIndex);
curvNameToShow = curvNames(p.CurveTypeIndex);
nCurvType = numel(curvTypeToShow);                 

sampTypeToShow = intTypes{p.SampleTypeIndex};
sampNameToShow = intNames{p.SampleTypeIndex};
nSampShow = numel(p.SampleChannelIndex);

nChanShow = numel(p.ChannelIndex);


nChanTot = nChanShow + nCurvType + nSampShow + expChans;

%Setup colors for various real and synthetic channels
chanCols = hsv(nChanTot);


%String for setting frame times. 
oneFrame= movieData.timeInterval_/(24*60*60); %Fraction of day corresponding to one frame sec, for datestr.m use
yearString = '2000-01-01 ';%The year/date doesn't seem to matter, so I just use this generic one
msString = '.000'; %The miliseconds portion of the date string




%% ------------------ Loading & Display ---------------------- %%

intAna = movieData.processes_{iMI}.loadChannelOutput(iProcChan);

for iFrame = 1:nFrames

    
    %------- Raw Fluorescence Images ----- %
    
    for iChan = 1:nChanShow
            
        currIm = stackRead([imDirs{iChan} filesep imNames{iChan}{iFrame}]);        
        
        currIm = make3DImageVoxelsSymmetric(currIm,pixXY,pixZ);
        
        if iFrame == 1 && iChan == 1
        
            %Initialize the volume data, using the resized image dimensions
            imSize = size(currIm);
            imSizeNM = imSize * pixXY;
                    
            volData = imarisApp.GetFactory.CreateDataSet;
            classDataSet=Imaris.tType.eTypeFloat; %Use float so we can show curvature too
            volData.Create(classDataSet ,imSize(1),imSize(2),imSize(3),nChanTot,nFrames); 


            %Set the pixel sizes
            volData.SetExtendMinX(0);
            volData.SetExtendMinY(0);
            volData.SetExtendMinZ(0);
            volData.SetExtendMaxX(imSizeNM(1));
            volData.SetExtendMaxY(imSizeNM(2));
            volData.SetExtendMaxZ(imSizeNM(3));
            volData.SetUnit('nm'); %Set units to nanometers
        
        end

        volData.SetDataVolumeAs1DArrayFloats(single(currIm(:)),iChan-1,iFrame-1);
        
        volData.SetChannelColorRGBA(iChan-1,iceConn.mapRgbaVectorToScalar([chanCols(iChan,1) chanCols(iChan,2) chanCols(iChan,3) 0]))            
        volData.SetChannelRange(iChan-1,min(currIm(:)),max(currIm(:)))
        volData.SetChannelName(iChan-1,chanNames{p.ChannelIndex(iChan)})
        
                    
    end
    
    % ------ Cortical Samples ------ %
    
    for j = 1:nSampShow
        sampIm = zeros(imSize);
        sampIm(intAna.branchProfiles(iFrame).surfPixInd) =  intAna.branchProfiles(iFrame).(sampTypeToShow)(:,p.SampleChannelIndex(j));

        volData.SetDataVolumeAs1DArrayFloats(single(sampIm(:)),nChanShow+j-1,iFrame-1);        
        volData.SetChannelColorRGBA(nChanShow+j-1,iceConn.mapRgbaVectorToScalar([chanCols(nChanShow+j,1) chanCols(nChanShow+j,2) chanCols(nChanShow+j,3) 0]))            
        volData.SetChannelRange(nChanShow+j-1,min(sampIm(:)),max(sampIm(:)))
        volData.SetChannelName(nChanShow+j-1,[chanNames{p.SampleChannelIndex(j)} ' ' sampNameToShow]) 
    end        
    
    
    %------- Curvature Samples ----- %       
    
    
    %Create the synthetic images showing the various curvatures
    for j = 1:nCurvType
        curvIm = ones(imSize) * (min(real(intAna.branchProfiles(iFrame).(curvTypeToShow{j})) .* curvConv(p.CurveTypeIndex(j)))-1);%So we can easily remove these from the display
        curvIm(intAna.branchProfiles(iFrame).surfPixInd) = real(intAna.branchProfiles(iFrame).(curvTypeToShow{j})) .* curvConv(p.CurveTypeIndex(j));
        
        
        volData.SetDataVolumeAs1DArrayFloats(single(curvIm(:)),nChanShow+nSampShow+j-1,iFrame-1);
        
        volData.SetChannelColorRGBA(nChanShow+nSampShow+j-1,iceConn.mapRgbaVectorToScalar([chanCols(nChanShow+nSampShow+j,1) chanCols(nChanShow+nSampShow+j,2) chanCols(nChanShow+nSampShow+j,3) 0]))            
        volData.SetChannelRange(nChanShow+nSampShow+j-1,min(curvIm(:)),max(curvIm(:)))
        volData.SetChannelName(nChanShow+nSampShow+j-1,curvNameToShow{j})                
    end
    
    %Get the time string for this frame
    if iFrame == 1
        %Datestr.m doesn't return the last portion if it's all zeros...
        secString = '00:00:00';        
    else    
        %Use datestr to convert to hr:min:sec
        secString = datestr(1+(iFrame-1)*oneFrame);
        secString = secString(13:end);
    end
    tString = [yearString secString msString];
    volData.SetTimePoint(iFrame-1,tString);
end

imarisApp.SetDataSet(volData)
volOb = imarisApp.GetFactory.CreateVolume;
imarisApp.GetSurpassScene.AddChild(volOb,-1)
imarisApp.GetSurpassCamera.Fit;


%% ----------- Correlation Exploration -------- %%
%NOTE: Only supports one curvature type and intensity type at a time.

if p.ExploreCorrelation
    
    % ----- Figure Display, Range Selection ----- %%
    for iFrame = 1:nFrames
        allInt{iFrame} = intAna.branchProfiles(iFrame).(sampTypeToShow)(:,p.SampleChannelIndex);
        allCurv{iFrame} = real(intAna.branchProfiles(iFrame).(curvTypeToShow{1})) .* curvConv(p.CurveTypeIndex);
    end
    allInt = vertcat(allInt{:});
    allCurv = vertcat(allCurv{:});    
    [N,C] = hist3([allCurv allInt],[nBins nBins]);
    histFig = figure;    
    histHan = imagesc(C{2},C{1},log10(N));
    hold on
    axHan = get(histHan,'Parent');
    xlabel([chanNames{p.SampleChannelIndex(j)} ' ' sampNameToShow])
    ylabel([curvNameToShow{1} ' ' curvUnits{p.CurveTypeIndex}])
    
    axis xy
    bPressed = 'Update';
    
    
    while any(strcmp(bPressed,{'Update','Save Figures'}))
        
        
        if strcmp(bPressed,'Update')
            if exist('rangeRect','var')
                rangeRect.delete;
            end
            rangeRect = imrect(axHan);
        end
        currRect = rangeRect.getPosition;
        %Because it returns the rectangle in this stupid fucking format....
        currRectCoord = [currRect(1) currRect(2) ; 
                        currRect(1) + currRect(3) currRect(2) ;
                        currRect(1) + currRect(3) currRect(2) + currRect(4);
                        currRect(1) currRect(2) + currRect(4)];
           
        if ~exist('showRect','var') || ~ishandle(showRect)
            showRect = fill(currRectCoord(:,1),currRectCoord(:,2),'w','FaceAlpha',.3,'EdgeColor','k');
        else
            set(showRect,'Vertices',currRectCoord)
        end
        
        for iFrame = 1:nFrames
            currSelect = zeros(imSize);
            isSelected = intAna.branchProfiles(iFrame).(sampTypeToShow)(:,p.SampleChannelIndex) >= currRect(1) & ...
                         intAna.branchProfiles(iFrame).(sampTypeToShow)(:,p.SampleChannelIndex) <= currRect(1) + currRect(3) & ...
                         real(intAna.branchProfiles(iFrame).(curvTypeToShow{1})) .* curvConv(p.CurveTypeIndex) >= currRect(2) & ...
                         real(intAna.branchProfiles(iFrame).(curvTypeToShow{1})) .* curvConv(p.CurveTypeIndex) <= currRect(2) + currRect(4);                     
            currSelect(intAna.branchProfiles(iFrame).surfPixInd(isSelected)) = 1;
            
            
            volData.SetDataVolumeAs1DArrayFloats(single(currSelect(:)),nChanShow+nSampShow+nCurvType,iFrame-1);
        
            volData.SetChannelColorRGBA(nChanShow+nSampShow+nCurvType,iceConn.mapRgbaVectorToScalar([chanCols(nChanShow+nSampShow+nCurvType+1,1) chanCols(nChanShow+nSampShow+nCurvType+1,2) chanCols(nChanShow+nSampShow+nCurvType+1,3) 0]))            
            volData.SetChannelRange(nChanShow+nSampShow+nCurvType,0, 2);
            volData.SetChannelName(nChanShow+nSampShow+nCurvType,'Current Selection');
            
            
        end
        
        bPressed = questdlg('What now?','Correlation Exploration','Update','Save Figures','Cancel','Cancel');
        
        if strcmp(bPressed,'Save Figures')
            [fileName,filePath] = uiputfile;
            saveThatShit([fileName ' histogram'],filePath)
            imarisApp.SaveSnapShot([filePath filesep fileName ' image.tif']);
        end
        
    end    
    
end

