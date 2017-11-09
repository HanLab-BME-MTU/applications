function [CCPstruct,nascentAdhInfo,focalAdhInfo,indivFAArea,shiftedFAarea,localCCPDensity,localCCPDensityNA,distToEdgeFA,distToEdgeNA] = getDistanceCCPToAdhesion(pathForTheMovieDataFile,bandwidth,iCCP,iPax,radiiMicron)

if nargin==1
    bandwidth = 5;
end
if nargin<3
    iCCP = 1; %assumed
end
if nargin<3
    iPax=3;
end
%% Load the MovieData
if isa(pathForTheMovieDataFile, 'MovieData')
    movieData = pathForTheMovieDataFile;
else
    movieDataPath = [pathForTheMovieDataFile '/movieData.mat'];
    movieData = MovieData.load(movieDataPath,false);
end
% run detectMovieNascentAdhesion
plotGraph = false;
[nascentAdhInfo,focalAdhInfo] = detectMovieNascentAdhesion(movieData,bandwidth,iPax,plotGraph);

%% --------------- Set up ---------------%%
% nChan=numel(movieData.channels_);
% psfSigma = 1.2; %hard coded for nascent adhesion
psfSigma = getGaussianPSFsigma(movieData.numAperture_, 1, movieData.pixelSize_*1e-9, movieData.getChannel(iCCP).emissionWavelength_*1e-9);

% Set up the output directories
outputFilePath = [movieData.outputDirectory_ filesep 'CCP_distanceAnalysis'];
imgPath = [outputFilePath filesep 'imgs'];
dataPath = [outputFilePath filesep 'data'];
tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];

if ~exist(figPath,'dir')
    mkdir(imgPath);
    mkdir(dataPath);
    mkdir(tifPath);
    mkdir(figPath);
end

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting CCPs...')

roiMask = movieData.getROIMask;
CCPstruct(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
    'amp',[], 'numberCCP',[],'cellArea',[],'CCPdensity',[],'distToAdh',[],'closeToNA',[],'idNA',[],'idFA',[]);

h1=figure;
jformat = ['%.' '3' 'd'];
% Changed it for isometric detection for nascent adhesion detection
pixSize = movieData.pixelSize_;
% iCCP = 1;
disp(['CCP channel was assumed to be in channel ' num2str(iCCP) '.'])
disp('Results will be saved under:')
disp(outputFilePath);

% radiiMicron=[0.05 0.1 0.2 0.5 1 2 4];
radii=radiiMicron*1000/pixSize;
numRadii=length(radii);
localCCPDensity=cell(1,numRadii);
shiftedFAarea = cell(1,numRadii);
localCCPDensityNA=cell(1,numRadii);
distToEdgeNA=cell(1,numRadii);
distToEdgeFA=cell(1,numRadii);
% localAdhDensity=cell(1,numRadii);

for j=1:movieData.nFrames_
    Iccp=double(movieData.channels_(iCCP).loadImage(j));
    maskProc = movieData.getProcess(movieData.getProcessIndex('MaskRefinementProcess'));
    % if there are masks for more than one channels, combine them.
    if sum(maskProc.checkChannelOutput)>1
        %Combine the the multiple masks to one
        maskEach = arrayfun(@(x) maskProc.loadChannelOutput(x,j),find(maskProc.checkChannelOutput),'UniformOutput',false);
        maskAll=reshape(cell2mat(maskEach),size(Iccp,1),size(Iccp,2),[]);
        mask = any(maskAll,3);
        % Select only one chunk of the mask
        %Label all objects in the mask
        labelMask = bwlabel(mask);
        %Get their area
        obAreas = regionprops(labelMask,'Area');       %#ok<MRPBW>
        p.ObjectNumber =1;
        %First, check that there are objects to remove
        if length(obAreas) > p.ObjectNumber 
            obAreas = [obAreas.Area];
            %Sort by area
            [dummy,iSort] = sort(obAreas,'descend'); %#ok<ASGLU>
            %Keep only the largest requested number
            mask = false(size(mask));
            for i = 1:p.ObjectNumber
                mask = mask | labelMask == iSort(i);
            end
        end
    elseif length(iChans)==1
        mask = maskProc.loadChannelOutput(iPax,j); % 1 is CCP channel
    end
    ultimateMask = roiMask(:,:,j) & mask;
    pstruct = pointSourceDetection(Iccp, psfSigma, 'Alpha',0.05,'Mask', ultimateMask);
    % filter out points where overall intensity is in the noise level
    pixelIntenMargin = Iccp(~ultimateMask);
    maxIntBg=quantile(pixelIntenMargin,0.99);
    psInt = pstruct.A+pstruct.c;
    idxSigCCP = psInt>maxIntBg;
    
    CCPstruct(j).xCoord = [round(pstruct.x(idxSigCCP)'), round(pstruct.x_pstd(idxSigCCP)')];
    CCPstruct(j).yCoord = [round(pstruct.y(idxSigCCP)'), round(pstruct.y_pstd(idxSigCCP)')];
    CCPstruct(j).amp = [pstruct.A(idxSigCCP)', pstruct.A_pstd(idxSigCCP)'];
    CCPstruct(j).numberCCP = length(pstruct.x(idxSigCCP));
    CCPstruct(j).cellArea = sum(ultimateMask(:))*(pixSize/1000)^2; % in um^2
    CCPstruct(j).CCPdensity = CCPstruct(j).numberCCP/CCPstruct(j).cellArea; % number per um2

    % Closest distance analysis: distance to NAs
    [idxNAfromCCP, distToNA] = KDTreeClosestPoint([nascentAdhInfo(j).xCoord(:,1), nascentAdhInfo(j).yCoord(:,1)],...
                                                    [CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1)]);
    CCPstruct(j).distToAdh = distToNA;
    CCPstruct(j).closeToNA = true(size(distToNA)); %true means it is close to NA than FA (false: close to FA more than NA)
    CCPstruct(j).idNA = idxNAfromCCP;
    
    % distance to FAs
    allBoundFAs=cell2mat(focalAdhInfo(j).boundFA);
    integratedFApoints = [focalAdhInfo(j).xCoord(:,1), focalAdhInfo(j).yCoord(:,1);...
                                                                             allBoundFAs(:,2), allBoundFAs(:,1)];
    [idxFAfromCCP, distToFA] = KDTreeClosestPoint(integratedFApoints,...
                                                    [CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1)]);
    CCPstruct(j).distToAdh(distToFA<distToNA) = distToFA(distToFA<distToNA);
    CCPstruct(j).closeToNA(distToFA<distToNA) = false; %(false: close to FA more than NA)
    CCPstruct(j).idFA = idxFAfromCCP;
    
    % plotting detected adhesions
    Ipax=double(movieData.channels_(iPax).loadImage(j));
    dIpax =double(Ipax)/max(max(Ipax)); 
    dIccp = double(Iccp)/max(max(Iccp));
    combI(:,:,1) = dIpax;
    combI(:,:,2) = dIccp;
    combI(:,:,3) = zeros(size(dIpax));
    imshow(combI,[]); hold on
    plot(CCPstruct(j).xCoord(:,1),CCPstruct(j).yCoord(:,1),'go')
    plot(nascentAdhInfo.xCoord(:,1),nascentAdhInfo.yCoord(:,1),'y.')
    nFA = focalAdhInfo.numberFA;
    for k = 1:nFA
        adhBoundary = focalAdhInfo.boundFA{k};
        plot(adhBoundary(:,2), adhBoundary(:,1), 'Color','r', 'LineWidth', 0.5) %adhesion boundary
    end
    %% Local density with peripheral shift
%         queryPtsNA = [curNAstruct.xCoord(:,1) curNAstruct.yCoord(:,1)];
    queryPtsFA = (cellfun(@(x) [x(:,2), x(:,1)],focalAdhInfo.boundFA,'UniformOutput',false));
    inputPoints = [CCPstruct.xCoord(:,1) CCPstruct.yCoord(:,1)];
    indivFAArea= []; %focalAdhInfo.area';
    % per each x,y, get the density
    queryPtsNA = [nascentAdhInfo.xCoord(:,1) nascentAdhInfo.yCoord(:,1)];

    for jj=numRadii:-1:1
        ccpDensityLabelsFA=zeros(size(focalAdhInfo.labelFA));
        % Nacent adhesion
%         curAreaNA=pi*radiiMicron(jj)^2; % in um2
        for pp=1:length(queryPtsNA)
            eachNAMask=false(size(mask));
            eachNAMask(sub2ind(size(mask),queryPtsNA(pp,2), queryPtsNA(pp,1)))=true;
            se = strel('disk',round(radii(jj))+1,0);
            eachNAMask=imdilate(eachNAMask,se);
            legitNAMask=mask & eachNAMask;
            insideIdx = maskVectors(inputPoints(:,1),inputPoints(:,2),legitNAMask);
            numCCParoundNA = sum(insideIdx);
            curAreaNA = sum(legitNAMask(:))*((pixSize/1000)^2); % in um2
            if curAreaNA>0
                ccpDensityLabelsFA(legitNAMask)=numCCParoundNA/curAreaNA;
                localCCPDensityNA{1,jj}=[localCCPDensityNA{1,jj}; numCCParoundNA/curAreaNA];

                bwdistMask = bwdist(~mask);
                curDistToEdgeNA = bwdistMask(round(queryPtsNA(pp,2)), round(queryPtsNA(pp,1)));
                distToEdgeNA{1,jj}=[distToEdgeNA{1,jj}; curDistToEdgeNA];
            end
        end

        for pp=1:numel(queryPtsFA)
            % CCP density around adhesions
            curFAMask = focalAdhInfo.labelFA==pp;
            % Expand this mask with radius
            curFAMaskExpand = bwmorph(curFAMask,'dilate', round(radii(jj)));
            curFAMaskShrunk = bwmorph(curFAMask,'erode', 1);
            curSearchMask = curFAMaskExpand & ~curFAMaskShrunk;
%             curSearchMask = curFAMaskExpand;
            % See how many ccp points are in this search mask
            insideIdx = maskVectors(inputPoints(:,1),inputPoints(:,2),curFAMaskExpand);
            numCCParound = sum(insideIdx);
            % Area around: should be excluding outside the cell mask
            legitFAMask=mask & curSearchMask;
            curArea = sum(legitFAMask(:))*((pixSize/1000)^2); % in um2
            repFAMask=mask & curFAMaskExpand;
            if curArea>0
                localCCPDensity{1,jj}=[localCCPDensity{1,jj}; numCCParound/curArea];
    %             % Adhesion density around CCP
    %             [idxAdh]= KDTreeBallQuery(queryPtsNA,inputPoints,radii(jj));
    %             numAdhAround=cellfun(@numel,idxAdh);
    %             curArea=pi*radiiMicron(jj)^2; % in um2
    %             localAdhDensity{ii,jj}=[localAdhDensity{ii,jj}; numAdhAround/curArea];
                ccpDensityLabelsFA(repFAMask)=numCCParound/curArea;

                bwdistMask = bwdist(~mask);
                curFACenter = regionprops(curFAMask,'Centroid');
                curDistToEdgeFA = bwdistMask(round(curFACenter.Centroid(2)), round(curFACenter.Centroid(1)));
                distToEdgeFA{1,jj}=[distToEdgeFA{1,jj}; curDistToEdgeFA];
                curFAMask=mask & curFAMask;
                curFAArea = sum(curFAMask(:))*((pixSize/1000)^2); % in um2
                curFAMaskExpand=mask & curFAMaskExpand;
                curFAAreaExpand = sum(curFAMaskExpand(:))*((pixSize/1000)^2); % in um2
                indivFAArea=[indivFAArea; curFAArea];
                shiftedFAarea{1,jj}=[shiftedFAarea{1,jj}; curFAAreaExpand];
            end
        end
        if jj==numRadii
            h3=figure;
        else
            figure(h3)
        end
        imshow(ccpDensityLabelsFA,[0 4]), colormap jet
        hC=colorbar;
        hC.Label.String='Local CCP Density (#/um^2)';
        ax = gca;
        axpos = ax.Position;
        cpos = hC.Position;
        cpos(3) = 0.2*cpos(3);
        cpos(4) = 0.5*cpos(4);
        hC.Position = cpos;
        ax.Position = axpos;
        hC.Label.Color='w';
        hC.Color='w';
        hC.FontSize=6;

        try
            hgexport(h3,strcat(tifPath,'/imgLocalCCPDensityAdh',num2str(j,jformat), 'R', num2str(radiiMicron(jj)*1000),'nm'),hgexport('factorystyle'),'Format','tiff')
            hgsave(h3,strcat(figPath,'/imgLocalCCPDensityAdh',num2str(j,jformat), 'R', num2str(radiiMicron(jj))*1000,'nm'),'-v7.3')
        catch
            disp('I''m saving no figures...')
        end
        hold off
    end

    figure(h1), hold on
    plot([CCPstruct(j).xCoord(CCPstruct(j).closeToNA,1) nascentAdhInfo(j).xCoord(CCPstruct(j).idNA(CCPstruct(j).closeToNA),1)]',...
        [CCPstruct(j).yCoord(CCPstruct(j).closeToNA,1) nascentAdhInfo(j).yCoord(CCPstruct(j).idNA(CCPstruct(j).closeToNA),1)]','y')
    plot([CCPstruct(j).xCoord(~CCPstruct(j).closeToNA,1) integratedFApoints(CCPstruct(j).idFA(~CCPstruct(j).closeToNA),1)]',...
        [CCPstruct(j).yCoord(~CCPstruct(j).closeToNA,1) integratedFApoints(CCPstruct(j).idFA(~CCPstruct(j).closeToNA),2)]','r')
    try
        hgexport(h1,strcat(tifPath,'/imgCCPtoAdh',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h1,strcat(figPath,'/imgCCPtoAdh',num2str(j,jformat)),'-v7.3')
    catch
        try
            hgexport(h1,strcat(tifPath,'/imgCCPtoAdh_',num2str(j,jformat)),hgexport('factorystyle'),'Format','tiff')
            hgsave(h1,strcat(figPath,'/imgCCPtoAdh_',num2str(j,jformat)),'-v7.3')
        catch
            disp('I''m saving no figures...')
        end
    end
    hold off
    % Show local CCP density for all cell interior
%     [Imask,Jmask] = ind2sub(size(mask),find(mask(:)));
%     queryPtsAllCell = [Jmask,Imask];
%     inputPoints = [CCPstruct.xCoord(:,1) CCPstruct.yCoord(:,1)];
%     jj=numRadii;
%     % CCP density around each pixel
%     [idx]= KDTreeBallQuery(inputPoints,queryPtsAllCell,radii(jj));
%     numCCParound=cellfun(@numel,idx);
%     [idxCellPointsAround]= KDTreeBallQuery(queryPtsAllCell,queryPtsAllCell,radii(jj));
%     curArea = cellfun(@numel,idxCellPointsAround)*((pixSize/1000)^2); % in um2
%     % Area around: should be excluding outside the cell mask
%     localCCPDensityPerPixel = numCCParound./curArea;
%     % Assign density to pixels
%     densityMap = zeros(size(mask));
%     densityMap(mask)=localCCPDensityPerPixel;
%     if jj==1
%         h2=figure;
%     else
%         figure(h2)
%     end
%     minNextDen = min(densityMap(densityMap>0));
%     maxDensity = max(densityMap(:));
%     imshow(densityMap,[minNextDen-5 maxDensity]), colormap jet
%     colorbar
%     hold on
%     plot(nascentAdhInfo.xCoord(:,1),nascentAdhInfo.yCoord(:,1),'r.')
%     for k = 1:nFA
%         adhBoundary = focalAdhInfo.boundFA{k};
%         plot(adhBoundary(:,2), adhBoundary(:,1), 'Color','r', 'LineWidth', 0.5) %adhesion boundary
%     end
% 
%     try
%         hgexport(h2,strcat(tifPath,'/imgLocalCCPDensity',num2str(j,jformat), 'R', num2str(radiiMicron(jj)),hgexport('factorystyle'),'Format','tiff')
%         hgsave(h2,strcat(figPath,'/imgLocalCCPDensity',num2str(j,jformat), 'R', num2str(radiiMicron(jj)),'-v7.3')
%     catch
%         disp('I''m saving no figures...')
%     end
    %
    % CCP Density plot for individual FAs
    % Make labels for ccp density per each FA
end
close(h1)
close(h3)
try
    save([dataPath filesep 'AdhInfo.mat'],'nascentAdhInfo','focalAdhInfo','CCPstruct','-v7.3');
catch
    try
        save([dataPath filesep 'AdhInfo_.mat'],'nascentAdhInfo','focalAdhInfo','CCPstruct','-v7.3');
    catch
        disp('I''m saving no data...')
    end
end

disp('Finished detecting and linking objects...')