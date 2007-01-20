function [resultList,plotData] = bilobeSpindles
%BILOBESPINDLES projects kinetochore signals onto the spindle axis in two-color images
%
% SYNOPSIS: bilobeSpindles
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 27-Oct-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spbCorrection_1 = [-0.0124 0.1647 0.2546]; % from diploid 1-7. CFP->GFP

spbCorrection_2 = [-0.0480 0.3796 0.3398]; % January 2007


% turn off exposure time/nd warnings
warningState = warning;
warning('off','R3DREADHEADER:exposureTimeChanged');
warning('off','R3DREADHEADER:ndFilterChanged');

plotGrid = false; % imaris will plot the 3D-grid
plotInt = false; % imaris will plot spb and ndc80
redoAll = false; % if false, code will not redo analysis

plotSigma = true; % plot also sigma of the intensities
useMax = true; % plot maximum intensity projection

% search files
cdBiodata;
[fileList,tokens]=searchFiles('(\w+)_R3D\.dv','log|DIC','ask');

% select files
selection = listSelectGUI(tokens);

if isempty(selection)
    resultList = [];
    plotData = [];
    return
end

% only keep selection
fileList = fileList(selection,:);
tokens = tokens(selection);

nData = length(selection);

% find what is common between all the names
overallName = tokens{1};
goodIdx = 1:length(overallName);
for i=2:nData
    goodIdx(goodIdx > length(tokens{i})) = [];
    goodIdx = find(overallName(1:length(goodIdx))==tokens{i}(goodIdx));
    overallName = overallName(goodIdx);
end

nameList(1:nData) = struct('rawMovieName',[],'filteredMovieName',[],'slistName',[],...
    'dataPropertiesName',[],'testRatiosName',[]);

resultList(1:nData) = struct('interpImg',[],'maxProj',[],'meanProj',[],...
    'e_spb',[],'e_perp',[],'e_3',[],'n_spb',[]);

status(1:nData,1) = -1;

% loop through selection and perform extraction
for iData = 1:nData



    % make directory if necessary
    if ~any(findstr(fileList{iData,2},tokens{iData}));
        mkdir(fileList{iData,2},tokens{iData})
        % find everything there is to move
        moveList = searchFiles(tokens{iData},[],fileList{iData,2});
        for i=1:length(moveList)
            movefile(fullfile(fileList{iData,2},moveList{i,1}),...
                fullfile(fileList{iData,2},tokens{iData},moveList{i,1}));
        end
        % update file path
        fileList{iData,2} = fullfile(fileList{iData,2},tokens{iData});
    end

    % set names - include directories already
    nameList(iData).rawMovieName = fullfile(fileList{iData,2},fileList{iData,1});
    nameList(iData).deconMovieName = ...
        [nameList(iData).rawMovieName(1:end-3),...
        '_D3D',nameList(iData).rawMovieName(end-2:end)];
    nameList(iData).filteredMovieName = ...
        [fileList{iData,2},filesep,'filtered_',tokens{iData},'.fim'];
    nameList(iData).movieHeaderName = ...
        [fileList{iData,2},filesep,'movieHeader_',tokens{iData},'.mat'];
    nameList(iData).filteredMovieName = ...
        [fileList{iData,2},filesep,'filtered_',tokens{iData},'.fim'];
    nameList(iData).slistName = ...
        [fileList{iData,2},filesep,'slist_',tokens{iData},'.mat'];
    nameList(iData).idlistName = ...
        [fileList{iData,2},filesep,'idlist_',tokens{iData},'.mat'];
    nameList(iData).idlist_LName = ...
        [fileList{iData,2},filesep,'idlist_L_',tokens{iData},'.mat'];
    nameList(iData).testRatiosName = ...
        [fileList{iData,2},filesep,'testRatios_',tokens{iData},'.mat'];
    nameList(iData).dataPropertiesName = ...
        [fileList{iData,2},filesep,'dataProperties_',tokens{iData},'.mat'];
    nameList(iData).dirName = fileList{iData,2};

    % check status

    if exist(nameList(iData).deconMovieName,'file')
        status(iData) = 0;
        if ~redoAll
            if exist(nameList(iData).filteredMovieName,'file')
                status(iData) = 1; % filtered
                if exist(nameList(iData).slistName,'file')
                    status(iData) = 2; % detected
                    if exist(nameList(iData).idlistName,'file')
                        status(iData) = 3; % done
                    end
                end
            end
        end
    else
        disp(sprintf(...
            'Deconvolved movie %s not found!',...
            nameList(iData).deconMovieName));
    end

    % check wheter we will do anything at all
    if status(iData) > -1

        if status(iData) == 0
            % load movie - use low-level r3dread
            rawMovie = r3dread(nameList(iData).rawMovieName);
            movieHeader = readr3dheader(nameList(iData).rawMovieName);

            % save movieHeader
            save(nameList(iData).movieHeaderName,'movieHeader');

            % select spb channel
            % spb channel is the first channel
            rawMovie = rawMovie(:,:,:,1);

            % make data properties for the first channel
            dataProperties = defaultDataProperties(movieHeader);
            dataProperties.waveIdx = 1; % spb channel
            % find filter parameters - only for the interesting wavelength!
            [FT_XY, FT_Z] = calcFilterParms(...
                dataProperties.WVL(dataProperties.waveIdx),...
                dataProperties.NA,1.51,'gauss',...
                dataProperties.sigmaCorrection, ...
                [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
            patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
            dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
            dataProperties.FT_SIGMA = [FT_XY,FT_XY,FT_Z];
            dataProperties.MAXSPOTS = 2;
            dataProperties.fitNPlusOne = false;
            dataProperties.amplitudeCutoff = 0; % set once we know the images better

            % save data Properties
            save(nameList(iData).dataPropertiesName,'dataProperties');


            % filter movie
            filteredMovie = filtermovie(rawMovie,dataProperties.FILTERPRM);

            % save filtered movie
            writemat(nameList(iData).filteredMovieName,...
                filteredMovie);

            status(iData) = 1;
        end

        if status(iData) == 1

            % detect two spots


            [slist, dataProperties, testRatios] = detectSpots(...
                nameList(iData).rawMovieName, ...
                nameList(iData).filteredMovieName, ...
                dataProperties,2); %#ok<NASGU>

            % save results
            save(nameList(iData).slistName,'slist');
            save(nameList(iData).dataPropertiesName,'dataProperties')
            save(nameList(iData).testRatiosName,'testRatios');

            status(iData) = 2;
        end
        if status(iData) == 2
            % make idlist, save
            idlist = linker(slist,dataProperties,1);
            save(nameList(iData).idlistName,'idlist');
            status(iData) = 3;
        end
    end % if analyze at all
end % loop iData

% loop through data, check whether we have the proper number of spots
for iData=1:nData
    if status(iData) == 3
        % load idlist
        load(nameList(iData).idlistName);
        % check for number of spots
        nSpots = size(idlist(1).linklist,1);

        if nSpots > 2
            filteredMovie = readmat(nameList(iData).filteredMovieName);
            load(nameList(iData).dataPropertiesName);
            lh = LG_loadAllFromOutside(filteredMovie,nameList(iData).dirName,...
                [],dataProperties,idlist,'idlist');
            uiwait(lh)
            idlist_L = LG_readIdlistFromOutside;
        else
            idlist_L = idlist;
        end

        if ~isempty(idlist_L)
            status(iData) = 4;
        end

        % save idlist_L
        save(nameList(iData).idlist_LName,'idlist_L');
    end
end


% loop through data, read intensities
for iData = 1:nData
    if status(iData) == 4
        % load idlist_L
        load(nameList(iData).idlist_LName);

        % select spbCorrection
        d=dir(nameList(iData).rawMovieName);
        if datenum(d.date) < datenum('1-Jan-2007 00:00:01')
            spbCorrection = spbCorrection_1;
        else
            spbCorrection = spbCorrection_2;
        end

        % find coordinates of two spots - already correct for x/y. Sort
        % linklist so that we measure bottom to top, and that we can get the
        % spotIdx afterwards
        idlist_L.linklist = sortrows(idlist_L.linklist,11);
        spbCoord = idlist_L.linklist(:,[10,9,11]);

        if size(spbCoord,1) ~=2
            % not two spb tags found
            status(iData) = -1;
        else

            % find vector - point from bottom to top
            v_spb = spbCoord(2,:) - spbCoord(1,:);
            [n_spb,e_spb] = normList(v_spb);

            % calculate angle from plane
            spbVectorAngle = 90 - acos(e_spb(3))*180/pi;

            % make perpendicular vector parallel to xy
            e_perp = [-e_spb(2),e_spb(1),0]./sum(e_spb(1:2).^2);

            % complete coordinate system
            e_3 = cross(e_spb,e_perp);

            % convert spb coords in pixels
            load(nameList(iData).dataPropertiesName);
            pix2mu = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];
            spbCoordNdc = (spbCoord + repmat(spbCorrection,2,1)) ./ repmat(pix2mu,2,1);
            spbCoord = spbCoord ./repmat(pix2mu,2,1);

            % spb-vector: 1/24th of spindle length. Perpendicular vectors: 1 pix
            e_spb = (e_spb * n_spb/24)./pix2mu;
            e_perp = (e_perp * dataProperties.PIXELSIZE_XY)./pix2mu;
            e_3 = (e_3 * dataProperties.PIXELSIZE_XY)./pix2mu;

            % load image - deconvolved movie
            rawMovie = r3dread(nameList(iData).deconMovieName);

            % arbitrary grid: 15 pix perpendicular, 24 bins parallel to spindle
            % axis (Gardner)
            % -> use a few more : 21 perp, bins from -3:27
            xLabels = -2.5:26.5;
            [xGrid,yGrid,zGrid] = arbitraryGrid(e_spb, e_perp, e_3, spbCoordNdc(1,:),...
                [xLabels(1) xLabels(end)], [-10 10], [-10 10], size(rawMovie));

            %     % show grid - (multiply by pix2mu for um)
            if plotGrid
                %             d1 = spbCoordNdc(1,:) + 15 * e_spb;
                %             d2 = spbCoordNdc(1,:) + 15 * e_perp;
                %             d3 = spbCoordNdc(1,:) + 15 * e_3;
                %pd(3).XYZ = [d1;d2;d3];pd(3).color = [0,0,1,0];pd(3).spotRadius =
                %1;
                pd(2).XYZ = [xGrid(:),yGrid(:),zGrid(:)];pd(2).color = [0,1,0,0];pd(2).spotRadius = 0.1;
                pd(1).XYZ = spbCoordNdc;pd(1).color = [1,0,0,0];pd(1).spotRadius = 1;
                imarisPlot3(pd,[],rawMovie(:,:,:,2))
            end
            % interpolate ndc80
            interpImg = interp3(rawMovie(:,:,:,2),yGrid,xGrid,zGrid,'linear',NaN);

            % read intensities from spb image
            [xGrid,yGrid,zGrid] = arbitraryGrid(e_spb, e_perp, e_3, spbCoord(1,:),...
                [xLabels(1) xLabels(end)], [-10 10], [-10 10], size(rawMovie));
            interpImgSpb = interp3(rawMovie(:,:,:,1),yGrid,xGrid,zGrid,'linear',NaN);




            % make maximum, average projections
            maxProj = nanmax(nanmax(interpImg,[],3),[],2)-min(interpImg(:));
            meanProj = nanmean(nanmean(interpImg,3),2)-min(interpImg(:));

            % store data

            resultList(iData).interpImg = interpImg;
            resultList(iData).interpImgSpb = interpImgSpb;
            resultList(iData).maxProj = maxProj;
            resultList(iData).meanProj = meanProj;
            resultList(iData).xLabels = xLabels;
            resultList(iData).e_spb = e_spb;
            resultList(iData).e_perp = e_perp;
            resultList(iData).e_3 = e_3;
            resultList(iData).n_spb = n_spb;
            resultList(iData).spbVectorAngle = spbVectorAngle;
            resultList(iData).name = nameList(iData).deconMovieName;

            % show two plots
            if plotInt
                imarisShowArray(cat(5,interpImg,interpImgSpb));
            end
        end
    end
end

% remove bad results
resultList(status==-1)=[];


% plotData: data needed for plotting
plotData(1:3) = struct('xall',[],'yall',[],'zall',[],'yTickLabels',[]);

% plot averages
figureNames = {'asymmetric','a, max=1','symmetric','s, max=1'};
for i=1:4

    boundaries = [0.9:0.025:1.9;1.1:0.025:2.1];
    meanBoundaries = mean(boundaries,1);
    nBoundaries = size(boundaries,2);
    nBins = length(xLabels)-4; %make compatible with bilobeDistribution
    xall = zeros(nBins,nBoundaries);
    yall=zeros(nBins,nBoundaries);
    zall=zeros(nBins,nBoundaries);
    zallS = zeros(nBins,nBoundaries);
    yTickLabels = cell(nBoundaries,1);
    nSpindles = zeros(nBoundaries,1);
    spindleLength = cat(1,resultList.n_spb);
    spindleIdx = cell(nBoundaries,1);

    switch i
        case {1,2}
            if useMax
                allMP = cat(2,resultList.maxProj);
            else
                allMP = cat(2,resultList.meanProj);
            end
        otherwise
            if useMax
                allMP = cat(2,resultList.maxProj);
            else
                allMP = cat(2,resultList.meanProj);
            end
            allMP = (allMP+allMP(end:-1:1,:))/2;
    end
    switch i
        case {4,2}
            allMP = allMP./repmat(max(allMP,[],1),size(allMP,1),1);
        case {1,3}
            allMP = allMP./repmat(sum(allMP,1),size(allMP,1),1);
    end

    % figure, hold on
    for ct = 1:nBoundaries,
        spindleIdx{ct} = (spindleLength>boundaries(1,ct) & spindleLength<boundaries(2,ct));
        if any(spindleIdx{ct})
            % make weighted average - points 0.1um away have no weight
            weights = (0.1-abs(spindleLength(spindleIdx{ct})-meanBoundaries(ct)))/0.1;
            [averageMP,dummy,sigmaMP] = weightedStats(allMP(:,spindleIdx{ct})',...
                weights,'w');
            %old: averageMP = mean(allMP(:,spindleIdx{ct}),2);
            switch i
                case {4,2}
                    sigmaMP = sigmaMP/max(averageMP);
                    averageMP = averageMP/max(averageMP);

                otherwise
                    sigmaMP = sigmaMP/sum(averageMP);
                    averageMP = averageMP/sum(averageMP);

            end
            %plot3(x,ct*ones(size(x)),z),
            nSpindles(ct) = sum(weights);
            %old: nSpindles(ct) = nnz(spindleIdx{ct});

            zall(:,ct)=averageMP(3:end-2);
            zallS(:,ct) = sigmaMP(3:end-2);
        else
            % there are no spindles this long
            zall(:,ct) = NaN;
            zallS(:,ct) = NaN;
        end
        xall(:,ct)=xLabels(3:end-2)/24;
        yall(:,ct)=meanBoundaries(ct);
        yTickLabels{ct}=sprintf('%1.1f/%1.2f', ...
            meanBoundaries(ct),nSpindles(ct));

    end

    % check zallS for inf
    zallS(~isfinite(zallS)) = NaN;

    figure('Name',[overallName,' ',figureNames{i}])
    if plotSigma
        ah = subplot(1,2,1);
    else
        ah = gca;
    end
    contourf(xall,yall,zall,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
    %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
    %     axis tight
    set(ah,'yTick',1:0.1:2,'yTickLabel',yTickLabels(1:4:end),'yGrid','on')

    switch i
        case {4,2}
            set(ah,'CLim',[0,1])
        otherwise
            set(ah,'CLim',[0,nanmax(zall(:))])
    end
    if ~plotSigma
        colorbar('peer',ah)
    end
    % add lines
    hold on
    for d=0.2:0.2:0.8
        line(d./meanBoundaries(meanBoundaries>=2*d),...
            meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
        line(1-d./meanBoundaries(meanBoundaries>=2*d),...
            meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
    end

    if plotSigma
        ah = subplot(1,2,2);
        contourf(xall,yall,zallS,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
        %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
        %     axis tight
        set(ah,'yTick',1:0.1:2,'yTickLabel',yTickLabels(1:4:end),'yGrid','on')

        switch i
            case {4,2}
                set(ah,'CLim',[0,1])
            otherwise
                set(ah,'CLim',[0,nanmax(zall(:))])
        end
        colorbar('peer',ah)
        % add lines
        hold on
        for d=0.2:0.2:0.8
            line(d./meanBoundaries(meanBoundaries>=2*d),...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
            line(1-d./meanBoundaries(meanBoundaries>=2*d),...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
        end
    end


    %     xlim([0,1])
    %     ylim([1,2])
    %view([0 90])
    plotData(i).xall = xall;
    plotData(i).yall = yall;
    plotData(i).zall = zall;
    plotData(i).zallS = zallS;
    plotData(i).yTickLabels = yTickLabels;
    plotData(i).spindleIdx = spindleIdx;


end

figure('Name',sprintf('%s individual; sum=1',overallName));
[sortedSpindleLength,sortIdx]=sort(spindleLength);
int = cat(2,resultList.maxProj);
int = int(:,sortIdx)';
int = int./repmat(sum(int,2),1,size(int,2));
intSymm = 0.5*(int + int(:,end:-1:1));
subplot(1,2,1),imshow([int,NaN(size(int,1),1),intSymm],[]),
colormap jet,
subplot(1,2,2),
plot(sortedSpindleLength,length(sortIdx):-1:1,'-+')

% loop to make histograms
stages = [1,1.2,1.6,2];
if useMax
    allMP = cat(2,resultList.maxProj);
else
    allMP = cat(2,resultList.meanProj);
end
% allMP = (allMP+allMP(end:-1:1,:))/2;
allMP = allMP./repmat(sum(allMP,1),size(allMP,1),1);

for ct = 1:length(stages)-1,
    figure('Name',sprintf('%s %1.1f -> %1.1f',overallName, stages(ct:ct+1)));
    sidx = find(spindleLength>stages(ct) & spindleLength<stages(ct+1));
    averageMP = mean(allMP(:,sidx),2);

    samp=sum(averageMP);
    averageMP = averageMP/samp;

    hold on
    for i=1:length(sidx)
        plot(xLabels(3:end-2)/24,allMP(3:end-2,sidx(i))/samp,'b');
    end
    plot(xLabels(3:end-2)/24,averageMP(3:end-2),'r','LineWidth',1.5)
end




% reset warnings
warning(warningState)



