function [resultList,plotData,plotStruct] = bilobesGfpKinetochores(project, normalize,doPlot,improveAlignment,redoAll,zOrientation,movieType)
%BILOBESPINDLES projects kinetochore signals onto the spindle axis in two-color images
%
% SYNOPSIS: bilobeSpindles
%
% INPUT project : type of projection: 'max'/{'sum'}
%       normalize : what to set to 1 {'sum'}/'max'/'none'
%           for kinetochores: project sum, norm sum
%           for tubulin: project sum, norm max
%       doPlot : whether to plot or not {1}/0
%       improveAlignment : whether or not to improve alignment. Default:
%                 true
%       redoAll : if 0, code will not perform detection if it has been done
%                  before. Default: false
%       zOrientation : orientation of the z-axis: {'perpendicular'} to the
%                      other two/'parallel' to image-z
%       moveType : {'raw'} or 'decon' (recommended for tub1)
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
align_1=true;
spbCorrection_2 = [-0.0029    0.3407    0.2784];
align_2=true;
%spbCorrection_2 = [-0.0480 0.3796 0.3398]; % January 2007

spbCorrection_3 = [-0.0417    0.1653    0.3145]; % February 2007 (2-7-7)
align_3=true;
spbCorrection_4 = [-0.0517    0.1898    0.3650]; % February 20 2007
align_4=false;

debug = false;
alignmentTolerance = 0.3;
maxIter = 10;

% turn off exposure time/nd warnings
warningState = warning;
warning('off','R3DREADHEADER:exposureTimeChanged');
warning('off','R3DREADHEADER:ndFilterChanged');

% PARAMETERS THAT CANNOT BE CHANGED FROM INPUT
plotGrid = false; % imaris will plot the 3D-grid over the aligned two-color data
plotInt = false; % imaris will plot spb and ndc80


% PARAMETERS THAT CAN BE CHANGED FROM INPUT
plotGallery = true; % plot gallery of projections
def_project = 'sum'; % maximum intensity projection
def_normalize = 'sum';
def_doPlot = true;
def_improveAlignment = true;
def_zOrientation = 'para'; % perp(endicular) / par(allel)
def_redoAll = false; % if false, code will not redo analysis
def_movieType = 'raw'; % or 'decon'

if nargin < 1 || isempty(project)
    project = def_project;
end
if nargin < 2 || isempty(normalize)
    normalize = def_normalize;
end
if nargin < 3 || isempty(doPlot)
    doPlot = def_doPlot;
end
if nargin < 4 || isempty(improveAlignment)
    improveAlignment = def_improveAlignment;
end
if nargin < 5 || isempty(redoAll)
    redoAll = def_redoAll;
end
if nargin < 6 || isempty(zOrientation)
    zOrientation = def_zOrientation;
end
if nargin < 7 || isempty(movieType)
    movieType = def_movieType;
end
if nargout == 3
    doPlot = true;
end

% search files
cdBiodata;
[fileList,tokens]=searchFiles('(\w+)_R3D\.dv|(\w+)\.r3d','log|DIC','ask');

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

if nData > 1 && (plotInt || plotGrid)
    selection = questdlg(sprintf(...
        'You selected to plot %i images in Imaris. Do you want to cancel?',...
        nData),'Warning','Continue','Cancel','Cancel');
    if strcmp(selection,'Cancel')
        error('cancelled by user')
    end
end

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
    % no problem if no decon if we don't use it
    if exist(nameList(iData).deconMovieName,'file') || strcmp(movieType,'raw')
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
            
            % correct a small mistake of Eugenio - the software said that
            % the lens was 10x
            if movieHeader.lensID == 10105 && movieHeader.pixelX > 0.6
                movieHeader.pixelX = movieHeader.pixelX/10;
                movieHeader.pixelY = movieHeader.pixelY/10;
            end

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
        elseif nSpots < 2
            status(iData) = -1;
        else
            idlist_L = idlist;
        end

        if ~isempty(idlist_L) && status(iData) > 0
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
            if isempty(improveAlignment)
                improveAlignment = align_1;
            end
        elseif datenum(d.date) < datenum('1-Feb-2007 00:00:01')
            spbCorrection = spbCorrection_2;
            if isempty(improveAlignment)
                improveAlignment = align_2;
            end
        elseif datenum(d.date) < datenum('20-Feb-2007 00:00:01')
            spbCorrection = spbCorrection_3;
            if isempty(improveAlignment)
                improveAlignment = align_3;
            end
        else
            spbCorrection = spbCorrection_4;
            if isempty(improveAlignment)
                improveAlignment = align_4;
            end
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
            inPlaneAngle = atan2(e_spb(2),e_spb(1))*180/pi;

            % make perpendicular vector parallel to xy
            e_perp = [-e_spb(2),e_spb(1),0];
            e_perp = e_perp / sqrt(sum(e_perp.^2));

            % complete coordinate system
            if strmatch(zOrientation,'perp')
                e_3 = cross(e_spb,e_perp);
            else
                e_3 = [0,0,1];
            end

            % convert spb coords in pixels
            load(nameList(iData).dataPropertiesName);
            pix2mu = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];
            %spbCoordNdc = (spbCoord + repmat(spbCorrection,2,1)) ./ repmat(pix2mu,2,1);
            spbCoord = spbCoord ./repmat(pix2mu,2,1);

            % spb-vector: 1/24th of spindle length. Perpendicular vectors: 1 pix
            e_spb = (e_spb * n_spb/24)./pix2mu;
            e_perp = (e_perp * dataProperties.PIXELSIZE_XY)./pix2mu;
            e_3 = (e_3 * dataProperties.PIXELSIZE_XY)./pix2mu;

            % load image - deconvolved movie
            switch movieType
                case 'decon'
                    rawMovie = r3dread(nameList(iData).deconMovieName);
                case 'raw'
                    rawMovie = r3dread(nameList(iData).rawMovieName);
            end

            % arbitrary grid: 15 pix perpendicular, 24 bins parallel to spindle
            % axis (Gardner)
            % -> use a few more : 21 perp, bins from -3:27
            xLabels = -2.5:26.5;

            % loop to find (somewhat) centered intensity distribution.
            % Idea: in axis direction, the two lobes should decrease to
            % "zero" symmetrically, in the perpendicular directions, the
            % maximum should be in the center
            ndcShift = spbCorrection./pix2mu;




            %shiftNorm = 99;
            %fh1=figure;
            %fh2=figure;
            %             fh3=figure;
            if debug
                fh4 = figure;
            end
            %             ah = gca; hold on;cm=jet(9);
            %             for i=1:9

            % improve alignment only on demand
            if improveAlignment
                ct = 1;
            else
                ct = 10; %run once only to get at intensities
            end
            done = false;
            store = zeros(maxIter,3);
            while  ~done % done if maxIter interations, shift<alignmentTolerance pix
                % read intensity
                spbCoordNdc = spbCoord + repmat(ndcShift,2,1);
                [xGrid,yGrid,zGrid] = arbitraryGrid(e_spb, e_perp, e_3, spbCoordNdc(1,:),...
                    [xLabels(1) xLabels(end)], [-10 10], [-10 10], size(rawMovie));
                % interpolate ndc80
                interpImg = interp3(rawMovie(:,:,:,2),yGrid,xGrid,zGrid,'*linear',NaN);

                % estimate background: Median of all 12 edges
                edgeMask = true(size(interpImg));
                edgeMask(2:end-1,2:end-1,2:end-1) = false;
                edgeMask(2:end-1,2:end-1,[1,end]) = false;
                edgeMask(2:end-1,[1,end],2:end-1) = false;
                edgeMask([1,end],2:end-1,2:end-1) = false;

                edgeList = interpImg(edgeMask);
                background = nanmedian(edgeList);

                % make maximum, average projections
                maxProj = nanmax(nanmax(interpImg,[],3),[],2)-background;
                meanProj = nanmean(nanmean(interpImg,3),2)-background;

                % make perpendicular projections
                p_perp = squeeze(nanmean(interpImg,1))-background;

                % plot - debug
                %                 figure(fh1)
                %                 subplot(3,3,i),imshow(p_perp,[]);
                %                 figure(fh2)
                %                 subplot(3,3,i),imshow(nanmean(interpImg,3),[]);
                %                 plot(ah,xLabels,(meanProj-min(meanProj))/max(meanProj-min(meanProj)),'Color',cm(i,:))
                if ct == 1 && debug
                    figure(fh4)
                    subplot(2,2,1),imshow(nanmean(interpImg,3),[]);
                    subplot(2,2,2),imshow(p_perp,[]);
                end


                % in axis direction: Take average of 3 end points, find
                % centroid - this is not stable
                % c_spb = sum([meanProj(1:3)'.*xLabels(1:3),...
                % meanProj(end-2:end)'.*xLabels(end-2:end)])/sum(meanProj([1:3,end-2:end]));
                % find by how much we need to shift along the axis
                %c_spb = (c_spb - 12)/2; % divide by 2 b/c of instability
                % instead, check where the curve goes below 10% of the maximum intensity;
                % count how many pixels below 10% are on either side
                pp=csaps(xLabels,(meanProj-min(meanProj))/max(meanProj-min(meanProj))-0.1);
                left=fnzeros(pp,xLabels([1,8]));
                right=fnzeros(pp,xLabels([end-7,end]));
                % there are 4 options
                % - below 10% on both sides
                % - below 10% on one side (2)
                % - never below 10% on either side
                % if below 10% on both sides, the center is right between
                % the zero-crossings
                % if below 10% on one side only, the center should be
                % shifted by as if crossing was right beyond border
                % if never below 10% at the borders, the center is at the
                % centroid of the data
                emptyLeft = isempty(left);
                emptyRight = isempty(right);
                switch emptyLeft + 2*emptyRight
                    case 0
                        % none empty
                        c_spb = (max(left(:)) + min(right(:)))/2 - 12;
                    case 1
                        % only right
                        c_spb = (xLabels(1)-1 + min(right(:)))/2 -12;
                    case 2
                        % only left
                        c_spb = (max(left(:)) + xLabels(end)+1)/2 -12;
                    case 3
                        % both empty
                        c_spb = sum(xLabels.*meanProj')/sum(meanProj) - 12;
                end


                % find centroid of perpendicular projection
                c_perp = centroid2D(p_perp);
                % find by how much we need to shift in the third dimension
                c_perp = c_perp([2,1]) - [11,11];
                store(ct,:) = [c_spb,c_perp];

                if ct == maxIter-1
                    % check for oscillation in the shift vector
                    signVec = sign(store(6:9,:));
                    for d=1:3
                        if all((signVec(:,d)+signVec([end,1:3],d))==0)
                            % update store. Difference/4 is half the amplitude
                            store(ct,d) = (store(ct,d)-store(ct-1,d))/4;
                        end
                    end

                end

                % transform the correction
                c_vec = [e_spb(:),e_perp(:),e_3(:)]*store(ct,:)';

                shiftNorm = norm([c_spb,c_perp]);

                % check if done (e.g. at ct=10, or very little shift)
                if shiftNorm < alignmentTolerance || ct >= maxIter;
                    done = true;
                else

                    % update ndc_shift
                    ndcShift = ndcShift + c_vec';

                    % debug
                    if debug
                        disp(sprintf('%1.3f -- %1.3f %1.3f %1.3f',shiftNorm,store(ct,:)))
                    end

                    % counter to make sure that we will exit the loop
                    % eventually
                    ct = ct + 1;
                end

            end

            % debug
            if debug
                figure(fh4)
                subplot(2,2,3),imshow(nanmean(interpImg,3),[]);
                subplot(2,2,4),imshow(p_perp,[]);
            end



            %     % show grid - (multiply by pix2mu for um)
            if plotGrid
                %                             d1 = spbCoordNdc(1,:) + 15 * e_spb;
                %                             d2 = spbCoordNdc(1,:) + 15 * e_perp;
                %                             d3 = spbCoordNdc(1,:) + 15 * e_3;
                %                 pd(3).XYZ = [d1;d2;d3];pd(3).color = [0,0,1,0];pd(3).spotRadius = 0.03;
                [xGrid,yGrid,zGrid] = arbitraryGrid(e_spb, e_perp, e_3, spbCoordNdc(1,:),...
                    [0.5,23.5], [-10 10], [-10 10], size(rawMovie));
                % make grid-box
                %                 xGrid(2:end-1,2:end-1,2:end-1) = NaN;
                %                 yGrid(2:end-1,2:end-1,2:end-1) = NaN;
                %                 zGrid(2:end-1,2:end-1,2:end-1) = NaN;

                % make three grid-planes
                xGrid(2:end,[1:10,12:end],[1:10,12:end]) = NaN;
                yGrid(2:end,[1:10,12:end],[1:10,12:end]) = NaN;
                zGrid(2:end,[1:10,12:end],[1:10,12:end]) = NaN;
                xGrid = xGrid(~isnan(xGrid));
                yGrid = yGrid(~isnan(yGrid));
                zGrid = zGrid(~isnan(zGrid));
                pd(2).XYZ = [xGrid(:)*pix2mu(1),yGrid(:)*pix2mu(2),zGrid(:)*pix2mu(3)];
                pd(2).color = [0.6403, 0.6513, 0.6464, 0];pd(2).spotRadius = 0.012;
                pd(1).XYZ = spbCoordNdc.*repmat(pix2mu,2,1);
                pd(1).color = [1,0,0,0];pd(1).spotRadius = 0.05;
                % shift the spb channel
                [xx,yy,zz] = ndgrid((1:size(rawMovie,1))-ndcShift(1),...
                    (1:size(rawMovie,2))-ndcShift(2),(1:size(rawMovie,3))-ndcShift(3));
                im = interp3(rawMovie(:,:,:,1),yy,xx,zz,'*linear',NaN);
                imarisPlot3(pd,pix2mu,cat(4,im,rawMovie(:,:,:,2)))
            end


            % read intensities from spb image
            [xGrid,yGrid,zGrid] = arbitraryGrid(e_spb, e_perp, e_3, spbCoord(1,:),...
                [xLabels(1) xLabels(end)], [-10 10], [-10 10], size(rawMovie));
            interpImgSpb = interp3(rawMovie(:,:,:,1),yGrid,xGrid,zGrid,'*linear',NaN);




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
            resultList(iData).inPlaneAngle = inPlaneAngle;
            resultList(iData).name = nameList(iData).deconMovieName;
            resultList(iData).ndcShiftPix = ndcShift;

            % show two plots
            if plotInt
                imarisShowArray(cat(5,interpImg,interpImgSpb));
            end
        end
    end
end

% remove bad results
resultList(status==-1)=[];

% get spindle length
spindleLength = cat(1,resultList.n_spb);
[dummy,sortIdx] = sort(spindleLength);
longestSpindle = max(spindleLength)*30/24;
% show gallery if selected
if plotGallery
    nData = length(resultList);
    nRows = floor(sqrt(nData));
    nCols = ceil(nData/nRows);
    proj = {'maximum projection. length--angle(z/xy)--idx','mean projection. length--angle(z/xy)--idx'};
    for p=1:2
        fh = figure('Name',proj{p});
        for d=sortIdx'
            % gfp: green, spb: red

            zeroSize= size(resultList(d).interpImg);
            zeroImage = zeros(zeroSize(1:2));
            switch p
                case 1 %max
                    red = nanmax(resultList(d).interpImgSpb,[],3);
                case 2 %mean
                    red = nanmean(resultList(d).interpImgSpb,3);
            end
            red = red - nanmin(red(:));
            red = red./nanmax(red(:));
            switch p
                case 1 %max
                    green = nanmax(resultList(d).interpImg,[],3);
                case 2 %mean
                    green = nanmean(resultList(d).interpImg,3);
            end
            green = green - nanmin(green(:));
            green = green./nanmax(green(:));
            image = cat(3,red,...
                green,...
                zeroImage);
            h=subplot(nRows,nCols,find(d==sortIdx));
            % turn the images to the right and center - it's easier on the eyes
            subimage((-14.5:14.5)*norm(resultList(d).e_spb.*pix2mu),...
                (-10:10)*dataProperties.PIXELSIZE_XY,permute(image,[2,1,3]));
            set(h,'XLim',[-longestSpindle/2,longestSpindle/2])

            %subimage turns visibility on - but we need it if we want to
            %show title. Therefore, make everything the same color as
            %figure
            set(h,'XTickLabel','','YTickLabel','','Color',get(fh,'Color'),'XColor',get(fh,'Color'),'YColor',get(fh,'Color'))
            movieName = resultList(d).name;
            fsidx = findstr(movieName,filesep);
            movieName = movieName(fsidx(end)+1:end); % later: replace _ with -, cut _R3D...
            movieName = regexprep(movieName,'_','-');
            endIdx = findstr(lower(movieName),'-r3');
            movieName = movieName(1:endIdx-1);
            title(h,{sprintf('%1.2f--%2.0f--%i',...
                resultList(d).n_spb, ...
                resultList(d).spbVectorAngle,...  %resultList(d).inPlaneAngle,...
                d);...
                sprintf('%s',movieName)});
        end
    end
end




switch project
    case 'max'
        projectedIntensities = cat(2,resultList.maxProj);
    case 'sum'
        projectedIntensities = cat(2,resultList.meanProj);
end
% normalize projectedIntensities so that every movie has the same total
% intensity
switch normalize
    case 'sum'
        projectedIntensities = projectedIntensities./...
            repmat(sum(projectedIntensities,1),size(projectedIntensities,1),1);
    case 'max'
        projectedIntensities = projectedIntensities./...
            repmat(max(projectedIntensities,[],1),size(projectedIntensities,1),1);
    case 'none'
        % normalize only so that the overall max is one
        projectedIntensities = projectedIntensities/max(projectedIntensities(:));
end

% plot data
plotData = {spindleLength,projectedIntensities(3:end-2,:)};
if doPlot
    plotStruct = bilobePlot(plotData);
end


% % plotData: data needed for plotting
% plotData(1:3) = struct('xall',[],'yall',[],'zall',[],'yTickLabels',[]);
%
% % plot averages
% figureNames = {'asymmetric','a, max=1','symmetric','s, max=1'};
% for i=1:4
%
%     boundaries = [0.9:0.025:1.9;1.1:0.025:2.1];
%     meanBoundaries = mean(boundaries,1);
%     nBoundaries = size(boundaries,2);
%     nBins = length(xLabels)-4; %make compatible with bilobeDistribution
%     xall = zeros(nBins,nBoundaries);
%     yall=zeros(nBins,nBoundaries);
%     zall=zeros(nBins,nBoundaries);
%     zallS = zeros(nBins,nBoundaries);
%     yTickLabels = cell(nBoundaries,1);
%     nSpindles = zeros(nBoundaries,1);
%     spindleLength = cat(1,resultList.n_spb);
%     spindleIdx = cell(nBoundaries,1);
%
%     switch i
%         case {1,2}
%             if useMax
%                 allMP = cat(2,resultList.maxProj);
%             else
%                 allMP = cat(2,resultList.meanProj);
%             end
%         otherwise
%             if useMax
%                 allMP = cat(2,resultList.maxProj);
%             else
%                 allMP = cat(2,resultList.meanProj);
%             end
%             allMP = (allMP+allMP(end:-1:1,:))/2;
%     end
%     % normalize data.
%     switch i
%         case {4,2}
%             allMP = allMP./repmat(max(allMP,[],1),size(allMP,1),1);
%         case {1,3}
%             allMP = allMP./repmat(sum(allMP,1),size(allMP,1),1);
%     end
%
%     % figure, hold on
%     for ct = 1:nBoundaries,
%         spindleIdx{ct} = (spindleLength>boundaries(1,ct) & spindleLength<boundaries(2,ct));
%         if any(spindleIdx{ct})
%             % make weighted average - points 0.1um away have no weight
%             weights = (0.1-abs(spindleLength(spindleIdx{ct})-meanBoundaries(ct)))/0.1;
%             [averageMP,dummy,sigmaMP] = weightedStats(allMP(:,spindleIdx{ct})',...
%                 weights,'w');
%             %old: averageMP = mean(allMP(:,spindleIdx{ct}),2);
%             switch i
%                 case {4,2}
%                     sigmaMP = sigmaMP/max(averageMP);
%                     averageMP = averageMP/max(averageMP);
%
%                 otherwise
%                     sigmaMP = sigmaMP/sum(averageMP);
%                     averageMP = averageMP/sum(averageMP);
%
%             end
%             %plot3(x,ct*ones(size(x)),z),
%             nSpindles(ct) = sum(weights);
%             %old: nSpindles(ct) = nnz(spindleIdx{ct});
%
%             zall(:,ct)=averageMP(3:end-2);
%             zallS(:,ct) = sigmaMP(3:end-2);
%         else
%             % there are no spindles this long
%             zall(:,ct) = NaN;
%             zallS(:,ct) = NaN;
%         end
%         xall(:,ct)=xLabels(3:end-2)/24;
%         yall(:,ct)=meanBoundaries(ct);
%         yTickLabels{ct}=sprintf('%1.1f/%1.2f', ...
%             meanBoundaries(ct),nSpindles(ct));
%
%     end
%
%     % check zallS for inf
%     zallS(~isfinite(zallS)) = NaN;
%
%     figure('Name',[overallName,' ',figureNames{i}])
%     if plotSigma
%         ah = subplot(1,2,1);
%     else
%         ah = gca;
%     end
%     contourf(xall,yall,zall,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
%     %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
%     %     axis tight
%     set(ah,'yTick',1:0.1:2,'yTickLabel',yTickLabels(1:4:end),'yGrid','on')
%
%     switch i
%         case {4,2}
%             set(ah,'CLim',[0,1])
%         otherwise
%             set(ah,'CLim',[0,nanmax(zall(:))])
%     end
%     if ~plotSigma
%         colorbar('peer',ah)
%     end
%     % add lines
%     hold on
%     for d=0.2:0.2:0.8
%         line(d./meanBoundaries(meanBoundaries>=2*d),...
%             meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
%         line(1-d./meanBoundaries(meanBoundaries>=2*d),...
%             meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
%     end
%
%     if plotSigma
%         ah = subplot(1,2,2);
%         contourf(xall,yall,zallS,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
%         %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
%         %     axis tight
%         set(ah,'yTick',1:0.1:2,'yTickLabel',yTickLabels(1:4:end),'yGrid','on')
%
%         switch i
%             case {4,2}
%                 set(ah,'CLim',[0,1])
%             otherwise
%                 set(ah,'CLim',[0,nanmax(zall(:))])
%         end
%         colorbar('peer',ah)
%         % add lines
%         hold on
%         for d=0.2:0.2:0.8
%             line(d./meanBoundaries(meanBoundaries>=2*d),...
%                 meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
%             line(1-d./meanBoundaries(meanBoundaries>=2*d),...
%                 meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
%         end
%     end
%
%
%     %     xlim([0,1])
%     %     ylim([1,2])
%     %view([0 90])
%     plotData(i).xall = xall;
%     plotData(i).yall = yall;
%     plotData(i).zall = zall;
%     plotData(i).zallS = zallS;
%     plotData(i).yTickLabels = yTickLabels;
%     plotData(i).spindleIdx = spindleIdx;
%
%
% end
%
% figure('Name',sprintf('%s individual; sum=1',overallName));
% [sortedSpindleLength,sortIdx]=sort(spindleLength);
% int = cat(2,resultList.maxProj);
% int = int(:,sortIdx)';
% int = int./repmat(sum(int,2),1,size(int,2));
% intSymm = 0.5*(int + int(:,end:-1:1));
% subplot(1,2,1),imshow([int,NaN(size(int,1),1),intSymm],[]),
% colormap jet,
% subplot(1,2,2),
% plot(sortedSpindleLength,length(sortIdx):-1:1,'-+')
%
% % loop to make histograms
% stages = [1,1.2,1.6,2];
% if useMax
%     allMP = cat(2,resultList.maxProj);
% else
%     allMP = cat(2,resultList.meanProj);
% end
% % allMP = (allMP+allMP(end:-1:1,:))/2;
% allMP = allMP./repmat(sum(allMP,1),size(allMP,1),1);
%
% for ct = 1:length(stages)-1,
%     figure('Name',sprintf('%s %1.1f -> %1.1f',overallName, stages(ct:ct+1)));
%     sidx = find(spindleLength>stages(ct) & spindleLength<stages(ct+1));
%     averageMP = mean(allMP(:,sidx),2);
%
%     samp=sum(averageMP);
%     averageMP = averageMP/samp;
%
%     hold on
%     for i=1:length(sidx)
%         plot(xLabels(3:end-2)/24,allMP(3:end-2,sidx(i))/samp,'b');
%     end
%     plot(xLabels(3:end-2)/24,averageMP(3:end-2),'r','LineWidth',1.5)
% end




% reset warnings
warning(warningState)



