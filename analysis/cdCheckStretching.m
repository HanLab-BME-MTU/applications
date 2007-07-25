function cdCheckStretching(type)
%CDCHECKSTRETCHING attempts to quantify fluorescent tag stretching
%
% SYNOPSIS: cdCheckStretching
%
% INPUT type: 1 - movie
%             2 - bigImage
%
% OUTPUT -- not defined yet
%
% REMARKS
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: jdorn
% DATE: 26-Feb-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warningState = warning;
warning off MATLAB:divideByZero

% make list of idlists - don't use cen1*-movies
idlistList = loadIdlistList([],'~any(strcmp(idlist(1).stats.labelcolor,''cen1*''))');
nIdlists = length(idlistList);

% idlistList contains the latest idlists. If it's an idlisttrack we should
% adjust intensities. If it is an idlist, we should adjust both position
% and intensities (b/c tags could have been removed)

% if it's a movie, fitting is straightforward. If it's a bigImage, we need
% to identify individual cells first before fitting.

switch type

    case 1

        % movie

        for iIdlist = 1:nIdlists

            % inputStructure has idlist, dataProperties, rawMovieName, movieDir
            % therefore: find the raw movie
            [dummy,dummy,loadStruct] = cdLoadMovie('corr/raw',...
                idlistList(iIdlist).dirName,-1);

            % make slist from idlist (!! this is a subfunction, not the
            % outdated function with the same name!!)
            out = idlist2slist(idlistList(iIdlist).idlist,...
                idlistList(iIdlist).slist,...
                idlistList(iIdlist).dataProperties);

            fitStruct = struct('slist',out,...
                'dataProperties',idlistList(iIdlist).dataProperties,...
                'movieDir',idlistList(iIdlist).dirName,...
                'rawMovieName',loadStruct.movieName);
            
            % for isolated spots, position should always be ok (and, in
            % principle, amplitude, though it will depend on sigma). To get
            % correct sigma, fit width of isolated spots first
            sigmaCorrectionOld = fitStruct.dataProperties.sigmaCorrection;
            originalSigma = fitStruct.dataProperties.FT_SIGMA([1,3])./sigmaCorrectionOld;
            spotIdx = catStruct(1,'out.sp.idxL');
            [dummy,spbIdx] = ismember({'spb1';'spb2'},idlistList(iIdlist).idlist(1).stats.labelcolor);
            spbSpotIdxL = ismember(spotIdx,spbIdx);
            done = false;
            while ~done
            gaussFit = cdFitGauss(fitStruct,{'a','sxy','s3','b'}); %#ok<NASGU>
            % only care about spb
             
             s = catStruct(1,'gaussFit.coords(:,5:6)'); 
             s=s(spbSpotIdxL,:);
            sigmaCorrection = robustMean(s(~isnan(s(:,1)),:))./originalSigma;
            disp(sigmaCorrection);
            fitStruct.dataProperties.sigmaCorrection = ...
                sigmaCorrection;
            fitStruct.dataProperties = defaultDataProperties(...
                fitStruct.dataProperties);
            if all((abs(sigmaCorrection-sigmaCorrectionOld)./sigmaCorrectionOld)<0.05)
                done = true;
            else
                sigmaCorrectionOld = sigmaCorrection;
            end
            end % while
            
            
            % check type. If we come from idlist_L, we want to adjust
            % positions, too
            if findstr(idlistList(iIdlist).idlist(1).stats.idname,'track_L')
                % only fit amplitude
                gaussFit = cdFitGauss(fitStruct,{'a','b'});
            else
                % if it was idlisttrack: find idlist_L
                if findstr(idlistList(iIdlist).idlist(1).stats.idname,'track')
                    if strcmp(idlist(iIdlist).type(end),'2')
                        load(idlist(iIdlist).dataFileName,'idlist_L2')
                        idlistList(iIdlist).idlist = idlist_L2;
                        fitStruct.idlist = idlist_L2;
                    else
                        load(idlist(iIdlist).dataFileName,'idlist_L')
                        idlistList(iIdlist).idlist = idlist_L;
                        fitStruct.idlist = idlist_L;
                    end
                end
                % fit amplitude and position
                gaussFit = cdFitGauss(fitStruct,{'x1','x2','x3','a','b'});
            end

            pix2mu = [fitStruct.dataProperties.PIXELSIZE_XY,...
                fitStruct.dataProperties.PIXELSIZE_XY,...
                fitStruct.dataProperties.PIXELSIZE_Z];

            % write gaussFit back into idlist - transform!
            for t=1:length(idlistList(iIdlist).idlist)
                % positions - if not fitted, they should have stayed the same
                if ~isempty(idlistList(iIdlist).idlist(t).linklist)
                    spotIdx = idlistList(iIdlist).idlist(t).linklist(:,2);
                    idxL = spotIdx>0; % good spotIdx
                    n = sum(idxL);
                    idlistList(iIdlist).idlist(t).linklist(idxL,9:11) = ...
                        gaussFit(t).coords(:,1:3).*repmat(pix2mu,n,1);
                    % amplitudes
                    idlistList(iIdlist).idlist(t).linklist(idxL,8) = gaussFit(t).coords(:,4);
                    % chi2
                    idlistList(iIdlist).idlist(t).linklist(idxL,12) = gaussFit(t).chi2(:);
                    % add background
                    idlistList(iIdlist).idlist(t).linklist(idxL,16) = gaussFit(t).coords(:,end);
                    % Q
                    if findstr(idlistList(iIdlist).idlist(1).stats.idname,'track')
                        % we only fitted amplitudes
                    else
                        % we should also update coordinate uncertainty
                    end
                end
            end

            % read intensities. We need the info in the fitStructure.
            % Therefore, update
            fitStruct.idlist = idlistList(iIdlist).idlist;
            intensities = cdCheckStretching_readIntensities(fitStruct);

            nTimepoints = length(intensities);
            spotIntensities=cat(2,intensities.spotIntensities)';
            nSpotsMax = length(intensities(1).spotIntensities);
            tt = repmat((1:nTimepoints)',1,nSpotsMax);
            [c{1:nSpotsMax}] = deal(ones(nTimepoints,1)); %#ok<AGROW>
            A = [blkdiag(c{:}),tt(:)];
            Y = spotIntensities(:);
            goodRows = isfinite(Y) & Y>0;

            % fit intensities and plot
            [xFit, stdX, goodRows, intAx] = robustExponentialFit2(Y(goodRows),A(goodRows,:),1);

            % correct spotIntensities for bleaching
            spotIntensitiesC = spotIntensities./repmat(exp(xFit(end)*(1:nTimepoints)'),1,nSpotsMax);
            figure('Name',sprintf('Bleach-corr int for %s',fitStruct.dataProperties.name))
            if nSpotsMax == 4
                colorOrder = [255 0 0;255 140 0;0 0 255;0 191 255]/255;
                legendCell = {'spb1' 'cen1' 'spb2' 'cen2'};
            else
                colorOrder = [255 0 0;50 205 50;0 0 255]/255;
                legendCell = {'spb1' 'cen' 'spb2'};
            end
            % plot corrected spotIntensities
            ah = subplot(2,1,1);
            set(ah,'ColorOrder',colorOrder,'NextPlot','add'),
            ph=plot(spotIntensitiesC./repmat(nanmedian(reshape(spotIntensitiesC(:,[1,3]),[],1)),nTimepoints,nSpotsMax));
            ylim([0,1.5])
            for i=1:length(ph)
                set(ph(i),'DisplayName',legendCell{i});
            end
            ah = subplot(2,1,2);
            %ah = axes;
            set(ah,'ColorOrder',colorOrder,'NextPlot','add'),
            ph=plot(spotIntensitiesC./repmat((nanmean(spotIntensitiesC(:,[1,3]),2)),1,nSpotsMax));
            ylim([0,1.5])
            for i=1:length(ph)
                set(ph(i),'DisplayName',legendCell{i});
            end
            title('above: normed by med spb int. Below: each frame normed by avg spb int')
            legend toggle

            [n_spindleVector, cenPosNorm,goodTimes] = ...
                cdProjectPositions(fitStruct.idlist);

            % cpn is cenPosNorm*spindleLength
            cpn = [zeros(size(cenPosNorm,1),1),cenPosNorm,ones(size(cenPosNorm,1),1)];% -0.5;
            cpn = cpn .* repmat(n_spindleVector,1,4);
            if nSpotsMax == 3
                cpn(:,3) = [];
            else
                % cpn is now s1 c1 s2 c2
                cpn = cpn(:,[1,2,4,3]);
            end

            tmp = nan(size(cpn));
            tmp(goodTimes,:)=cpn;
            cpn = tmp;

            % make two kymographs. One for spotIntensities alone, one with
            % spotIntensities plus in-between intensities. Only do the
            % in-between intensities if there is data between all spots.
            % For the spotIntensity kymograph, we just use detected spot
            % intensities. For the rest, we have to calculate projections

            % for spots: make psf-size
            [gaussXY] = calcFilterParms(...
                fitStruct.dataProperties.WVL,...
                fitStruct.dataProperties.NA,1.51,'gauss',...
                fitStruct.dataProperties.sigmaCorrection);

            % make 10nm bins
            xMatrix = repmat(min(cpn(:,1))-0.1:0.01:max(cpn(:,3))+0.1,nTimepoints,1);
            yMatrix = repmat((1:nTimepoints)',1,size(xMatrix,2));
            zMatrix = zeros(size(xMatrix));
            for t=goodTimes'
                for i=1:nSpotsMax
                    zMatrix(t,:) = zMatrix(t,:) + ...
                        exp(-(xMatrix(t,:)-cpn(t,i)).^2/(2*gaussXY^2)) *...
                        spotIntensitiesC(t,i)/mean(spotIntensitiesC(t,[1,3]));
                end
            end
            % remove rows that have no data
            zMatrix(zMatrix == 0) = NaN;
            % correct for possible fitting errors
            zMatrix(zMatrix < 0) = 0;

            %             figure,contourf(xMatrix,yMatrix,zMatrix,...
            %                 'LineStyle','none','LevelList',linspace(0,1,100));
            %             % green colormap
            %             cm=hsl2rgb([linspace(0.3,0.4,64)',linspace(0.7,0.8,64)',linspace(0.2,1,64)']);
            %             colormap(cm);
            % make black axes background for NaNs
            %             set(gca,'Color','k')

            % second kymograph. This one will be an RGB figure
            % green: intensity of spots
            [red,green,blue,white] = deal(zeros(size(zMatrix)));

            fieldList = {'s1c1' 'c1c2' 's2c2'};
            if nSpotsMax == 3;
                fieldList = fieldList([1,3]);
                cpnList = [1,2];
            else
                cpnList = [1,2,4];
            end

            for t=goodTimes'
                % check whether all data is there. If not, skip
                allDataThere=  ~isempty(intensities(t).s1c1Int) && ...
                    ~isempty(intensities(t).s2c2Int) &&...
                    ~isempty(intensities(t).c1c2Int);

                redList = [];
                greenList = [];
                whiteList = [];
                xList = [];

                if allDataThere(1:end-(4-nSpotsMax))

                    % read intensities and the corresponding values along
                    % the axis
                    cpnCt = 1;
                    for f = fieldList
                        % intensities are nanMedians (hopefully, this will be
                        % robust-y to too hot pixels
                        redInt = nanmedian(nanmedian(intensities(t).([f{1},'Int']),2),3);
                        redIntBg = nanmedian(intensities(t).([f{1},'Bg']),2);
                        redInt = redInt - redIntBg; % this takes care of uneven backgrounds

                        greenInt = nanmean(nanmean(intensities(t).([f{1},'Sp']),2),3);

                        whiteInt = nanmedian(nanmedian(intensities(t).([f{1},'IntRaw']),2),3);
                        whiteIntBg = nanmedian(intensities(t).([f{1},'BgRaw']),2);
                        whiteInt = whiteInt - whiteIntBg;

                        % x-values: take norm of first vector
                        xVec = intensities(t).([f{1},'Vec'])(1,:);
                        [xStep,xVecN] = normList(xVec);
                        cosGamma = dot(xVecN,intensities(t).s1s2Vec);
                        % projected step: step*cos(g). zero is at -0.5 pix
                        xValues = ((1:length(redInt))'-0.5)*xStep*cosGamma;

                        % add to correct value in cpn and store
                        xValues = xValues + cpn(t,cpnList(cpnCt));

                        redList = [redList;redInt];
                        greenList = [greenList;greenInt];
                        whiteList = [whiteList;whiteInt];
                        xList = [xList;xValues];

                        % up counter
                        cpnCt = cpnCt + 1;

                    end

                    % add to red, green
                    red(t,:) = interp1(xList,redList,xMatrix(t,:),'linear',NaN);
                    green(t,:) = interp1(xList,greenList,xMatrix(t,:),'linear',NaN);
                    white(t,:) = interp1(xList,whiteList,xMatrix(t,:),'linear',NaN);

                end
            end




            rgb=cat(3,red,green,blue);
            maxRgb = max(rgb(:));
            rgb = rgb/maxRgb;
            rgb = rgb./repmat(exp(xFit(end)*(1:nTimepoints)'),[1,size(red,2),3]);
            % reorder for proper display
            rgb = rgb(end:-1:1,:,:);

            % show image
            %             figure,imshow(rgb);
            g=green/maxRgb./repmat(exp(xFit(end)*(1:nTimepoints)'),[1,size(red,2)]);
            g(g==0) = NaN;
            g(g<0) = 0;
            figure('Name',...
                sprintf('%s: theoretical int.',idlistList(iIdlist).dataProperties.name)),
            contourf(xMatrix,yMatrix,g,...
                'LineStyle','none','LevelList',linspace(0,1,100));
            cmg=isomorphicColormap('green');
            colormap(cmg)
            set(gca,'Color','k')
            
            w=white/max(white(:))./repmat(exp(xFit(end)*(1:nTimepoints)'),[1,size(red,2)]);
            w(w==0) = NaN;
            w(w<0) = 0;
            figure('Name',sprintf('%s: raw int.',idlistList(iIdlist).dataProperties.name)),
            contourf(xMatrix,yMatrix,w,...
                'LineStyle','none','LevelList',linspace(0,1,100));
            cmw=isomorphicColormap('blue');
            colormap(cmw)
            set(gca,'Color','k')

            r=red/maxRgb./repmat(exp(xFit(end)*(1:nTimepoints)'),[1,size(red,2)]);
            r(r==0) = NaN;
            r(r<0) = 0;
            figure('Name',sprintf('%s: residual int.',idlistList(iIdlist).dataProperties.name)),
            contourf(xMatrix,yMatrix,r,...
                'LineStyle','none','LevelList',linspace(0,1,100));
            cmr=isomorphicColormap('red');
            colormap(cmr)
            set(gca,'Color','k')

           

            disp('stopping. Press key to continue')
            keyboard

            % figure,hold on
            % for i=1:nSpotsMax,
            %     stem3(cpn(:,i),goodTimes,spotIntensitiesC(goodTimes,i),'-','Color',extendedColors(i))
            %     plot3(cpn(:,i),goodTimes,spotIntensitiesC(goodTimes,i),'-','Color',extendedColors(i))
            % end
            % box on
            % grid on
            % zlim([0,nanmax(spotIntensitiesC(:))])

            % xx = cpn(:,repeatEntries(1:nSpotsMax,3));
            % zz = spotIntensitiesC(:,repeatEntries(1:nSpotsMax,3));
            % yy = repmat((1:nTimepoints)',1,3*nSpotsMax);
            %
            % zz(:,3:3:end) = NaN;
            % xx(:,1:3:end) = xx(:,1:3:end) - 0.04;
            % xx(:,2:3:end) = xx(:,2:3:end) + 0.04;
            % xx(:,3:3:end) = xx(:,3:3:end) + 0.06;
            %
            % figure,surf(xx,yy,zz,'FaceColor','interp','FaceLighting','phong','LineStyle','none')
            %  cm=hsl2rgb([linspace(0.6,0.7,64)',linspace(0.7,0.8,64)',linspace(0.2,1,64)']);
            %     colormap(cm)

            %             % for spots: make psf-size
            %             [gaussXY] = calcFilterParms(...
            %                 fitStruct.dataProperties.WVL,...
            %                 fitStruct.dataProperties.NA,1.51,'gauss');
            %
            %             % make 10nm bins
            %             xMatrix = repmat(min(cpn(:,1))-0.1:0.01:max(cpn(:,3))+0.1,nTimepoints,1);
            %             yMatrix = repmat((1:nTimepoints)',1,size(xMatrix,2));
            %             zMatrix = zeros(size(xMatrix));
            %             for t=goodTimes'
            %                 for i=1:nSpotsMax
            %                     zMatrix(t,:) = zMatrix(t,:) + ...
            %                         exp(-(xMatrix(t,:) - cpn(t,i)).^2/(2*gaussXY^2)) * spotIntensitiesC(t,i);
            %                 end
            %             end
            %             % this is pretty, but we will have to make a kymograph instead
            %             figure,surf(xMatrix,yMatrix,zMatrix,'FaceColor','interp','FaceLighting','phong','LineStyle','none')
            %             cm=hsl2rgb([linspace(0.6,0.7,64)',linspace(0.7,0.8,64)',linspace(0.2,1,64)']);
            %             colormap(cm)
            %             zlim([0,nanmax(spotIntensitiesC(:))])
            %             view([1.5 50])
            %             light('Position',[-1,-5,0],'Style','local')
            %             light('Position',[3,-5,0],'Style','local')
            %             light('Position',[1,100,3],'Style','local')

        end % loop idlists

    case 2
        % bigImage

        % do cell-wise fitting. Potentially only load part of the image
        for iIdlist = 1:nIdlists

            % calculate distances
            coords = idlistList(iIdlist).idlist(1).linklist(:,9:11);
            distList = pdist(coords); % default is euclidean

            % cluster with complete linkage
            links = linkage(distList,'complete');
            % make clusters if the longest distance between two points is less than
            % 2.2 microns
            labels = cluster(links,'cutoff',2.2,'criterion','distance');

            % loop through individual cells and fit
            nCells = max(labels);

            for c = 1:nCells

                cIdx = labels == c;
                nTags = length(cIdx);

                % only work with 3 or more tags
                if nTags > 2 && nTags < 5
                    % reduce the linklist in the idlist
                    idlistCell = idlistList(iIdlist).idlist;
                    idlistCell.linklist = idlistCell.linklist(cIdx);
                    idlistCell.stats.labelcolor = idlistCell.stats.labelcolor(cIdx);

                    slistCell = idlistList(iIdlist).slist;
                    slistCell.sp = slistCell.sp(cIdx);


                    % inputStructure has idlist, dataProperties, rawMovieName, movieDir
                    % therefore: find the raw movie
                    [dummy,dummy,loadStruct] = cdLoadMovie('corr/raw',...
                        idlistList(iIdlist).dirName,struct('noMovie',true));
                    fitStruct = struct('slist',slistCell,...
                        'dataProperties',idlistList(iIdlist).dataProperties,...
                        'movieDir',idlistList(iIdlist).dirName,...
                        'rawMovieName',loadStruct.movieName);

                    % fit amplitude and position
                    gaussFit = cdFitGauss(fitStruct,{'x','y','z','a','b'});

                    % positions - if not fitted, they should have stayed the same
                    idlistList(iIdlist).idlist.linklist(cIdx,9:11) = gaussFit.coords(:,1:3);
                    % amplitudes
                    idlistList(iIdlist).idlist.linklist(cIdx,8) = gaussFit.coords(:,4);
                    % chi2
                    idlistList(iIdlist).idlist.linklist(cIdx,12) = gaussFit.chi2(:);
                    % Q

                    % -- TBD (if we need it)
                end % if metaphase cell

                % need to cut out local intensity from bigImg!!

            end % loop fit cells

        end % loop idlists
end % switch bigImg vs movie



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ INTENSITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning(warningState)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% idlist2slist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = idlist2slist(idlist,slist,dataProperties)

% make a new slist from th eidlist
pix2mu = [dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

out(1:length(slist)) = struct('sp',[]);

for t=1:length(slist)
    if ~isempty(idlist(t).linklist)
        % keep the order of spots of the idlist; don't worry about 'true'
        % spot numbers
        spotIdx = find(idlist(t).linklist(:,2)>0);
        bgList = [slist(t).sp.bg];
        % tracker may assign new spot numbers. Therefore initialize with a
        % common bg
        background = mean(bgList);
        for i=1:length(spotIdx)
            out(t).sp(i).cord = idlist(t).linklist(spotIdx(i),9:11)./pix2mu;
            out(t).sp(i).amp = idlist(t).linklist(spotIdx(i),8);
            out(t).sp(i).bg = background;
            out(t).sp(i).idxL = spotIdx(i); % row of idlist corresponding to the spot
        end
    end
end