function idlist = linker(slist,dataProperties,verbose,DEBUG,maxSpots,movMax)
%LINKER is a function of the chromdyn project to link tags between frames
%
%
% OUTPUT idlist : 1-by-nTimepoints structure with fields:
%          .linklist 1: time, 2: spot#, 3: spotFlag, 4: tagIdx, 5: tagFlag,
%              6: linkup, 7: linkdown, 8: amp, 9-11: xyz, 12: sigma0
%                   spotFlags: 0: normal
%                              1: estimated spot
%                              2: found by tracker
%                              3: primary fusion spot
%                              4: secondary fusion spot
%                              5: primary fusion spot found by tracker
%
%                   tagFlags : 0: normal
%                              1: intensity has been corrected
%                              2: tag has only one appearance (this one)
%                              3: tag has only one appearance (not this
%                                 one)
%          .info
%          (1).stats
%
%
% recalc    -1 start from slist
%            0 full recalc, start from idlist
%            1 don't go through first linking
%            2 don't go thorugh source-target linking
%           r+10: also estimate fusions

%============
%% TEST INPUT
%============

if nargin < 2
    error('Not enough input arguments!')
end
if nargin < 3 || isempty(verbose)
    verbose = 1;
end
if nargin < 4 || isempty(DEBUG)
    DEBUG = 0;
end
if nargin < 5 || isempty(maxSpots)
    % do nothing
else
    dataProperties.MAXSPOTS = maxSpots;
end

% test for new slist
if ~isfield(slist,'COM') && ~isfield(slist,'linklist')
    error('please run latest version of spotDetector first!')
end

% assign defaults, read dataProperties
defaults = {...
    'MAXSPOTS','MAXSPOTS',4;...
    'numSources','linker_numSources',4;...
    'relativeMaxDistance','linker_relativeMaxDistance',-1;...
    'absoluteMaxDistance','linker_absoluteMaxDistance',-1;...
    'relAmpWeight','linker_relAmpWeight',1/1.5;...
    'diffusionFactor','linker_diffusionFactor',3;...
    'replaceAmplitudes','linker_replaceAmplitudes',0;...
    'useCOM','linker_useCOM',1;...
    'fuseRatio','linker_fuseRatio',1.5;...
    };
nDefaults = size(defaults,1);
% assign defaults
constants = cell2struct(defaults(:,3),defaults(:,1),1);
% overwrite defaults with input from dataProperties.
% All linker_ fields are optional
for iDefault = 1:nDefaults
    if isfield(dataProperties,defaults{iDefault,2})
        constants.(defaults{iDefault,1}) = ...
            dataProperties.(defaults{iDefault,2});
    end
end



% backwards compatibility
if isfield(dataProperties,'linkerProperties') && ~isempty(dataProperties.linkerProperties)
    fnames = fieldnames(dataProperties.linkerProperties);
    for name = fnames'
        constants.(name{1}) = dataProperties.linkerProperties.(name{1});
    end
    % backwards compatibility
    if isfield(dataProperties.linkerProperties,'maxDistance')
        constants.relativeMaxDistance = ...
            dataProperties.linkerProperties.maxDistance;
    end
end


% constant for how many sigmas will be used by the tracker as searchRadius
dataProperties.trackerRadiusMultiplicator = 5;
constants.lostChi2 = 1/dataProperties.trackerRadiusMultiplicator^2;

% pixel to microns
constants.pix2mu = [dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];

% remember name
constants.name = dataProperties.name;

% check if recalc
if isfield(slist,'stats') &&...
        isfield(slist(1).stats,'recalc') &&...
        ~isempty(slist(1).stats.recalc)
    recalc = slist(1).stats.recalc{1};

    % use defaults if not otherwise supplied
    if length(slist(1).stats.recalc) > 2 && ...
            ~isempty(slist(1).stats.recalc{3})
        constants.relAmpWeight = 1/slist(1).stats.recalc{3};
    end
    if length(slist(1).stats.recalc) > 1 && ...
            ~isempty(slist(1).stats.recalc{2})
        constants.MAXSPOTS = slist(1).stats.recalc{2};
    end
else
    recalc = -1;
end

if recalc > 5
    recalc = recalc - 10;
    doFusions = 1;
else
    doFusions = 0;
end

%============

%=======================
%% READ DATA FROM SLIST
%=======================

% data is read into idlist(t).linklist-structure
% [1: time, 2: spot#, 3: spotFlag, 4: tagIdx, 5: tagFlag,
%  6: linkup, 7: linkdown, 8: amp, 9-11: xyz, 12: sigma0 (detector)]
% assigns idlist.centroid (center of mass - weights=1), nSpots, ampList
% (sum of amplitudes). As temporary fields: distMatAmp, distMatXyz,
% sourceIdxList
if recalc == -1
    [idlist, nTimepoints, nSpots, ampList, goodTimes, goodIdx] = ...
        linkerReadSlist(slist, constants);
else
    [idlist, nTimepoints, nSpots, ampList, goodTimes, goodIdx, ...
        goodTimesM, t1t2, tagIndices, maxTagIdx, intList] = ...
        linkerReadIdlist(slist, constants);
    verbose = 0;
    intAx = [];
end
%=======================


%=======================
%% EXIT IF NO LINKING
%=======================

% for the moment, only break at 1 tp - allow setting to be changed via
% constants later
if nTimepoints < 2
    % add help, start history, expand Q-matrices
    maxTagIdx = nSpots; % every spot is a tag
    idlist(goodTimes).linklist(:,4) = idlist(goodTimes).linklist(:,2);
    % no debug possibility here
    idlist = finishIdlist(idlist,maxTagIdx,dataProperties,goodTimes,recalc,0);
    return
end
%=======================


%=======================
%% INITIALIZE LINKING
%=======================

% here, we
% (1) find bleaching (we know drift already)
% (2) build distMat of corrected intensities and corrected distances
% (3) calculate sigmaInt, sigmaPos from corrected data
%  Furthermore, we get goodTimesM (source-timepoints), t1t2 (list of
%  pairs of consecutive sources), and potentially intAx, the handle to the
%  intensity axes
if recalc < 1
    [idlist, goodTimesM, t1t2, intAx] = ...
        linkerInitializeLap(idlist, nSpots, ampList, goodIdx,...
        nTimepoints, constants, verbose);
end

%=======================


%=======================
%% SOLVE LAP
%=======================

% To ensure continuity, we first link all frames with MAXSPOTS tags (or the
% maximum number of tags). In a second step, all frames with fewer tags are
% linked from several frames with MAXSPOTS tags.

% tagIndices are the indices of all the tags in the frame, intList is the
% corresponding list of intensities, maxTagIdx the maximum tag index
if recalc < 1
    [idlist, tagIndices, intList, maxTagIdx] = ...
        linkerFirstLAP(idlist, nTimepoints, nSpots, t1t2, verbose, intAx);
end
%=======================


%=======================
%% POST-PROCESSING
%=======================

% loop again through idlist, and link the spots in frames where there is
% not the "good" number of tags. Then, sort and augment idlist and add
% linkup and linkdown.
if recalc < 2
    [idlist,maxTagIdx] = ...
        linkerSecondLAP(idlist, goodTimes, goodTimesM, nSpots, intList,...
        tagIndices, maxTagIdx, constants, verbose, intAx);
end



% for every tag that has not been linked: Make a positional estimate on the
% assumption of Brownian motion. Also, estimate intensities and a
% searchradius.
if recalc < 3
    [idlist] = linkerEstimatePositions(...
        idlist, maxTagIdx, goodTimes, constants, verbose, intAx);
end

% fill Q-matrices. Don't put in the uncertainties of the estimated tags yet
idlist = linkerWriteQmatrices(idlist,goodTimes);

% once fusions are implemented, we need this. (linkerFuseTags will need to
% call linkerEstimatePositions again!)
if doFusions
    [idlist] = linkerFuseTags(idlist, dataProperties, ...
        maxTagIdx, constants, verbose, intAx);
end

%=======================

%% DEBUG
if DEBUG
    for im=1:ceil(nTimepoints/25)
        figure
        for t=(im-1)*25+1:im*25,
            if t<=nTimepoints
                iplot = t-(im-1)*25;
                subplot(5,5,iplot),
                frame = movMax(:,:,:,:,t);
                imshow(frame,[]),
                hold on
                text(-1,-1,num2str(t))
                for i=1:maxTagIdx
                    if ~isempty(idlist(t).linklist)
                        if idlist(t).linklist(i,5) ~= 3
                            th=text(idlist(t).linklist(i,9)/constants.pix2mu(1),...
                                idlist(t).linklist(i,10)/constants.pix2mu(1),...
                                num2str(idlist(t).linklist(i,4)));
                            if idlist(t).linklist(i,2) > 0
                                set(th,'Color','r')
                            else
                                set(th,'Color','b')
                            end
                            switch idlist(t).linklist(i,5)
                                case 0
                                    % no marking
                                case 1
                                    % green dot
                                    plot(idlist(t).linklist(i,9)/constants.pix2mu(1),...
                                        idlist(t).linklist(i,10)/constants.pix2mu(1),'hm')
                                case 2
                                    % green x
                                    plot(idlist(t).linklist(i,9)/constants.pix2mu(1),...
                                        idlist(t).linklist(i,10)/constants.pix2mu(1),'xg')
                            end
                        end
                        plot(idlist(t).centroid(1)/constants.pix2mu(1),...
                            idlist(t).centroid(2)/constants.pix2mu(2),'*g')
                    end
                end

            end
        end
    end
end


%=======================
%% WRITE IDLIST
%=======================

% add help, start history, expand Q-matrices
% finishIdlist; % nested subfunctions are nice, but they make debugging
% hell
if recalc == -1

            % add help. Sort and expand Q-matrices where applicable. Remove unnecessary
            % fields
            %write explanation for linklist
            idlist(1).stats.help{1,1}='Columns of linklist:';
            idlist(1).stats.help{2,1}='1: timespot #';
            idlist(1).stats.help{3,1}='2: spot # (0: tag not found)';
            idlist(1).stats.help{4,1}='3: unused';
            idlist(1).stats.help{5,1}='4: tag color';
            idlist(1).stats.help{6,1}='5: flag';
            idlist(1).stats.help{7,1}='6: linkup to spot #';
            idlist(1).stats.help{8,1}='7: linkdown to spot #';
            idlist(1).stats.help{9,1}='8: intensity';
            idlist(1).stats.help{10,1}='9-11: x/y/z-coordinates in um (Image Coordinates!)';
            idlist(1).stats.help{11,1}='12: chi^2 of the spot on which the tag is located';
            idlist(1).stats.help{12,1}='';
            idlist(1).stats.help{13,1}='Flags in col 3:';
            idlist(1).stats.help{14,1}='1: estimated spot';
            idlist(1).stats.help{15,1}='2: MTM spot (found by tagTracker)';
            idlist(1).stats.help{16,1}='3/4: primary/secondary fusion spot';
            idlist(1).stats.help{17,1}='Flags in col 5:';
            idlist(1).stats.help{18,1}='1: adjusted intensity';
            idlist(1).stats.help{19,1}='2: only occurence of tag';
            idlist(1).stats.help{20,1}='3: estimated position of only occurence';

            %write labellist
            idlist(1).stats.labellist{1,1}='spb1';
            idlist(1).stats.labellist{2,1}='cen1';
            idlist(1).stats.labellist{3,1}='spb2';
            idlist(1).stats.labellist{4,1}='cen2';
        end

        if recalc == -1
            %write color<->label as all '?'
            idlist(1).stats.labelcolor(1:maxTagIdx,1)=cellstr('?');
        else
            % write color<->label and keep previous labels. MAXSPOTS could have
            % changed. If it increased, we will just remember another '?', if it
            % decreased, the good labels hopefully don't switch position
            oldLabels = ...
                idlist(1).stats.labelcolor(1:min(dataProperties.MAXSPOTS,...
                length(idlist(1).stats.labelcolor)));
            idlist(1).stats.labelcolor = {}; % necessary if maxspots decreased
            idlist(1).stats.labelcolor(1:maxTagIdx,1)=cellstr('?');
            idlist(1).stats.labelcolor(1:length(oldLabels)) = oldLabels;
        end


        % write maxColor (=maxIdx)
        idlist(1).stats.maxColor = maxTagIdx;

        if recalc == -1

            % write name, created, history
            idlist(1).stats.name = dataProperties.name;
            idlist(1).stats.created = date;
            idlist(1).stats.status = {[date ': idlist created']};
        else
            idlist(1).stats.status{end+1,1} = ...
                sprintf('%s : recalc %i',date,recalc);
        end

        % remove unnecessary fields
        if ~DEBUG
            if isfield(idlist,'distMatAmp')
                idlist = rmfield(idlist,...
                    {'distMatAmp','distMatXyz','distMat'});
                if isfield(idlist,'sourceIdxList')
                    idlist = rmfield(idlist,'sourceIdxList');
                end
            end
        end
        if isfield(idlist(1).stats,'recalc')
            stats = idlist(1).stats;
            stats = rmfield(stats, 'recalc');
            idlist(1).stats = stats;
        end





        % expand the Q-matrices. For the lost tags, add boundaries for tracker
        % (square them, so that we don't need to distinguish between lost tags and
        % found tags in the tracker). Since Q is diagonal, we could also make it a
        % sparse matrix



        for t=goodTimes'


            % write chi2 - account for possible fusions
            for i = 1:length(idlist(t).info.chi)
                spotIdx = idlist(t).linklist(:,2)==i;
                idlist(t).linklist(spotIdx,12) = idlist(t).info.chi(i);
            end

            % loop through trackInit and write boundaries (admittedly not the most
            % beautiful way, but it's Friday)
            for i=1:size(idlist(t).trackInit,1)
                % read tagIdx
                tagIdx = idlist(t).trackInit(i,1);
                % fill Q-matrix with sqare of the boundary
                idlist(t).info.detectQ_Pix(...
                    (tagIdx-1)*3+1:tagIdx*3,(tagIdx-1)*3+1:tagIdx*3) =...
                    diag(idlist(t).trackInit(i,2:end).^2);
            end

            % remove chi
            info = idlist(t).info;
            info = rmfield(info,{'chi'});
            idlist(t).info = info;

        end



% ADD SAVE OPTION HERE


%=======================




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NESTED SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function finishIdlist
% 
%         if recalc == -1
% 
%             % add help. Sort and expand Q-matrices where applicable. Remove unnecessary
%             % fields
%             %write explanation for linklist
%             idlist(1).stats.help{1,1}='Columns of linklist:';
%             idlist(1).stats.help{2,1}='1: timespot #';
%             idlist(1).stats.help{3,1}='2: spot # (0: tag not found)';
%             idlist(1).stats.help{4,1}='3: unused';
%             idlist(1).stats.help{5,1}='4: tag color';
%             idlist(1).stats.help{6,1}='5: flag';
%             idlist(1).stats.help{7,1}='6: linkup to spot #';
%             idlist(1).stats.help{8,1}='7: linkdown to spot #';
%             idlist(1).stats.help{9,1}='8: intensity';
%             idlist(1).stats.help{10,1}='9-11: x/y/z-coordinates in um (Image Coordinates!)';
%             idlist(1).stats.help{11,1}='12: chi^2 of the spot on which the tag is located';
%             idlist(1).stats.help{12,1}='';
%             idlist(1).stats.help{13,1}='Flags in col 3:';
%             idlist(1).stats.help{14,1}='1: estimated spot';
%             idlist(1).stats.help{15,1}='2: MTM spot (found by tagTracker)';
%             idlist(1).stats.help{16,1}='3/4: primary/secondary fusion spot';
%             idlist(1).stats.help{17,1}='Flags in col 5:';
%             idlist(1).stats.help{18,1}='1: adjusted intensity';
%             idlist(1).stats.help{19,1}='2: only occurence of tag';
%             idlist(1).stats.help{20,1}='3: estimated position of only occurence';
% 
%             %write labellist
%             idlist(1).stats.labellist{1,1}='spb1';
%             idlist(1).stats.labellist{2,1}='cen1';
%             idlist(1).stats.labellist{3,1}='spb2';
%             idlist(1).stats.labellist{4,1}='cen2';
%         end
% 
%         if recalc == -1
%             %write color<->label as all '?'
%             idlist(1).stats.labelcolor(1:maxTagIdx,1)=cellstr('?');
%         else
%             % write color<->label and keep previous labels. MAXSPOTS could have
%             % changed. If it increased, we will just remember another '?', if it
%             % decreased, the good labels hopefully don't switch position
%             oldLabels = ...
%                 idlist(1).stats.labelcolor(1:min(dataProperties.MAXSPOTS,...
%                 length(idlist(1).stats.labelcolor)));
%             idlist(1).stats.labelcolor = {}; % necessary if maxspots decreased
%             idlist(1).stats.labelcolor(1:maxTagIdx,1)=cellstr('?');
%             idlist(1).stats.labelcolor(1:length(oldLabels)) = oldLabels;
%         end
% 
% 
%         % write maxColor (=maxIdx)
%         idlist(1).stats.maxColor = maxTagIdx;
% 
%         if recalc == -1
% 
%             % write name, created, history
%             idlist(1).stats.name = dataProperties.name;
%             idlist(1).stats.created = date;
%             idlist(1).stats.status = {[date ': idlist created']};
%         else
%             idlist(1).stats.status{end+1,1} = ...
%                 sprintf('%s : recalc %i',date,recalc);
%         end
% 
%         % remove unnecessary fields
%         if ~DEBUG
%             if isfield(idlist,'distMatAmp')
%                 idlist = rmfield(idlist,...
%                     {'distMatAmp','distMatXyz','distMat'});
%                 if isfield(idlist,'sourceIdxList')
%                     idlist = rmfield(idlist,'sourceIdxList');
%                 end
%             end
%         end
%         if isfield(idlist(1).stats,'recalc')
%             stats = idlist(1).stats;
%             stats = rmfield(stats, 'recalc');
%             idlist(1).stats = stats;
%         end
% 
% 
% 
% 
% 
%         % expand the Q-matrices. For the lost tags, add boundaries for tracker
%         % (square them, so that we don't need to distinguish between lost tags and
%         % found tags in the tracker). Since Q is diagonal, we could also make it a
%         % sparse matrix
% 
% 
% 
%         for t=goodTimes'
% 
% 
%             % write chi2 - account for possible fusions
%             for i = 1:length(idlist(t).info.chi)
%                 spotIdx = idlist(t).linklist(:,2)==i;
%                 idlist(t).linklist(spotIdx,12) = idlist(t).info.chi(i);
%             end
% 
%             % loop through trackInit and write boundaries (admittedly not the most
%             % beautiful way, but it's Friday)
%             for i=1:size(idlist(t).trackInit,1)
%                 % read tagIdx
%                 tagIdx = idlist(t).trackInit(i,1);
%                 % fill Q-matrix with sqare of the boundary
%                 idlist(t).info.detectQ_Pix(...
%                     (tagIdx-1)*3+1:tagIdx*3,(tagIdx-1)*3+1:tagIdx*3) =...
%                     diag(idlist(t).trackInit(i,2:end).^2);
%             end
% 
%             % remove chi
%             info = idlist(t).info;
%             info = rmfield(info,{'chi'});
%             idlist(t).info = info;
% 
%         end
% 
%     end % nested function
% end % main function