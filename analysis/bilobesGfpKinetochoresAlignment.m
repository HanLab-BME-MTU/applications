function [spbCorrection,alignmentVectors] = bilobeSpindlesAlignment
%BILOBESPINDLES calculates the aligment between two channels
%
% SYNOPSIS: bilobeSpindlesAlignment
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

% search files
[fileList,tokens]=searchFiles('(\w+)_R3D\.dv','log','ask');

% select files
selection = listSelectGUI(tokens);

% only keep selection
fileList = fileList(selection,:);
tokens = tokens(selection);

nData = length(selection);

nameList(1:nData) = struct('rawMovieName','','filteredMovieName','','slistName','',...
    'dataPropertiesName','','testRatiosName','');

% resultList(1:nData) = struct('interpImg','','maxProj','','meanProj','',...
%     'e_spb','','e_perp','','e_3','','n_spb','');

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
    nameList(iData).movieHeaderName = ...
        [fileList{iData,2},filesep,'movieHeader_',tokens{iData},'.mat'];
    nameList(iData).filteredMovieName{1} = ...
        [fileList{iData,2},filesep,'filtered_w1_',tokens{iData},'.fim'];
    nameList(iData).slistName{1} = ...
        [fileList{iData,2},filesep,'slist_w1_',tokens{iData},'.mat'];
    nameList(iData).idlistName{1} = ...
        [fileList{iData,2},filesep,'idlist_w1_',tokens{iData},'.mat'];
    nameList(iData).idlist_LName{1} = ...
        [fileList{iData,2},filesep,'idlist_L_w1_',tokens{iData},'.mat'];
    nameList(iData).dataPropertiesName{1} = ...
        [fileList{iData,2},filesep,'dataProperties_w1_',tokens{iData},'.mat'];
    nameList(iData).filteredMovieName{2} = ...
        [fileList{iData,2},filesep,'filtered_w2_',tokens{iData},'.fim'];
    nameList(iData).slistName{2} = ...
        [fileList{iData,2},filesep,'slist_w2_',tokens{iData},'.mat'];
    nameList(iData).idlistName{2} = ...
        [fileList{iData,2},filesep,'idlist_w2_',tokens{iData},'.mat'];
    nameList(iData).idlist_LName{2} = ...
        [fileList{iData,2},filesep,'idlist_L_w2_',tokens{iData},'.mat'];
    nameList(iData).dataPropertiesName{2} = ...
        [fileList{iData,2},filesep,'dataProperties_w2_',tokens{iData},'.mat'];
    nameList(iData).dirName = fileList{iData,2};

    % load movie - use low-level r3dread
    rawMovieFull = r3dread(nameList(iData).rawMovieName);
    movieHeader = readr3dheader(nameList(iData).rawMovieName);

    % save movieHeader
    save(nameList(iData).movieHeaderName,'movieHeader');


    % loop wavelengths
    for iw = 1:2

        % if slist already exists, don't filter etc. again
        if exist(nameList(iData).slistName{iw},'file')
            load(nameList(iData).slistName{iw})
            load(nameList(iData).dataPropertiesName{iw})
            
            %--------- temporary bugfix --- add dataProperties.name
            dataProperties.name = tokens{iData};
            save(nameList(iData).dataPropertiesName{iw},'dataProperties');
            %------------------------
            
            
            filteredMovie = readmat(nameList(iData).filteredMovieName{iw});
        else

            % select spb channel
            % spb channel is the first channel
            rawMovie = rawMovieFull(:,:,:,iw);

            % make data properties for the first channel
            dataProperties = defaultDataProperties(movieHeader);
            dataProperties.waveIdx = iw; % spb channel
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
            dataProperties.name = tokens{iData};

            % save data Properties
            save(nameList(iData).dataPropertiesName{iw},'dataProperties');


            % filter movie
            filteredMovie = filtermovie(rawMovie,dataProperties.FILTERPRM);

            % save filtered movie
            writemat(nameList(iData).filteredMovieName{iw},...
                filteredMovie);

            % detect two spots


            [slist, dataProperties, testRatios] = detectSpots(...
                nameList(iData).rawMovieName, ...
                nameList(iData).filteredMovieName{iw}, ...
                dataProperties,2);

            % save results
            save(nameList(iData).slistName{iw},'slist');
            save(nameList(iData).dataPropertiesName{iw},'dataProperties')
            %save(nameList(iData).testRatiosName,'testRatios');

        end

        if exist(nameList(iData).idlist_LName{iw},'file')
            load(nameList(iData).idlist_LName{iw})
            idlist = idlist_L;
        else

            % make idlist, save
            idlist = linker(slist,dataProperties,1);
            save(nameList(iData).idlistName{iw},'idlist');
        end


        % check idlist. Change to load both channels, spots in Imaris later
        % if necessary
        if ~isempty(idlist.linklist)
            % if empty idlist, just don't save a idlist_L. 
        lh = LG_loadAllFromOutside(filteredMovie,nameList(iData).dirName,...
            [],dataProperties,idlist,'idlist');
        uiwait(lh)
        idlist_L = LG_readIdlistFromOutside;

        % Remove these lines if we check for number of spots later
        %         idlist_L = idlist;
        save(nameList(iData).idlist_LName{iw},'idlist_L');
        end

    end

end

alignmentVectors = nan(nData*2,3);
% loop through data, find offset
for iData = 1:nData
    % load idlist_L - remember coordinate flip!
    try
    load(nameList(iData).idlist_LName{1});
    coord1 = idlist_L.linklist(:,[10,9,11]);
    load(nameList(iData).idlist_LName{2});
    coord2 = idlist_L.linklist(:,[10,9,11]);
    % find correspondence
    [dist,distVec] =  distMat2(coord1,coord2);
    [oneTwo] = lap(dist,[],[],1);
    nc1 = size(coord1,1);
    nc2 = size(coord2,1);
    % there will always be at least 2 vectors, even if they're nans
    tmpVectors = nan(max(nc1,2),3);
    tmpDist = inf(max(nc1,2),1);
    for i=1:nc1
        % read out links
        if oneTwo(i) <= nc2
            tmpVectors(i,:) = squeeze(distVec(i,oneTwo(i),:))';
            tmpDist(i) = dist(i,oneTwo(i));
        end
    end
    % take two shortest vectors
    [dummy,sortIdx] = sort(tmpDist);
    alignmentVectors((iData-1)*2+1:iData*2,:) = tmpVectors(sortIdx(1:2),:);
    catch
    end
end % loop through data to find short vectors
alignmentVectors(any(~isfinite(alignmentVectors),2),:) = [];
% calculate spbCorrection
spbCorrection = robustMean(alignmentVectors);





