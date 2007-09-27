function sigmaCorrections = bilobesGfpKinetochoresFitPsf
% fit of Gauss to PSF of raw/decon image pairs of SPBs



% search files
[fileList,tokens]=searchFiles('(\w+)_R3D\.dv','log','ask');

% select files
selection = listSelectGUI(tokens);

% only keep selection
fileList = fileList(selection,:);
tokens = tokens(selection);

nData = length(selection);

nameList(1:nData,2) = struct('rawMovieName','','filteredMovieName','','slistName','',...
    'dataPropertiesName','','testRatiosName','');

sigmaCorrections(1:nData) = struct('rawGfp',[],'rawCfp',[],...
    'deconGfp',[],'deconCfp',[]);
scNames = reshape(fieldnames(sigmaCorrections),2,2);
% add more fields
sigmaCorrections(1).dataName = [];
for i=1:4
    sigmaCorrections(1).([scNames{i} 'Img']) = [];
end

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
    % raw movie
    nameList(iData,1).rawMovieName = fullfile(fileList{iData,2},fileList{iData,1});
    nameList(iData,1).movieHeaderName = ...
        [fileList{iData,2},filesep,'movieHeader_',tokens{iData},'.mat'];
    nameList(iData,1).filteredMovieName{1} = ...
        [fileList{iData,2},filesep,'filtered_w1_',tokens{iData},'.fim'];
    nameList(iData,1).slistName{1} = ...
        [fileList{iData,2},filesep,'slist_w1_',tokens{iData},'.mat'];
    nameList(iData,1).idlistName{1} = ...
        [fileList{iData,2},filesep,'idlist_w1_',tokens{iData},'.mat'];
    nameList(iData,1).idlist_LName{1} = ...
        [fileList{iData,2},filesep,'idlist_L_w1_',tokens{iData},'.mat'];
    nameList(iData,1).dataPropertiesName{1} = ...
        [fileList{iData,2},filesep,'dataProperties_w1_',tokens{iData},'.mat'];
    nameList(iData,1).filteredMovieName{2} = ...
        [fileList{iData,2},filesep,'filtered_w2_',tokens{iData},'.fim'];
    nameList(iData,1).slistName{2} = ...
        [fileList{iData,2},filesep,'slist_w2_',tokens{iData},'.mat'];
    nameList(iData,1).idlistName{2} = ...
        [fileList{iData,2},filesep,'idlist_w2_',tokens{iData},'.mat'];
    nameList(iData,1).idlist_LName{2} = ...
        [fileList{iData,2},filesep,'idlist_L_w2_',tokens{iData},'.mat'];
    nameList(iData,1).dataPropertiesName{2} = ...
        [fileList{iData,2},filesep,'dataProperties_w2_',tokens{iData},'.mat'];
    nameList(iData,1).dirName = fileList{iData,2};

    % decon movie
    % find decon movie name
    deconMovieName = searchFiles('D3D.dv',[],fileList{iData,2});
    nameList(iData,2).rawMovieName = fullfile(fileList{iData,2},deconMovieName{1});
    nameList(iData,2).movieHeaderName = ...
        [fileList{iData,2},filesep,'movieHeader_',tokens{iData},'_decon.mat'];
    nameList(iData,2).filteredMovieName{1} = ...
        [fileList{iData,2},filesep,'filtered_w1_',tokens{iData},'_decon.fim'];
    nameList(iData,2).slistName{1} = ...
        [fileList{iData,2},filesep,'slist_w1_',tokens{iData},'_decon.mat'];
    nameList(iData,2).idlistName{1} = ...
        [fileList{iData,2},filesep,'idlist_w1_',tokens{iData},'_decon.mat'];
    nameList(iData,2).idlist_LName{1} = ...
        [fileList{iData,2},filesep,'idlist_L_w1_',tokens{iData},'_decon.mat'];
    nameList(iData,2).dataPropertiesName{1} = ...
        [fileList{iData,2},filesep,'dataProperties_w1_',tokens{iData},'_decon.mat'];
    nameList(iData,2).filteredMovieName{2} = ...
        [fileList{iData,2},filesep,'filtered_w2_',tokens{iData},'_decon.fim'];
    nameList(iData,2).slistName{2} = ...
        [fileList{iData,2},filesep,'slist_w2_',tokens{iData},'_decon.mat'];
    nameList(iData,2).idlistName{2} = ...
        [fileList{iData,2},filesep,'idlist_w2_',tokens{iData},'_decon.mat'];
    nameList(iData,2).idlist_LName{2} = ...
        [fileList{iData,2},filesep,'idlist_L_w2_',tokens{iData},'_decon.mat'];
    nameList(iData,2).dataPropertiesName{2} = ...
        [fileList{iData,2},filesep,'dataProperties_w2_',tokens{iData},'_decon.mat'];
    nameList(iData,2).dirName = fileList{iData,2};
    
    sigmaCorrections(iData).dataName = tokens{iData};
    
    % loop raw, decon
    for movieType = 1:2

        % load movie - use low-level r3dread
        rawMovieFull = r3dread(nameList(iData,movieType).rawMovieName);
        movieHeader = readr3dheader(nameList(iData,movieType).rawMovieName);

        % save movieHeader
        save(nameList(iData,movieType).movieHeaderName,'movieHeader');


        % loop wavelengths
        for iw = 1:2

            % if slist already exists, don't filter etc. again
            if exist(nameList(iData,movieType).slistName{iw},'file')
                load(nameList(iData,movieType).slistName{iw})
                load(nameList(iData,movieType).dataPropertiesName{iw})

                %--------- temporary bugfix --- add dataProperties.name
                dataProperties.name = tokens{iData};
                save(nameList(iData,movieType).dataPropertiesName{iw},'dataProperties');
                %------------------------


                %filteredMovie = readmat(nameList(iData,movieType).filteredMovieName{iw});
            else

                
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
                save(nameList(iData,movieType).dataPropertiesName{iw},'dataProperties');


                % filter movie
                filteredMovie = filtermovie(rawMovie,dataProperties.FILTERPRM);

                % save filtered movie
                writemat(nameList(iData,movieType).filteredMovieName{iw},...
                    filteredMovie);

                % detect two spots


                [slist, dataProperties, testRatios] = detectSpots(...
                    nameList(iData,movieType).rawMovieName, ...
                    nameList(iData,movieType).filteredMovieName{iw}, ...
                    dataProperties,2);

                % save results
                save(nameList(iData,movieType).slistName{iw},'slist');
                save(nameList(iData,movieType).dataPropertiesName{iw},'dataProperties')
                %save(nameList(iData,movieType).testRatiosName,'testRatios');

            end

            % we may need an idlist for display purposes
            if exist(nameList(iData,movieType).idlist_LName{iw},'file')

            else

                % make idlist, save
                idlist = linker(slist,dataProperties,1);
                save(nameList(iData,movieType).idlistName{iw},'idlist');
            end
            
            % figure out the current wavelength
            % make it simple for right now: >500 is GFP
            currentScName = scNames{(dataProperties.WVL(iw)<0.5)+1,movieType};
            
            % look at slist. Only keep spot with highest intensity
            if length(slist.sp) > 1
                amp = catStruct(1,'slist.sp.amp');
                [a,idx] = max(amp);
                slist.sp = slist.sp(idx);
            end
            
            % fit PSF with Gauss
            fitStruct = struct(...
                'slist',slist,...
                'dataProperties',dataProperties,...
                'movieDir',fileList{iData,2},...
                'rawMovieName',nameList(iData,movieType).rawMovieName,...
                'waveIdx',iw);
            
            
            [fitStruct,fitImages] = cdFitPsf(fitStruct);
            
            % read sigmaCorrection
            sigmaCorrections(iData).(currentScName) = ...
                fitStruct.dataProperties.sigmaCorrection;
            sigmaCorrections(iData).([currentScName, 'Img']) = fitImages;
            
            
                


        end % loop color

    end % loop raw/decon
    
    % display: gfpRaw,cfpRaw,gfpDecon,cfpDecon
    % top: raw, bottom av+fit
    figure('Name',sigmaCorrections(iData).dataName);
    for i=1:4
        % read fitImage
        fitImages = sigmaCorrections(iData).([scNames{i}, 'Img']);
        % top: raw
        if ~isempty(fitImages)
        topIm = cat(1,nanmean(fitImages{1},3),...
            squeeze(nanmean(fitImages{1},1))');
        red = cat(1,nanmean(fitImages{2},3),...
            squeeze(nanmean(fitImages{2},1))');
        maxRed = max(red(:));
        green = cat(1,nanmean(fitImages{3},3),...
            squeeze(nanmean(fitImages{3},1))');
        blue = zeros(size(topIm));
        botIm = cat(3,red,green,blue)/maxRed;
        subplot(2,4,i)
        imshow(topIm,[]);
        title(scNames{i})
        subplot(2,4,i+4)
        imshow(botIm)
        end
        
    end
        

end % loop data