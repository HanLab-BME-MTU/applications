% /project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All

% Input Path
allCellsFileName = '12-Nov-2017_LBP_dLBP_1.mat';
InFilePath =  '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData';
loadMatPath = fullfile(InFilePath, allCellsFileName);


matFile_name = 'FullTimeSeries_Gen2n3_12Nov2017';

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Chop up Assaf''s MDs');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

load(loadMatPath);%, 'allCellsMovieData'); 

randset = {};

if isempty(randset)
    allCellsSet = allCellsMovieData;
    % Output path
    omeTiffDir = '/work/bioinformatics/shared/dope/data/OMETIFF/fulltime_Gen2n3_12-Nov-2017';
    if ~exist(omeTiffDir, 'dir')
        mkdir(omeTiffDir);
    end
else
    %     cellDataSet = cell2mat(allCellsMovieData(1,randsample(length(allCellsMovieData),randset,false)));
    allCellsSet = allCellsMovieData(1,randsample(length(allCellsMovieData),randset,false));
    % Output path
    omeTiffDir = ['/work/bioinformatics/shared/dope/data/OMETIFF/test/R' num2str(randset)];
    if ~exist(omeTiffDir, 'dir')
        mkdir(omeTiffDir);
    end

end

javaaddpath('/home2/s170480/matlab/extern/bioformats/bioformats_package.jar','-end')
MList = cell(1, length(allCellsSet));

% pixelSize = ome.units.quantity.Length(java.lang.Double(.325), ome.units.UNITS.MICROMETER);
PIXS = @(y) ome.units.quantity.Length(java.lang.Double(y), ome.units.UNITS.MICROMETER);
pixelSizepf = arrayfun(@(x) PIXS(x), repmat(.325,[1 length(allCellsSet)]), 'UniformOutput', false);


% How many frame to include in the movie chops...
% timeChop = 10?




checkFailList = cell(1,length(allCellsSet));
parfor i = 1:length(allCellsSet)

    iCell = allCellsSet{i};
    
    xs = iCell.xs;
    ys = iCell.ys;
    ts = iCell.ts;
    
    startTime = ts;
    ntime = length(ts);
    endTime = ts(1) + ntime;

    iCell.key_noTime = iCell.key;
    iCell.key = [iCell.key '_t' num2str(endTime)]

%     disp(iCell.key);
    
    try 
        MD = load(iCell.MD, 'MD');
        MD = MD.MD;

        
        
        % first check if already created
        movieFileOutMat = [omeTiffDir filesep iCell.key filesep iCell.key '_CellX.mat'];
        
        allCellsSet{i}.cellMD = MD.getFullPath();
        allCellsSet{i}.key = iCell.key;
        allCellsSet{i}.key_noTime = iCell.key_noTime;
        
        if exist(movieFileOutMat,'file') == 2
            % check if passes sanity check
            MD = load(movieFileOutMat,'MD');
            MD = MD.MD;
            try
                MD.sanityCheck();
                try
                    SDCindx = MD.getProcessIndex('EfficientSubpixelRegistrationProcess');
                    if ~isempty(SDCindx) && MD.processes_{SDCindx}.success_
                        checkFailList{i} = {'done', allCellsSet{i}};
                    else
                        checkFailList{i} = {'SDCincomplete', allCellsSet{i}};
                    end
                catch
                    checkFailList{i} = {'getProcessFail', allCellsSet{i}};
                end                
            catch
                checkFailList{i} = {'fail', allCellsSet{i}};
                %                 disp('running MD:' mov
            end 
%             checkFailList{i} = {'done',allCellsSet{i}}
        else
            checkFailList{i} = {'noExist',allCellsSet{i}};
        end
    catch
        checkFailList{i} = {'badMD',allCellsSet{i}};
    end
end

matFile = '/work/bioinformatics/shared/dope/data/OMETIFF/fulltime_Gen2n3_12-Nov-2017/FullTimeSeries_Gen2n3_12Nov2017.mat'
load(matFile, 'cellDataSet');
omeTiffDir = '/work/bioinformatics/shared/dope/data/OMETIFF/fulltime_Gen2n3_12-Nov-2017';
parfor i = 1:length(cellDataSet)

    iCell = cellDataSet{i};
    
    xs = iCell.xs;
    ys = iCell.ys;
    ts = iCell.ts;
    
    startTime = ts;
    ntime = length(ts);
    endTime = ts(1) + ntime;

%     iCell.key_noTime = iCell.key;
%     iCell.key = [iCell.key '_t' num2str(endTime)]
    
    % first check if already created
    movieFileOutMat = [omeTiffDir filesep iCell.key filesep iCell.key '_CellX.mat'];

    cellDataSet{i}.cellMD = movieFileOutMat;
%     cellDataSet{i}.key = iCell.key;
%     cellDataSet{i}.key_noTime = iCell.key_noTime;   
    
end

matFile2 = '/work/bioinformatics/shared/dope/data/OMETIFF/fulltime_Gen2n3_12-Nov-2017/FullTimeSeries_Gen2n3_12Nov2017_corrected.mat'
save(matFile2, 'cellDataSet');


load('./checkList.mat');
javaaddpath('/home2/s170480/matlab/extern/bioformats/bioformats_package.jar','-end')
parfor i = 1:length(allCellsSet)

    iCell = allCellsSet{i};
    
    xs = iCell.xs;
    ys = iCell.ys;
    ts = iCell.ts;
    
    startTime = ts;
    ntime = length(ts);
    endTime = ts(1) + ntime;

    iCell.key_noTime = iCell.key;
    iCell.key = [iCell.key '_t' num2str(endTime)]

%     disp(iCell.key);
     
    try 
        MD = load(iCell.MD, 'MD');
        MD = MD.MD;
        if ~strcmp(checkFailList{i}{1},'done')

            %% HARD CODED!
            pixelSize_ = 0.325;
            FOVRadius = round(35/pixelSize_);
            movieM = zeros(2*FOVRadius+1, 2*FOVRadius+1, ntime);

            %% Configure OME-TIFF metadata
            metadata = createMinimalOMEXMLMetadata(movieM, 'XYTZC');
            metadata.setPixelsPhysicalSizeX(pixelSizepf{i}, 0);
            metadata.setPixelsPhysicalSizeY(pixelSizepf{i}, 0);
            metadata.setImageDescription(iCell.key, 0);
            metadata.setExperimenterGroupDescription(iCell.expStr,0);
            metadata.setExperimenterGroupID(iCell.date,0);
            metadata.setDatasetName(iCell.expStr,0);
            metadata.setDatasetID(iCell.MD,0);
            metadata.setDatasetDescription([...
                'notes=''for Gen2and3 Nov02nd 2017-ARJ''', 'key=''' iCell.key ''', key_noTime=''' iCell.key_noTime ''',date=''' iCell.date ''',Celltype=''' iCell.cellType ''',metEff=' num2str(iCell.metEff) ',locationStr=' num2str(iCell.locationStr)],0);

            for itime = 1:ntime
                curTime = ts(itime);
                curI = MD.getChannel(1).loadImage(curTime);
                movieM(:,:,itime) = curI(round((ys(itime)-FOVRadius)):(round(ys(itime)+FOVRadius)),round((xs(itime)-FOVRadius)):round((xs(itime)+FOVRadius)));
            end

            disp(['Creating movieData for cell ID key: ' iCell.key]);
            movieFileOut = [omeTiffDir filesep iCell.key filesep iCell.key '_CellX.ome.tiff'];
            disp(['Making OME-TIFF : ' movieFileOut]);

            % Save as OME-TIFF 
            mkdir(fileparts(movieFileOut));
            bfsave(movieM, movieFileOut, 'metadata', metadata);

            % MovieData Validation
            disp('Loading with MovieData BF reader')
            MD = MovieData(movieFileOut, true, fileparts(movieFileOut));
            MD.sanityCheck();
            MD.save();

            extProc = ExternalProcess(MD, 'LBPfeatures');
            MD.addProcess(extProc);
            MD.processes_{1}.setParameters(iCell);

            MD.addProcess(EfficientSubpixelRegistrationProcess(MD));
            SDCindx = MD.getProcessIndex('EfficientSubpixelRegistrationProcess');
            MD.processes_{SDCindx}.run()
        end
            % MDpath = MD.getFullPath();
            % Save individual cell movieData
            allCellsSet{i}.cellMD = MD.getFullPath();
            allCellsSet{i}.key = iCell.key;
            allCellsSet{i}.key_noTime = iCell.key_noTime;
            MList{i} = MD;
            checkFail{i} = true;
    catch ME
        disp(['FAILED MD{i} Load : ' num2str(i)]);
        disp(['FAILED MD Load : ' iCell.MD]);
        disp(ME);
        checkFail{i} = false;
    end
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Save cell array with cell MD info');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% save([omeTiffDir filesep 'clickFuryCellData.mat'], 'allCellsSet');

annotationSet = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells (by local master index)
annotationSet('Blebbing') = [NaN];
annotationSet('Balled') = [NaN];
annotationSet('Spread') = [NaN];
annotationSet('Smooth') = [NaN];
annotationSet('Short extension') = [NaN];
annotationSet('Long extension') = [NaN];
annotationSet('Migration') = [NaN];
annotationSet('Elongated') = [NaN];
annotationSet.keys

% Cell classes / labels / annotations:
% Blebbing- membrane ruffling 
% Balled- prominent halo, circular
% Spread- increased cell area, increased cytoplasm?
% Smooth- no discernible membrane ruffles
% Short extension- protrusion close to cell body
% Long extension- protrusion extending far from cell body
% Migration - cell movement that results in a change in xy
% Elongated- high eccentricity, long, thin
% Exclusion classifier:
% Dead / crap 
% Out of focus
% Multiple cells
% undefined

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Save .mat file with expected data structures and names (cell array and dictionary)');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

cellDataSet = allCellsSet;
if isempty(randset)
    matFileName = [omeTiffDir filesep matFile_name '.mat'];
else
    matFileName = [omeTiffDir filesep matFile_name '_R' num2str(randset) '.mat'];
end

save(matFileName, 'cellDataSet', 'annotationSet');

% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Examples to acccess metadata via MovieData/MovieList');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% uiwait

% MList{3}.getReader.formatReader.getMetadataStore.getDatasetID(0)
% MList{3}.getReader.formatReader.getMetadataStore.getDatasetName(0) 
% MList{3}.getReader.formatReader.getMetadataStore.getImageDescription(0) 
% MList{3}.getReader.formatReader.getMetadataStore.getExperimenterGroupDescription(0)
% MList{3}.getReader.formatReader.getMetadataStore.getDatasetDescription(0)

% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% How to run cellXplorerDR and dopeAnnotator');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% disp('Opt 1: load aboved saved .mat file by providing the file when the program is called');
% cellXplorerDR(matFileName);
% 
% uiwait
% % ff = findall(0,'Tag','cellXplore'); delete(ff);
% disp('Opt 2: call function immediately and be prompted (via UI) to manually select file to load');
% cellXplorerDR
% 
% uiwait
% 
% disp('Same applies for dopeAnnotator')
% dopeAnnotator(matFileName);
% uiwait
% dopeAnnotator

% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% [optional & slow..] Create MovieList');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% % ML = MovieList([MList{:}], omeTiffDir, 'movieListFileName_', 'movieListCells.mat');
% % ML.sanityCheck();
% % ML.save();



