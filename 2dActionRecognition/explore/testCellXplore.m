% % script/steps to test/run cellXploreDR
% 
% 
% %% how to close out figures -- 
% disp('Run the following to close all figures')
disp('ff = findall(0,''Type'', ''Figure''); close(ff);');
% 
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Test cellXplorer');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Gen & Load Fake Data');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% 
% load('mockData.mat');
% genFakeMovies;
% 
% % playMovie(movies{1})
% cellXploreDR(data, 'extra', [], 'movies', movies);
% 
% % Example how to add new DR view choices
% data.DR.rand = rand([150, 2])


% get data from assaf (DONE)
% select sub-set of interest
% generate struct array
% make movie cell based on struct array. (optional to create MDs)
% ---
% input into cellX
% cellX parse new data struct.
% annotate at will
% update struct data structure.
% calculate DR based on features

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Get Assaf''s Data');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

allCellsFileName = '28-Mar-2017_LBP_dLBP_1.mat';
InFilePath =  '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData';
loadMatPath = fullfile(InFilePath,allCellsFileName);

% My output path
% devOpsDir = '/work/bioinformatics/s170480/Data/LCH/DevOps';

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Load Assaf''s Data');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

load(loadMatPath, 'allCellsMovieData'); 
% extractMoviesToCellExplorer

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Select Subset');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

randset = 50;

if isempty(randset)
    % Get master movie list
    exprMD = cellfun(@(x) x.expStr, allCellsMovieData,'UniformOutput',false);
    uExprMD = unique(exprMD);
    sum(ismember(cellfun(@(x) x.expStr, allCellsMovieData,'UniformOutput',false), uExprMD{2}));

    % Grab a single experiment
    singleExprMD = allCellsMovieData(1,ismember(exprMD, uExprMD{1:2}));
    disp(['Experiment movie: ' singleExprMD{1}.expStr])
%     cellDataSet = cell2mat(singleExprMD);
else
    
%     cellDataSet = cell2mat(allCellsMovieData(1,randsample(length(allCellsMovieData),randset,false)));
    singleExprMD = allCellsMovieData(1,randsample(length(allCellsMovieData),randset,false));
end
    
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Create struct array');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


% cellDataSet = cellDataSet(1:3);
% t=table({cellDataSet.key}')


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Load Assaf''s MDs');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% load movies into memory



nCells = numel(singleExprMD);
cellMovies = cell(1,nCells);

for i = 1:nCells

    iCell = singleExprMD{i};
    disp(iCell.key);
  
    xs = iCell.xs;
    ys = iCell.ys;
    ts = iCell.ts;
    
    ntime = length(ts);
    
    MD = load(iCell.MD, 'MD');
    MD = MD.MD;
    
    %% HARD CODED!
    pixelSize_ = 0.325;
    FOVRadius = round(35/pixelSize_);
    movie = zeros(2*FOVRadius+1, 2*FOVRadius+1, ntime);
      
%     R = MD.getReader();
%     cR = CellReader(R);
%     curI = cR(1,ts);
    
    parfor itime = 1:ntime
        curTime = ts(itime);
        curI = MD.getChannel(1).loadImage(curTime);

        movie(:,:,itime) = curI(round((ys(itime)-FOVRadius)):(round(ys(itime)+FOVRadius)),round((xs(itime)-FOVRadius)):round((xs(itime)+FOVRadius)));
    end
    cellMovies{i} = movie;
end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Crop Assaf''s master MDs & create cellMDs');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% 
% 
% load('/work/bioinformatics/s170480/Data/LCH/DevOps/testMovies');
%   
save(['/work/bioinformatics/s170480/Data/LCH//DevOps/testMoviesRand_' num2str(randset)], 'cellMovies','cellDataSet');

cellDataSet = cell2mat(singleExprMD);
% cellXploreDR(cellDataSet, cellMovies);
cellXploreDR(cellDataSet, cellMovies, 'annotations', cellDataSet(1).annotations)

% %% how to close out figures -- 
% disp('Run the following to close all figures')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('ff = findall(0,''Type'', ''Figure''); close(ff);');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');






% 
% 
% 
% %%% Description Script for testing new plusTipsTrackerPackage3D
% 
% if lower(string(computer('arch'))) == 'glnxa64'  % matlab is getting more pythonic!
%     % Analysis Output Directory
%     analysis_dir = '/work/bioinformatics/s170480/Analysis/3D/Tracking/A1_HeLa_EB1/scriptTest';
%     % Raw Input Image paths ('Original Data set here: /project/bioinformatics/Danuser_lab/shared/proudot/3d-vis/utrackPackage/scriptBased/UTrack-QD-v1/data/Ce...' <Preview truncated at 128 characters>)
%     ch0 = '/work/bioinformatics/s170480/Data/3D/Tracking/A1_HeLa_Cells_EB1/ch0';
%     ch1 = '/work/bioinformatics/s170480/Data/3D/Tracking/A1_HeLa_Cells_EB1/ch1';
% 
% elseif  lower(string(computer('arch'))) == 'win64'
%     % Analysis Output Directory
%     analysis_dir = 'C:\Users\Andrew\Analysis\3D\Tracking\A1_HeLa_EB1\scriptTest';
%     % Raw Input Image paths ('Original Data set here: /project/bioinformatics/Danuser_lab/shared/proudot/3d-vis/utrackPackage/scriptBased/UTrack-QD-v1/data/Ce...' <Preview truncated at 128 characters>)
%     ch0 = 'C:\Users\Andrew\Data\raw\3D\Tracking\A1_HeLa_EB1\ch0';
%     ch1 = 'C:\Users\Andrew\Data\raw\3D\Tracking\A1_HeLa_EB1\ch1';
% end
% 
% %----------------------------------------
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Initializing Channels                %%');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% % Retrieve current location
% % Channel creation
% % Create a channels object
% channel(1) = Channel(ch0);
% % channel.fluorophore_='alexa647';
% % channel.emissionWavelength_=name2wavelength('alexa647')*1e9;
% % channel.imageType_='TIRF';
% 
% channel(2) = Channel(ch1);
% % channel(2).fluorophore_='';
% % channel(2).emissionWavelength_=name2wavelength('egfp')*1e9;
% % channel(2).imageType_='';
% %----------------------------------------
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp('% Creating MovieData                  %%');
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% % MovieData creation
% outputDir = analysis_dir; %fileparts(analysis_dir);
% MD = MovieData(channel, outputDir);
% 
% movieDataFileName = 'movieDataScriptTest.mat';
% MD.setPath(outputDir);
% MD.setFilename(movieDataFileName);
% 
% % Imaging/Microscope/Movie Parameters
% % MD.numAperture_= 1.49; 
% MD.pixelSize_= 100;
% MD.pixelSizeZ_ = 235;
% MD.camBitdepth_= 16;
% MD.timeInterval_ = 1;
% MD.notes_= 'EB3 MT tracking 3D testing movie with new plusTipsTrackerPackage3D'; 
% MD.sanityCheck;
% MD.save;
% 
% % Load the movie
% clear MD;
% MD = MovieData.load(fullfile(outputDir, movieDataFileName));
% MD.reset();
% 
% 
% 
% 
% 
% 
% @(x,y)PlusTipTrackerPackage3D.AmiraExportExtProcess(x,y, ...
%                            'RefProcess', x.processes_{x.getProcessIndex('CreateReferenceFrameProcess')}, ... % note will ask user if ambiguous
%                            'filename', [filesep 'amiraVertexLabRef' filesep  'detectionLabRef.am']), ... 
% 
%        
% function extProc = AmiraExportExtProcess(owner, outputDir, varargin)
% %%%%%%%%%%%%%%%%%
% ip = inputParser;
% ip.addRequired('owner', @(x) isa(x,'MovieData'));
% ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
% ip.addParameter('RefProcess', owner.processes_{owner.getProcessIndex('CreateReferenceFrameProcess')}, @(x) isa(x,'Process'));
% ip.addParameter('filename', 'amiraVertex.am', @ischar);
% 
% ip.parse(owner, outputDir, varargin{:})
% 
% outputDir = ip.Results.outputDir;
% filename = ip.Results.filename;
% MD = ip.Results.owner;
% 
% funParams = struct();
% funParams.OutputDirectory = [outputDir  filesep 'Amira'];
% funParams.filename = [funParams.OutputDirectory filesep filename];
% 
% % load movieInfo 
% funParams.movieInfoPath = MD.processes_{MD.getProcessIndex('CreateReferenceFrameProcess')}.outFilePaths_;
% funParams.movieInfo = load(MD.processes_{MD.getProcessIndex('CreateReferenceFrameProcess')}.outFilePaths_);            
% funParams.movieInfo = funParams.movieInfo.detectionsStageRef;
% 
% dataAnisotropy = [MD.pixelSize_ MD.pixelSize_ MD.pixelSizeZ_];
% dataIsotropy = [MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
% funParams.scales = dataIsotropy;
% funParams.prop = {};
% 
% % Create External Process
% extParams = struct();
% extParams.name = 'AmiraExport';
% extParams.fun = @(pr, varargin) amiraWriteMovieInfo(pr.getParameters().filename, ...
%                                                     pr.getParameters().movieInfo, ...
%                                                     varargin{:});
% extParams.parameters = funParams;
% extProc = ExternalProcess(owner, extParams);
% 
% extProc.setInFilePaths(funParams.movieInfoPath); 
% extProc.setOutFilePaths(funParams.filename);
% 
% %%%%%%%%%%%%%%%%%
