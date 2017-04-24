 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% ====== Test cellXplorer =============');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Run stage drift correction example');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% first create the movie data...
% MD = MovieData.load('/work/bioinformatics/s170480/Data/LCH/DevOps/160808_m498_m634_m610_m481/28-Mar-2017_m610_s13_t95_x1162_y1122_CellX.mat')
% MD.addProcess(EfficientSubpixelRegistrationProcess(MD))
% MD.processes_{1}.run()
% movieViewer(MD)

% MD=MovieData.load('/work/bioinformatics/s170480/Data/LCH/DevOps/160808_m498_m634_m610_m481/28-Mar-2017_m481_s17_t95_x1344_y985_CellX.mat')
% movieViewer(MD)
% MD.addProcess(EfficientSubpixelRegistrationProcess(MD))
% MD.processes_{1}.run()
% movieViewer(MD)

% Via the GUI
% MD.addProcess(EfficientSubpixelRegistrationProcess(MD))
% MD.addPackage(GenericPackage(MD))
% MD
% h = MD.packages_{1}.GUI(MD);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Get Assaf''s Data');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

allCellsFileName = '28-Mar-2017_LBP_dLBP_1.mat';
InFilePath =  '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData';
loadMatPath = fullfile(InFilePath,allCellsFileName);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Load Assaf''s Data');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

load(loadMatPath, 'allCellsMovieData'); 

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Select (Random) Subset');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

randset = 50;
cellDataSubSet = allCellsMovieData(1,randsample(length(allCellsMovieData),randset,false));
    
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Load Assaf''s MDs');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% load movies into memory

nCells = numel(cellDataSubSet);
cellMovies = cell(1,nCells);

parfor i = 1:nCells

    iCell = cellDataSubSet{i};
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
          
    for itime = 1:ntime
        curTime = ts(itime);
        curI = MD.getChannel(1).loadImage(curTime);

        movie(:,:,itime) = curI(round((ys(itime)-FOVRadius)):(round(ys(itime)+FOVRadius)),round((xs(itime)-FOVRadius)):round((xs(itime)+FOVRadius)));
    end
    cellMovies{i} = movie;
end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('% Save prepared data set to .mat file?');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% save(['/work/bioinformatics/s170480/Data/LCH/DevOps/testMoviesRand_' num2str(randset)]);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Run cellXplorere!');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

disp('Here are a few examples...')

disp('%%%%%% Load feature data without movies')4
disp('% cellXploreDR(cellDataSubSet);')
% cellXploreDR(cellDataSubSet);

disp('%%%%%% Load feature data and movies (but not the pre-defined annotation set)')
disp('% cellXploreDR(cellDataSubSet, cellMovies);')
% cellXploreDR(cellDataSubSet, cellMovies);

disp('%%%%%% Included your predefined annotation set (keys only)');
disp('cellXploreDR(cellDataSubSet, cellMovies, ''annotations'', cellDataSubSet{1}.annotations)');
% cellXploreDR(cellDataSubSet, cellMovies, 'annotations', cellDataSubSet{1}.annotations);

disp('%%%%%% Load in a pre-computed DR view')
disp('% First prepare data structure like so:')
disp('% DR.line = [1:50; y_DR]'' % (where DR.tSNE1 is a nxD matrix n==#cellMmovies, D==DR dims)');
disp('cellXploreDR(cellDataSubSet, cellMovies, ''annotations'', cellDataSubSet{1}.annotations, ''DR'', DR)');
DR.line = [1:50; 1:50]';
cellXploreDR(cellDataSubSet, cellMovies, 'annotations', cellDataSubSet{1}.annotations, 'DR', DR);


disp('Export saved annotations to the workpace by clicking the "export cellData with.."');
disp('Give the variable a name, multiple variables will be created in the workspace')
disp(' 1) the cell array with annotations included')
disp(' 2) the hashmap with tag to cellkey')
disp(' you can also export the internal data stucutre by clicking: celldata export');
disp('then you can take a look at the annotation info by the struct: cellXdata.meta.anno')

disp('The you can re-load the cell array')
disp(' e.g., cellXploreDR(cellXData_12Apr2017_0110, cellMovies )');

% %% how to close out figures -- 
disp('Run the following to close all figures')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('ff = findall(0,''Tag'',''cellXplore''); delete(ff);');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

