function allCellMovies = extractMoviesToCellExplorer(allCellsFname)

addpath(genpath('/home2/azaritsky/code/common'));
addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition/'));
addpath(genpath('/home2/azaritsky/code/extern'));

if nargin < 1    
    allCellsFname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/CellExplorerData/28-Mar-2017_LBP_dLBP_1.mat';
end

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

allCellsFnameOut = [allCellsFname(1:end-4) '_movies.mat'];

%% HARD CODED!
pixelSize = 0.325;
FOVRadius = round(35/pixelSize);
% curFOVRadius = round(params.FOVRadius*curScale);
%%

load(allCellsFname); % allCellsMovieData
nCells = length(allCellsMovieData);


allCellMovies = cell(1,nCells);

for icell = 1 :nCells
    curCellData = allCellsMovieData{icell};
    
    xs = curCellData.xs;
    ys = curCellData.ys;
    ts = curCellData.ts;
    
    ntime = length(ts);
    
    MD =  MovieData.load(curCellData.MD);
    %     [sizeY,sizeX] = size(MD.getChannel(1).loadImage(ts(1)));
    
    movie = zeros(2*FOVRadius+1,2*FOVRadius+1,ntime);
    
    for itime = 1 : ntime
        curTime = ts(itime);
        curI = MD.getChannel(1).loadImage(curTime);
        % For feature extraction - try to keep the BB constant! (here it is dynamic..)
        movie(:,:,itime) = curI((ys(itime)-FOVRadius):(ys(itime)+FOVRadius),(xs(itime)-FOVRadius):(xs(itime)+FOVRadius));
    end
    allCellMovies{icell} = movie;
end

save(allCellsFnameOut,'allCellMovies');

end


