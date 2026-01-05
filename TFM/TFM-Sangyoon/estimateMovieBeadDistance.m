function [] = estimateMovieBeadDistance(ML)
%function [] = estimateMovieBeadDistance(ML) reads ML and MDs, pulls up one
%bead image, then run estimateBeadDensity and quantify the statistics about
%bead distance and density
%   input:
%       ML:         MovieList
%   output:
%       Graphs and data are stored where ML is located.
% Sangyoon Han. May 2020

%% Output set up
folderPath = [ML.getPath filesep 'BeadDistance_' ML.getFilename];
folderPath = folderPath(1:end-4);

%% Read ML
% ML=MovieList.load(ML.getFullPath);
numMovies = numel(ML.movieDataFile_);
%% Per each movieData, read a bead image and run the function
beadDistAllCell=cell(numMovies,1);
beadDensityAll=zeros(numMovies,1);

for ii=1:numMovies
    %Read
    curMD=MovieData.load(ML.movieDataFile_{ii});
    %See if TFM has run
    iTFM = curMD.getPackageIndex('TFMPackage');
    if isempty(iTFM)
        disp(['No TFM was done in this movie: ' curMD.getFilename '.'])
        disp('Using the first channel for bead image');
        %continue
        iBeadChan = 1;
    else
        tPack = curMD.getPackage(iTFM);
        % Get the channel
        dispCalProc = tPack.getProcess(2);
        iBeadChan = dispCalProc.checkChannelOutput;
    end
    beadChan = curMD.getChannel(iBeadChan);
    [beadDistAll,beadDensity] = estimateBeadDistance(beadChan.loadImage(1),curMD.pixelSize_,beadChan.psfSigma_);
    % save the results
    beadDistAllCell{ii}=beadDistAll;
    beadDensityAll(ii)=beadDensity;    
end
%% Plotting
% bead distance, then bead density
%% bead distance
if ~exist(folderPath,'dir')
    mkdir(folderPath)
end
beadDistGroup = cell2mat(beadDistAllCell);
h1=figure(1); hHist=histogram(beadDistGroup);
xlabel('Distance (\mum)')
ylabel('Number of occurrence')
[maxCount,iMaxCount]=max(hHist.BinCounts);
text(hHist.BinEdges(iMaxCount+1)+0.5*hHist.BinWidth,maxCount*0.95,['mean: ' num2str(mean(beadDistGroup),2) ' \mum']);
text(hHist.BinEdges(iMaxCount+1)+0.5*hHist.BinWidth,maxCount*0.85,['std: ' num2str(std(beadDistGroup),2) ' \mum']);

hgexport(h1,strcat(folderPath,'/beadDist'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(folderPath,'/beadDist'),'-v7.3')
print(h1,strcat(folderPath,'/beadDist.tif'),'-dtiff')

tableBeadDistGroup=table(beadDistGroup);
writetable(tableBeadDistGroup,strcat(folderPath,'/beadDist.csv'))

%% bead density
beadDensityAllCell{1}=beadDensityAll;
h2=figure(2); boxPlotCellArray(beadDensityAllCell,{'All movies'},1,0,1);
ylabel('Bead density (#/\mum^2)')

hgexport(h2,strcat(folderPath,'/beadDensity'),hgexport('factorystyle'),'Format','eps')
hgsave(h2,strcat(folderPath,'/beadDensity'),'-v7.3')
print(h2,strcat(folderPath,'/beadDensity.tif'),'-dtiff')

tableBeadDensityGroup=table(beadDensityAllCell);
writetable(tableBeadDensityGroup,strcat(folderPath,'/beadDensity.csv'))
