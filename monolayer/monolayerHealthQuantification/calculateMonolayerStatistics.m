function [SD,stats,centroids] = calculateMonolayerStatistics(filePath,binList,magList,nameList,varargin)
%CALCULATEMONOLAYERSTATISTICS calculates cellular statistics including 
%   circularity, extent, area, and density using the segmented images saved
%   as labels.mat from the segmentCellMonolayer function. Additional inputs 
%   are used to group the segmented images by condition, as well as provide
%   names for each condition.
%
%INPUTS:
%   filePath:   Filepath of the folder containing the saved labels.mat
%       file, this will be identical to the filePath input into
%       segmentCellMonolayer function
%
%   binList:    Integer array of bin numbers for each image being analyzed.
%       For example; [1, 1, 1, 2, 2, 2, 2] would group the first three
%       images into one bin/condition and the next four images into a
%       second bin/condition. Indices in binList must equal the length of
%       labels.mat (aka total number of images)
%
%   magList:    Integer array of the magnification amount used to take the
%       microscopic images. For example; [10, 10, 10, 4, 4, 4, 4] would
%       establish that 10x magnification was used to take the first three
%       images and 4x magnification was used to take the last four. The
%       magnification amounts do not need to be consistent across or within
%       bins.
%
%   nameList:   Cell array which provides the desired names as strings for 
%       each bin/condition to be used when plotting the results. The total
%       number of names must be equal to the total number of different
%       bins. For example; if the binList is [1, 1, 1, 2, 2, 2, 2] then two
%       names must be given in the nameList such as {'Pat. 47','Pat. 109'}
%
%   doPlot:     Optional binary input which defines if the results should 
%       be plotted. Zero would mean no plotting is done. (default = 1)
%
%OUTPUTS:
%   SD:         Segmentation data structure that contains the stats and
%       centroids variables, the structure is also populated further in
%       downstream functions and serves to transfer data between functions.
%
%   stats:      Cell array of same length as labels which contains all the
%       calculated stats in a structure array for each image.
%
%   centroids:  Cell array of same length as labels which contains lists of
%       all cellular centroids for each image. This variable is used for
%       the collective order calculation, however, it is also saved to the
%       folder specified in filePath for use downstream. Additionally,
%       concatenated at the end of this variable is the image size, be
%       mindful of this fact when using the centroids variable outside of
%       the function processes.
%
%EXAMPLE:
%   [SD,stats,centroids] = calculateMonolayerStatistics('/storage/disk1/sehaarma/GP065 SET',[1,1,2,2],[10,10,10,10],{'Con. 1','Con. 2'})
%       Analyzing the four images generated from the previous example 
%       in segmentCellMonolayer, binning the first two images into bin 1 
%       and the second two into bin 2. Additionally specifying a 
%       magnification amount of 10x for all iamges and labeling the two 
%       bins Con. 1 and Con. 2 respectively.

%% PARSE INPUTS
disp('Running calculateMonolayerStatistics')
fprintf('Parsing Inputs...');
p = inputParser;
defaultDoPlot = 1;
validBinary = @(x) x == 1 || x == 0;
validBinList = @(x) length(x) == length(magList);
validMagList = @(x) length(x) == length(binList);
validString = @(x) iscell(x);
addRequired(p,'filePath',@(x)exist(x,'dir'));
addRequired(p,'binList',validBinList);
addRequired(p,'magList',validMagList);
addRequired(p,'nameList',validString);
addOptional(p,'doPlot',defaultDoPlot,validBinary);
parse(p,filePath,binList,magList,nameList,varargin{:});
doPlot = p.Results.doPlot;
disp(' done!')

%% ** LOAD SD STRUCTURE **
SD = load(strcat(filePath,filesep,'SD.mat')); %load cell array of image labels
imSize = cellfun(@size,SD.labels,'UniformOutput',false);
SD.imSize = imSize;

%% ** CALCULATE STATISTICS **
fprintf('Calculating statistics...')
props = {'Area','Centroid','Circularity','Eccentricity','Extent',...
    'MajorAxisLength','MinorAxisLength','Orientation'}; %define properties to analyze

f = @(a) regionprops(a,props); %define anonymous function for regionprops with defined props
stats = cellfun(f,SD.labels,'UniformOutput',false); %call anonymous function in cellfun to run regionprops on each cell in labels

centroids{length(stats)} = []; %initialize centroids variable
for ii = 1:length(stats) %loop through each image
    centList = [stats{ii}.Centroid]; centList = centList'; %get list of centroids and convert to vertical array
    centroids{ii} = [centList(2:2:end),centList(1:2:end)]; %convert centroid array to Nx2 vector list
end
SD.centroids = centroids;

%calculate scaling based on image magnification, used for density calculations
scaleFactor = max(magList)./magList;

%Calculating average area, circularity, extent, and density for each image 
%and binning them according to the bin number specified in binList and
%associated name in nameList
totalBins = max(binList); %get total number of bins
cellArrayArea{totalBins} = []; %initialize area array
cellArrayCirc{totalBins} = []; %initialize circ array
cellArrayExte{totalBins} = []; %initialize exte array
cellArrayDens{totalBins} = []; %initialize dens array
for jj = 1:totalBins %loop through each bin
    indx = binList == jj; %get logical array with ones corresponding to the images assigned to the current bin
    meanArea = zeros(1,sum(indx)); %initialize meanarea
    meanCirc = zeros(1,sum(indx)); %initialize meancirc
    meanExte = zeros(1,sum(indx)); %initialize meanexte
    cellDens = zeros(1,sum(indx)); %initialize celldens
    ll = 1;
    for kk = find(indx,1):find(indx,1,'last') %loop through the images in the current bin (jj)
        meanArea(ll) = mean([stats{kk}.Area]); %calculate mean area for current image (kk)
        meanCirc(ll) = mean([stats{kk}.Circularity]); %calculate mean circularity for current image (kk)
        meanExte(ll) = mean([stats{kk}.Extent]); %calculate mean extrent for current image (kk)
        numCells = length([stats{kk}.Area]); %calculate total number of cells in current image (kk)
        %calculate cell density in current image (kk)
        cellDens(ll) = numCells / ((imSize{kk}(1) * scaleFactor(kk)) * (imSize{kk}(2) * scaleFactor(kk))); 
        ll = ll + 1;
    end
    cellArrayArea{jj} = meanArea; %store mean area in cell array for plotting
    cellArrayCirc{jj} = meanCirc; %store mean circ in cell array for plotting
    cellArrayExte{jj} = meanExte; %store mean exte in cell array for plotting
    cellArrayDens{jj} = cellDens; %store cell dens in cell array for plotting
end
disp(' done!')

%% ** PLOT RESULTS **
if doPlot
    %Area
    figure,
    boxPlotCellArray(cellArrayArea,nameList,1,0,0,'forceShowP',1,'forceTtest',1);
    title('Cell Area')
    %Circularity
    figure,
    boxPlotCellArray(cellArrayCirc,nameList,1,0,0,'forceShowP',1,'forceTtest',1);
    title('Cell Circularity')
    %Extent
    figure,
    boxPlotCellArray(cellArrayExte,nameList,1,0,0,'forceShowP',1,'forceTtest',1);
    title('Cell Extent')
    %Density
    figure,
    boxPlotCellArray(cellArrayDens,nameList,1,0,0,'forceShowP',1,'forceTtest',1);
    title('Cell Density')
end

%% ** SAVE RESULTS **
%Saving stats
SD.Area = cellArrayArea;
SD.Circularity = cellArrayCirc;
SD.Extent = cellArrayExte;
SD.Density = cellArrayDens;
%Saving bin information
SD.binList = binList;
SD.magList = magList;
SD.nameList = nameList;
save(strcat(filePath,filesep,'SD.mat'),'-struct','SD')

end

%add input checker which verifies that labels and binList have the same
%length
%
%add input checker that verifies the max number in bin list is equal to the
%length of nameList
%
%add input checker that verifies that binlist and maglist ahve the same
%length
%
%want an input variable that bins each of the images analyzed
%so if the users wants images 1-4 in bin 1 and then 5-7 in bin 2
%the input would look something like [1, 1, 1, 1, 2, 2, 2,];
%then i need to program binning of the stats into arrays which are then
%placed into data cell arrays
%
%also need an input which defines the names of the bins which will then be
%placed into a label cell array