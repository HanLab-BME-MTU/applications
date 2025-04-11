function [SD,grtemp,vrkbttemp,k] = calculateCollectiveOrder(filePath,dr,border,varargin)
%CALCULATECOLLECTIVEORDER takes the segmentation data structure saved to
%   SD.mat and calculates the collective order for each image in SD using
%   radius band width of dr and image boundary culling set by bounds.
%
%INPUTS:
%   filePath:   The folder path to the SD.mat file which contains the
%       necessary segmentation data structure.
%
%   dr:         The width of the radius band used in the calculation of the
%       radial distribution function. (default = 2)
%
%   border:     Number of pixels to ignore at each edge of the image. This
%       help reduce error as the edge centroids have poor radial 
%       distribution due to lacking some neighbors in the direction of the 
%       edge. (default = 20)
%
%   doPlot:     Optional binary input which defines if the results should 
%       be plotted. Zero would mean no plotting is done. (default = 1)
%
%OUTPUTS:
%   SD:         Segmentation data structure which contains various
%       variables necessary for the various functions in the segmeentation 
%       and analysis process.
%
%   gr:         Radial distribution function used to calculate vrkbt. This
%       variable is also saved to the SD structure.
%
%   vrkbt:      Normalized potential of mean force function used to
%       calculate the collective order parameter k. This variable is also
%       saved to the SD structure.
%   
%   k:          Collective order parameter stored in a cell array of size
%       1xNumImages. This variable is also saved to the SD structure.

%% ** LOAD SD STRUCTURE **
SD = load(strcat(filePath,filesep,'SD.mat')); %load segmentation data structure
minSize = min([SD.imSize{:}]); %quickly determine the smallest image dimension for input parsing purposes

%% PARSE INPUTS
disp('Running calculateCollectiveOrder')
fprintf('Parsing Inputs...');
p = inputParser;
defaultDoPlot = 1;
validBinary = @(x) x == 1 || x == 0;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosNumNonMax = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= minSize);
addRequired(p,'filePath',@(x)exist(x,'dir'));
addRequired(p,'dr',validScalarPosNum);
addRequired(p,'border',validScalarPosNumNonMax);
addOptional(p,'doPlot',defaultDoPlot,validBinary);
parse(p,filePath,dr,border,varargin{:});
doPlot = p.Results.doPlot;
disp(' done!')

%% ** CALCULATE CENTROID DISTANCES **
fprintf('Fitting Curves...')
centroids = SD.centroids; %separate out centroid list
imSize = SD.imSize; %grab image size
totalImages = length(imSize); %count total images
densityList = [SD.Density{:}]; %separate and flatten cell densities

ll = 1; %initialize loop parameter
gr{totalImages} = []; %initialize radial distribution function
vrkbt{totalImages} = []; %initialize normalized potential mean force
k{totalImages} = []; %initialize collective order parameter
for ii = 1:totalImages
    for jj = 1:length(centroids{ii}) %loop through all centroids
        x = centroids{ii}(jj,:); %grab current centroid from list
        %verify that the current centroid is not too close the image edge
        if x(1) > imSize{ii}(1)-border || x(1) < border
            if x(2) > imSize{ii}(2)-border || x(2) < border
                continue
            end
        end
        temp = centroids{ii}; %grab centroid list
        temp(jj,:) = []; %delete current centroid to avoid self-comparison
        for jj = 1:length(temp) %loop through all other centroids
            y = temp(jj,:);
            currImDists(ll) = norm(x - y); %distance between current centroid (x) and all other centroids (y)
            ll = ll + 1; %increment loop parameter
        end
    end

%% ** CALCULATE G(R) and V(R) **
    currSortedDists = sort(currImDists); %sort data by increasing distance between centroids
    currEdges = 0:dr:(currSortedDists(end)+dr); %generate bin edges using dr as the bin width
    [N,~] = histcounts(currSortedDists,currEdges); %binning distance data into dr sized bins
    r = 0:dr:currSortedDists(end); %generate radius array as independent variable
    shellArea = 2*pi*r*dr; %calculate shell area over increasing radii
    grtemp = N ./ (shellArea .* densityList(ii)); %normalize shell area to obtain g(r)
    indx = find(grtemp(1:round(length(grtemp)/2))==0,1,'last');
    grtemp = grtemp(indx:end);
    vrkbttemp = -1 .* log(grtemp); %inverted natrual log to obtain v(r)kbT
    
%% ** FIT V(R) and CALCULATE K **
    vx = 0.1:0.1:length(vrkbttemp)/10;
    platVal = mean(vrkbttemp(75:100));
    vrkbttemp = vrkbttemp - platVal;
    localMin = islocalmin(vrkbttemp,'MinSeparation',1000); %find minimum point 
    indMin=find(localMin,1); %get index of minimum
    leftPad = 10;
    rightPad = 10;
    for aa = 1:10
      if vrkbttemp(indMin(1)-leftPad) > 0
         leftPad = leftPad - 1;
      end
      if vrkbttemp(indMin(1)+rightPad) > 0
         rightPad = rightPad - 1; 
      end
    end
    xfit=vx((indMin(1)-leftPad):(indMin(1)+rightPad)); %extract x values for fitting
    vfit=vrkbttemp((indMin(1)-leftPad):(indMin(1)+rightPad)); %extract y values for fitting

    f=polyfit(xfit,vfit,2); %fit the quadratic
    fittedCurve = polyval(f,xfit); %generate discrete values from fit
    ktemp = abs(2 * f(1)); %calculate collective order
    
%% ** PLOT RESULTS **    
    if doPlot
        figure,
        plot(vx,vrkbttemp,'LineWidth',2, 'Color',"#D95319"), hold on
        plot(xfit,fittedCurve,'LineWidth',3, 'Color',"[0 0 0]"), hold on
        plot(vx(1:1000),zeros(1,1000),'k--'), hold on
        title('V(r) with Fitted Curve','FontSize',16)
        xlabel('Indx','FontSize',16)
        ylabel('V(r)','FontSize',16)
        xlim([0 10])
        text(vx(indMin)+1,vrkbttemp(indMin)+1,strcat('k = ',num2str(ktemp))), hold off
    end
    
    gr{ii} = grtemp;
    vrkbt{ii} = vrkbttemp;
    k{ii} = ktemp;
end
disp(' done!')

%% ** BOX PLOT RESULTS **
% %close all
% if doPlot
%     binList = SD.binList;
%     totalBins = max(binList);
%     cellArrayColl{totalBins} = [];
%     for ii = 1:totalBins
%        indx = binList == ii;
%        collOrder = zeros(1,sum(indx));
%        ll = 1;
%        for jj = find(indx,1):find(indx,1,'last')
%            collOrder(ll) = k{jj};
%            ll = ll + 1;
%        end
%        cellArrayColl{ii} = collOrder;
%     end
%     SD.cellArrayColl = cellArrayColl;
%     figure,
%     boxPlotCellArray(cellArrayColl,SD.nameList,1,0,0,'forceShowP',1);
%     title('Collective Order')
% end

%% ** SAVE RESULTS ** 
SD.gr = grtemp;
SD.vrkbt = vrkbttemp;
SD.k = k;
save(strcat(filePath,filesep,'SD.mat'),'-struct','SD')

end

