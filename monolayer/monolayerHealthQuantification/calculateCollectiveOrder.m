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
    vg = gradient(vrkbttemp); %gradient of vrkbt to form sine wave
    windowSize = round(length(vg) / 50); %calculate window size from signal length
    if windowSize < 6 %verify window size is not too small
        windowSize = length(vg);
    end
    movRMS = sqrt(movmean(vg.^2, windowSize)); %moving rms to find signal amplitude over time
    halfpt = round(length(movRMS) * 0.75); %get idx 3 quarters of the way into the signal for truncating later
    movRMS(isinf(movRMS)) = nan; %replace infinite values with nans
    movRMS = fillmissing(movRMS,'nearest'); %replace nans with nearest neighbor
    [kneeIDX,~] = knee_pt(movRMS(1:halfpt)); %calculate knee point using truncated signal
    if length(vrkbttemp) > (kneeIDX + 10) %verify vrkbt is long enough before adding 5 to index to push it further into the plateau region
        kneeIDX = kneeIDX + 5;
    end
    if length(vrkbttemp) > (kneeIDX + 10) && kneeIDX > 10 %ensure that plateau value averageing window does not exceed array boundaries
        platVal = mean(vrkbttemp(kneeIDX-10:kneeIDX+10));
    elseif kneeIDX < 10 && length(vrkbttemp) > (2*kneeIDX + 1)
        platVal = mean(vrkbttemp(1:(2*kneeIDX)));
    else
        platVal = vrkbttemp(kneeIDX);
    end
    vrkbttemp = vrkbttemp - platVal; %scale the curve such that the plateau is located around zero
    [~,indMin] = min(vrkbttemp(1:kneeIDX)); %find location of first local minima (aka global minima since it is almost guaranteed to be the same)
    padSize = min(10,indMin); %adaptive padding size to ensure window does not exceed array boundaries
    leftPad = padSize;
    rightPad = padSize;
    for aa = 1:padSize %determine edges of quadratic trough for fitting
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

function [res_x, idx_of_result] = knee_pt(y,x,just_return)
%function [res_x, idx_of_result] = knee_pt(y,x,just_return)
%Returns the x-location of a (single) knee of curve y=f(x)
%  (this is useful for e.g. figuring out where the eigenvalues peter out)
%
%Also returns the index of the x-coordinate at the knee
%
%Parameters:
% y (required) vector (>=3 elements)
% x (optional) vector of the same size as y
% just_return (optional) boolean
%
%If just_return is True, the function will not error out and simply return a Nan on
%detected error conditions
%
%Important:  The x and y  don't need to be sorted, they just have to
%correspond: knee_pt([1,2,3],[3,4,5]) = knee_pt([3,1,2],[5,3,4])
%
%Important: Because of the way the function operates y must be at least 3
%elements long and the function will never return either the first or the
%last point as the answer.

%set internal operation flags
use_absolute_dev_p = true;  %ow quadratic
%deal with issuing or not not issuing errors
issue_errors_p = true;
if (nargin > 2 && ~isempty(just_return) && just_return)
    issue_errors_p = false;
end
%default answers
res_x = nan;
idx_of_result = nan;
%check...
if (isempty(y))
    if (issue_errors_p)
        error('knee_pt: y can not be an empty vector');
    end
    return;
end
%another check
if (sum(size(y)==1)~=1)
    if (issue_errors_p)
        error('knee_pt: y must be a vector');
    end
    
    return;
end
%make a vector
y = y(:);
%make or read x
if (nargin < 2 || isempty(x))
    x = (1:length(y))';
else
    x = x(:);
end
%more checking
if (ndims(x)~= ndims(y) || ~all(size(x) == size(y)))
    if (issue_errors_p)
        error('knee_pt: y and x must have the same dimensions');
    end
    
    return;
end
%and more checking
if (length(y) < 3)
    if (issue_errors_p)
        error('knee_pt: y must be at least 3 elements long');
    end
    return;
end
%make sure the x and y are sorted in increasing X-order
if (nargin > 1 && any(diff(x)<0))
    [~,idx]=sort(x);
    y = y(idx);
    x = x(idx);
else
    idx = 1:length(x);
end
%the code below "unwraps" the repeated regress(y,x) calls.  It's
%significantly faster than the former for longer y's
%
%figure out the m and b (in the y=mx+b sense) for the "left-of-knee"
sigma_xy = cumsum(x.*y);
sigma_x  = cumsum(x);
sigma_y  = cumsum(y);
sigma_xx = cumsum(x.*x);
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mfwd = (n.*sigma_xy-sigma_x.*sigma_y)./det;
bfwd = -(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det;
%figure out the m and b (in the y=mx+b sense) for the "right-of-knee"
sigma_xy = cumsum(x(end:-1:1).*y(end:-1:1));
sigma_x  = cumsum(x(end:-1:1));
sigma_y  = cumsum(y(end:-1:1));
sigma_xx = cumsum(x(end:-1:1).*x(end:-1:1));
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mbck = flipud((n.*sigma_xy-sigma_x.*sigma_y)./det);
bbck = flipud(-(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det);
%figure out the sum of per-point errors for left- and right- of-knee fits
error_curve = nan(size(y));
for breakpt = 2:length(y)-1
    delsfwd = (mfwd(breakpt).*x(1:breakpt)+bfwd(breakpt))-y(1:breakpt);
    delsbck = (mbck(breakpt).*x(breakpt:end)+bbck(breakpt))-y(breakpt:end);
    %disp([sum(abs(delsfwd))/length(delsfwd), sum(abs(delsbck))/length(delsbck)])
    if (use_absolute_dev_p)
        % error_curve(breakpt) = sum(abs(delsfwd))/sqrt(length(delsfwd)) + sum(abs(delsbck))/sqrt(length(delsbck));
        error_curve(breakpt) = sum(abs(delsfwd))+ sum(abs(delsbck));
    else
        error_curve(breakpt) = sqrt(sum(delsfwd.*delsfwd)) + sqrt(sum(delsbck.*delsbck));
    end
end
%find location of the min of the error curve
[~,loc] = min(error_curve);
res_x = x(loc);
idx_of_result = idx(loc);
end