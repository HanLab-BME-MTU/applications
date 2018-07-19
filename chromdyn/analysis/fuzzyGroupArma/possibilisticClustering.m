function [centers, dataWeight, distances, eta, typicality, membership, diagnostics] = possibilisticClustering(data,initialGuess,centerFunction,distanceFunction,eta, pcOptions, centerFunctionParameters, distanceFunctionParameters)
%POSSIBILISTICCLUSTERING is a fuzzy clustering variant describing data points in relation to cluster centers
%
% SYNOPSIS: [centers, dataWeight, distances, eta, typicality, membership] = ...
%               possibilisticClustering(data,initialGuess,...
%               centerFunction,distanceFunction,eta,m, ...
%               centerFunctionParameters, distanceFunctionParameters)
%
%               possibilisticClustering(data,initialGuess,...
%               centerFunction,distanceFunction,eta,pcOptions, ...)
%
% INPUT data :          data to be clustered. Needs to be compatible with
%                       centerFunction, distanceFunction. size(data,1) has
%                       to return a scalar indicating how many data points
%                       there are
%		initialGuess:   1-by-n struct array with initial guess for the
%                       centers, or nData-by-1 array of group labels or
%                       nData-by-nCenters array of memberships.
%                       In order to allow the user to manually choose
%                       initial guesses, supply a cell array:
%                       {dm, referenceCoordinates, doFCM}
%                       dm: dissimilarity matrix between the data points.
%                       referenceCoordinates: (opt) nData-by-2 array of
%                          reference coordinates to orient the MDS plot.
%                       doFCM: (opt),{0}/1 whether or not to run FCM on the
%                          hard partition from the ROIs. Alternatively, the
%                          initial guesses for the centers are calculated
%                          using only the points within the ROIs.
%                       This approach opens a figure for the user to select
%                       ROIs to specify potential clusters.
%		centerFunction: function to calculate cluster prototypes from the
%                       data. CenterFunction should take four input
%                       arguments:
%                       data, membership, typicality, m,
%                       centerFunctionParameters.
%                       centerFunctionParameters are arbitrary user defined
%                       values. For information on membership and
%                       typicality, see below. The funcion should be able
%                       to handle typicalities of all 1
%                       As output, it should return a 1-by-nCenters struct
%                       array with the characteristics of the cluster
%                       prototypes
%		distanceFunction: function to calculate the distance between the
%                       individual data points and the cluster centers.
%                       CenterFunction should take three input arguments:
%                       data, centers, distanceFunctionParameters.
%                       distanceFunctionParameters are arbitrary user
%                       defined values.
%                       As output, it should return a nData-by-nCenters
%                       array of distances, where an element (i,j) contains
%                       the distance of data point i to the cluster
%                       prototype j.
%                       IF YOU USE 2-NORMS, SQUARE THE DISTANCE ALREADY!
%       eta, m          (opt) Clustering parameters. eta is about equal
%                       to the square of the diameter of the cluster (it is
%                       where the cluster membership becomes 0.5. If not
%                       supplied, it will be estimated. eta can either be
%                       scalar or contain different values for each center
%                       m is the fuzzyness-parameter (has to be >1). Values
%                       close to 1 result in hard clustering (distances
%                       above eta have weight 0, below eta, they have
%                       weight 1). The higher m, the fuzzier the clusters.
%                       Default: 1.5
%       pcOptions       (opt) In place of m, you can submit an options
%                       structure with the following optional fields
%                       .m - m, as above
%                       .algorithm - 1: original Krishnapuram&Keller 1993,
%                                    2: Zhang & Leung 2004
%                                    3: Dorn 2007
%                       .debug - 1: stop and calculate objective function
%                                   inside the main clustering loop
%       centerFunctionParameters, distanceFunctionParameters: optional
%                       parameters for distance function. Note: Even if
%                       your distance function doesn't need additional
%                       parameters, it has to be able to accept them as
%                       (empty) inputs.
%
% OUTPUT centers :      1-by-nCenters struct array with the characteristics
%                       of the final cluster prototypes as determined by
%                       centerFunction.
%        dataWeight:    nData-by-nCenters array of weights, where
%                       element (i,j) describes the weight of data point i
%                       in the calculation of center j.
%                       dataWeight = typicality^m*membership^m/sum
%        distances :    nData-by-nCenters array as output by
%                       distanceFunction
%        eta :          1-by-nCenters array of final estimates for eta.
%        typicality:    nData-by-nCenters array of typicalities
%                       (possibilistic clustering)
%        membership:    nData-by-nCenters array of memberships
%                       (probabilistic clustering)
%
% REMARKS  fcmCenterFunction and fcmDistanceFunction (subfunctions of this
%          function) can be used to perform possibilistic c-means
%          clustering. P-FCM assumes spherical clusters.
%          Centers are characterized by .mean
%          gkCenterFunction and gkDistanceFunction (subfunctions of this
%          function) are implementations of the Gustafsson-Kessel algorithm
%          that allows arbitrarily oriented elliptical shapes.
%          Centers are characterized by .mean and .covariance (both fields
%          are required in the initial guess!).
%          Both algorithms take [x y ...] coordinate lists as input.
%
%
%          For 2-d data, plot results with possibilisticClusteringPlot
%
%          see Krishnapuram & Keller (1993) A possibilistic approach to
%             clustering. IEEE Transactions on Fuzzy Systems 1:98
%
%             Fuzzy cluster analysis: methods for classification, data
%             analysis and image recognition. By: Hoeppner, Kruse, Klawonn
%             and Runkler. Wiley 1999
%
%             Zhang & Leung (2004) Improved Possibilistic C-Means
%             Clustering Algorithms. IEEE Transactions on Fuzzy Systems
%             12:2
%
%             The code is different from the above cited works in that it
%             applies a density correction to the cluster typicalities.
%             This density correction will avoid coincident cluster
%             prototypes in the case of overlapping clusters
%
%          It may seem a bit annoying to always have to calculate the
%          mth power of the membership, but that's how the algorithms
%          are usually constructed.
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 21-Oct-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default options
def_algorithm = 3;
def_m = 1.5;
def_debug = 0;
diagnostics = struct;


%==========================
%% TEST INPUT
%==========================

% we require 4 input arguments (including the inital guess)
if nargin < 4 || isempty(data) || isempty(initialGuess) || ...
        isempty(centerFunction) || isempty(distanceFunction)
    error('POSSIBILISTICCLUSTERING:INPUTERROR','possibilisticClustering requires 4 nonempty input arguments')
end

% we don't really care about data - that's for the user to determine.
% However, we need to be able to get the length
nData = size(data,1);
if nData == 1 && size(data,2)>1
    data = data';
    nData=size(data,1);
end

% check whether distanceFunction, centerFunction exist
if isempty(which(distanceFunction))
    error('POSSIBILISTICCLUSTERING:INPUTERROR','distanceFunction %s not found',distanceFunction)
end
if isempty(which(centerFunction))
    error('POSSIBILISTICCLUSTERING:INPUTERROR','centerFunction %s not found',centerFunction)
end
% check for parameters
if nargin < 7
    centerFunctionParameters = [];
end
if nargin < 8
    distanceFunctionParameters = [];
end


if nargin < 6 || isempty(pcOptions)
    m = def_m;
    algorithm = def_algorithm;
    debug = def_debug;
else
    if isstruct(pcOptions)
        % pcOptions structure
        if isfield(pcOptions,'m')
            m = pcOptions.m;
        else
            m = def_m;
        end
        if m<=1
            warning('POSSIBILISTICCLUSTERING:ADJUSTM',...
                'm was below 1 (%1.5f). It was set to 1+eps',m)
            m = max(m,1+eps);
        end
        if isfield(pcOptions,'algorithm')
            algorithm = pcOptions.algorithm;
        else
            algorithm = def_algorithm;
        end
        if isfield(pcOptions,'debug')
            debug = pcOptions.debug;
        else
            debug = def_debug;
        end

    else
        % m, not pcOptions is given
        m = pcOptions;
        % make sure m is above 1
        if m<=1
            warning('POSSIBILISTICCLUSTERING:ADJUSTM',...
                'm was below 1 (%1.5f). It was set to 1+eps',m)
            m = max(m,1+eps);
        end
        algorithm = def_algorithm;
        debug = def_debug;
    end
end

% check for proper initialGuess - do that at the end, because of possible
% center calculation
% init = 1: centers given
% init = 2: fuzzy memberships given
% init = 3: hard memberships given
% init = 4: manual initial guess
init = 0;
if isstruct(initialGuess)
    % it's a center array
    init = 1;
    centers = returnRightVector(initialGuess,1,'r');
    nCenters = length(centers);
    membership = 0;

elseif isnumeric(initialGuess)
    if any(size(initialGuess)==1)
        % it's a label array - convert to membership

        % find unique entries
        [dummy,dummy,colIdx] = unique(initialGuess);
        % convert to indices
        nCenters = max(colIdx);
        idxList = sub2ind([nData,nCenters],(1:nData)',colIdx);

        membership = zeros(nData,nCenters);
        membership(idxList) = 1;

        init = 3;

    else
        % it's already a membership array
        [nMemberships,nCenters] = size(initialGuess);
        if nMemberships == nData
            membership = initialGuess;
        else
            error('POSSIBILISTICCLUSTERING:INPUTERROR','every data point needs a membership!')
        end

        % if there are not more than 3 different membership levels,
        % consider it a hard membership
        if length(unique(membership(:))) < 4
            init = 3;
        else
            init = 2;
        end
    end

elseif iscell(initialGuess)
    % manual initial guess. Deal with it here, because we need nCenters
    [init, membership, centers] = ...
        manualInitialGuess(initialGuess,data,m,...
        centerFunction,centerFunctionParameters);
    % init has all the relevant information for below, we just need to
    % determine nCenters
    nCenters = size(membership,2);


else
    error('POSSIBILISTICCLUSTERING:INPUTERROR',...
        'possibilisticClustering cannot understand initialGuess')
end


% typicality, eta need nCenters, which needs initial guess

% initialize typicality
typicality = ones(nData,nCenters);

% check eta
if nargin < 5 || isempty(eta)
    % estimate eta
    estimateEta = true;
    eta = 1;
elseif isscalar(eta)
    estimateEta = false;
    eta = repmat(eta,1,nCenters);
    % assume that the user can count both the number of centers and the
    % number of etas
else
    estimateEta = false;
end

% check for reestimateEta
if any(eta<0)
    reestimateEta = false;
    eta = abs(eta);
else
    reestimateEta = true;
end

w = warning;
warning('off','MATLAB:divideByZero')

%==========================================================================


%=======================
%% INIT LOOP
%=======================

% if we have to estimate eta, or if we have hard memberships, do three rounds
% of probabilistic clustering (could be extended to run every time, or on
% demand)
if estimateEta || init == 3

    % check whether we can calculate a center in the first round
    if init == 1
        getCenter = false;
    else
        getCenter = true;
    end

    for i=1:3

        % calculate centers
        if getCenter
            centers = feval(centerFunction, ...
                data, membership, typicality, m, centerFunctionParameters);
        else
            % if we come here the next time, membership is defined
            getCenter = true;
        end

        % calculate distances
        distances = feval(distanceFunction, ...
            data, centers, distanceFunctionParameters);

        % calculate "memberships"
        % don't use noise-clustering (maybe with the robustMean, we can
        % find outlier distances to determine the noise cutoff in the
        % future)
        % membership = 1/{sum(j=1:nCenters)((dxk^2/dxj^2)^(1/(m-1)))}
        % Theorem 1.11 in the book
        % Eq.1 in Krishnapuram
        [nData,nCenters] = size(distances);
        membership = zeros(nData,nCenters);
        for c = 1:nCenters
            membership(:,c) = 1./sum(...
                (repmat(distances(:,c),1,nCenters)./distances).^(1/(m-1)),2);
        end
        % werever the distance was 0, we get NaN. Replace with 1
        membership(isnan(membership)) = 1;
    end

    % estimate eta if necessary
    % Eq. 1.10 in the book
    % Eq. 9 in Krishnapuram
    if estimateEta
        eta = sum(membership.^m.*distances,1)./sum(membership.^m,1);
    end

    %possibilisticClusteringPlot(data,membership,centers),hold on, plot([0,1],[0,1],'.b')

    % update initialization status
    init = 4;

end

%=======================


%==========================
%% POSSIBILISTIC CLUSTERING
%==========================

% loop until the norm of membership doesn't change anymore
done = false;
if init == 1
    getCenter = false;
else
    getCenter = true;
end


% switch algorithm
%     case {1,2}
%         minIter = 2;
%     case 3
%         minIter = 11;
% end
minIter = 2;

% loop twice. Re-estimate eta using a membership cutoff of 0.25 (higher?)
doneOptimize= false;
doneRepeat = false;
iRepeat = 1;
while ~ doneRepeat
    iter = 0;
    while ~doneOptimize
        iter = iter + 1;

        % calculate centers
        if getCenter
            centers = feval(centerFunction, ...
                data, membership, typicality, m, centerFunctionParameters);



        else
            % if we come here the next time, membership is defined
            getCenter = true;
        end

        % calculate distances
        distances = feval(distanceFunction, ...
            data, centers, distanceFunctionParameters);
        if debug
            keyboard
            % objective function
            J = @(distances) sum(sum(membership.^m.*(typicality.^m.*distances + ...
                repmat(eta,nData,1) * (1 - typicality).^m)));

            for x = 1:21
                for y = 1:21
                    deltaX = (x-11)/100;
                    deltaY = (y-11)/100;
                    cTmp = centers(1).mean + [deltaX,deltaY];
                    distances = feval(distanceFunction, ...
                        data, cTmp, distanceFunctionParameters);
                    obj(x,y) = J(distances);
                end
            end

            figure,surf(obj)
            colormap(isomorphicColormap('green'));
        end

        % calculate typicality
        % typicality = 1/{1+((dxk^2/eta)^(1/(m-1)))}
        % Theorem 1.12 in the book
        % Eq.8 in Krishnapuram

        % remember old typicality
        oldTypicality = typicality;

        % If there are two overlapping clusters, the density of points in
        % the overlap region will make two close-by prototypes move on top
        % of each other. Therefore, correct the membership for overlapping
        % density.
        % see also Zhang and Leung

        typicality = ...
            1./(1+(distances./repmat(eta,nData,1)).^(1/(m-1)));
        membership = ones(nData,nCenters);


        switch algorithm
            case 1
                % no correction for typicality
            case 2
                % Zhang and Leung 2004
                for c = 1:nCenters
                    membership(:,c) = 1./sum(...
                        (repmat(typicality(:,c),1,nCenters)./typicality).*...
                        (repmat(distances(:,c),1,nCenters)./distances).^(1/(m-1)),2);
                end
            case 3
                % my version
                % to correct the typicalities for overlaps, we really should be using the
                % typicalities themselves, because they are more informative about the
                % cluster shape! This approach will in practice decrease the low
                % typicalities far beyond the cutoff eta by a bit, which is
                % not so bad, as it merely diminishes the already weak influence of
                % outliers.

                % try to get a bit more robust - do first optim w/o size
                % correction
                if iRepeat == 1
                    nMembers = 1;
                else
                    nMembers = repmat(sum(typicality>0.8,1),nData,1);
                end
                wTypicality = nMembers.*typicality;
                membership = wTypicality./repmat(sum(wTypicality,2),1,nCenters);
            otherwise
                error('unknown algorithm!')


        end





        %         cutMembership = membership;
        %         % multiplying by the relative membership will lead to long-range
        %         % repulsion. Therefore, only use relative membership up to 2*eta
        %         cutMembership(distances > 2* repmat(eta,nData,1)) = 0;
        %         densityEstimation = cutMembership ./ repmat(sum(cutMembership,2),1,nCenters);
        %         densityEstimation(isnan(densityEstimation)) = 1;
        %         membership = membership.* densityEstimation;

        %        membership = zeros(nData,nCenters);
        %                 for c=1:nCenters
        %         membership(:,c) = ...
        %             1./(1+(distances(:,c)./repmat(eta(:,c),nData,1)).^(1/(m-1))).*...
        %         1./sum((repmat(distances(:,c),1,nCenters)./distances).^(1/(m-1)),2);
        %                 end
        %         end

        % debug-plots for 2-d data
        %                 figure,subplot(2,2,1),plot(distances),subplot(2,2,2),plot(membership),
        %                 possibilisticClusteringPlot(subplot(2,2,3),data,membership,centers)
        %                 possibilisticClusteringPlot(subplot(2,2,4),data,distances,centers,[],2)
        %                 colormap(gray(256));

        % terminate if typicalities change by less than 0.01
        % -> check for better criterion!
        % try: typicalities that are above 1e-3 change by less than 1%
        goodTypicalityIdxL = typicality(:)>1e-3;
        %if iter>=minIter && all(abs(oldTypicality(:)- typicality(:)) < 0.01)
        if all(abs((oldTypicality(goodTypicalityIdxL)-typicality(goodTypicalityIdxL))./oldTypicality(goodTypicalityIdxL))<0.01)
            doneOptimize = true;
            %             % recalculate membership without correction
            %             membership = ...
            %             1./(1+(distances./repmat(eta,nData,1)).^(1/(m-1)));

            % return diagnostics
            diagnostics.(sprintf('iter%i',iRepeat)) = iter;

        end

    end
    switch algorithm
        case {1,2}
            if iRepeat==1 && reestimateEta

                % re-estimate eta - if change here, change below!
                eta = sum((membership.*typicality).^m.*distances,1)./sum((membership.*typicality).^m,1);

                %     % re-estimate eta (eq.10 in Krishnapuram)
                %     goodMemberIdxL = membership > 0.25;
                %     goodDistances = distances;
                %     goodDistances(~goodMemberIdxL) = 0;
                %     eta = sum(goodDistances.^2,1)./sum(goodMemberIdxL,1);
                % restart fitting
                doneOptimize = false;
                iRepeat = 2;
            else
                % exit loop
                doneRepeat = true;
            end
        case 3
            switch iRepeat
                case 1
                    % next turn: use cluster size
                    iRepeat = 2;
                    %doneOptimize = false;
                    doneRepeat = true;
                case 2
                    % check whether reestimate eta
                    if reestimateEta
                        % re-estimate eta
                        eta = sum((membership.*typicality).^m.*distances,1)./sum((membership.*typicality).^m,1);
                        iRepeat = 3;
                        doneOptimize = false;
                    else
                        doneRepeat = true;
                    end
                case 3
                    % that's it
                    doneRepeat = true;
            end
    end

end

% assign dataWeights for output
dataWeight = (membership .* typicality).^m;
dataWeight = dataWeight./max(dataWeight(:));


% reset warnings
warning(w);


%=================================

%=================================
%% INTERNAL SUBFUNCTIONS
%=================================

%% Manual initial guess
function [init, membership, centers] = manualInitialGuess(initialGuess,data,m,centerFunction,centerFunctionParameters)

% calculate distance matrix between data, show plot and let the user
% select polygons

% read distance matrix from initial guess
dm = initialGuess{1};
nData = length(dm);

% get MDS coords
try
    [scaledCoords, stress, disparities] = mdscale(dm, 2,'Criterion','stress');
    crit = 'stress';
catch
    [scaledCoords, stress, disparities] = mdscale(dm, 2,'Criterion','sstress');
    crit = 'squared stress';
end

% check for reference coordinates
if length(initialGuess) < 2
    % no reference coords
else
    % transform
    [dummy,dummy,transform] = procrustes(initialGuess{2},scaledCoords);
    scaledCoords = scaledCoords*transform.T;
end

% make scaled distance matrix for Shepard plot
scaledDm = pdist(scaledCoords);
scaledDm = squareform(scaledDm);

% Maybe we need to change that in the future, but try to plot MDS with
% Shepard plot and VAT in the same figure



fh = figure('Name','Initial guess for PCM');

% MDS plot
ah1 = subplot(2,3,[1,2,4,5]);
plot(scaledCoords(:,1),scaledCoords(:,2),'.');
text(scaledCoords(:,1),scaledCoords(:,2),num2str((1:nData)'));
set(get(ah1,'Title'),'String',sprintf('MDS plot (crit: %s)\nPlease Select ROIs here',crit));
xlabel('MDS distance')
ylabel('MDS distance')
axis equal

% Shepard plot
ah2 = subplot(2,3,3);
[dummy,sIdx] = sort(dm(:));
plot(dm(:),scaledDm(:),'.',dm(sIdx),disparities(sIdx),'-r')
set(get(ah2,'Title'),'String',sprintf('Stress: %1.3f',stress));
xlabel('dissimilarities')
ylabel('MDS distances')
axis square

% VAT
%   [sortedIdx, dmSorted] = clusterTendencyVAT(dm);
%     ah3 = subplot(2,3,6);
%     imshow(dmSorted,[0,max(dm(:))]);
%     % create shifted y-labels
%     yTickLabel = '';
%     for i=nData:-1:1
%         if isEven(i)
%             yTickLabel(i,:) = sprintf('     %3.0f',sortedIdx(i));
%         else
%             yTickLabel(i,:) = sprintf('%3.0f     ',sortedIdx(i));
%         end
%     end
%     set(ah3,'xtick',1:nData,'xtickLabel',[],...
%         'ytick',1:nData,'ytickLabel',yTickLabel,'Visible','on','TickDir','in','Position',[0.1,0.1,0.8,0.8]);
%     %set(get(ah2,'Title'),'String','sorted DM');
%     set(get(ah3,'Title'),'String',sprintf('sorted DM %2.3f - %2.3f',min(dm(:)),max(dm(:))));

% select ROIs
membership = selectRois(ah1);

if length(initialGuess) < 3 || isempty(initialGuess{3}) || ~initialGuess{3}
    % calculate centers directly
    centers = feval(centerFunction, ...
        data, membership, ones(size(membership)), m, centerFunctionParameters);
    init = 1;
else
    % do FCM for three rounds just to get a somewhat better estimate of
    % the centers (or just calculate centers once)
    init = 3;
    centers = [];
end
% close figure
close(fh)

%===================================
%% BASIC DISTANCE/CENTER FUNCTIONS
%===================================

% algoritms require data as n-by-d arrays. N data points, d dimensions

%------ FUZZY C-MEANS ---------

% distances: euclidean distance from center
% center: weighted mean of points centers.mean

%% FCM-Distance
function distances = fcmDistanceFunction(data,centers,dummy)

% distance is simple distMat
distances = distMat2(data,cat(1,centers.mean)).^2;

% for testing: p-values as distances. Note that smaller m (e.g. 1.1, but
% not 1.01) help, and that it is important to get reasonably good etas from
% the start to distinguish inliers from outliers. Eta depends on the value
% by which the distances are divided.
% distances = distMat2(data,cat(1,centers.mean));
% distances = distances./sqrt(0.01);
% distances = -log10(1-tcdf(distances,1));

%% FCM-Mean
function centers = fcmCenterFunction(data,membership,typicality,m,dummy)

% centers are the weighted means of the data points
nCenters = size(membership,2);
nDims = size(data,2);
centers(1:nCenters) = struct('mean',[]);
for c=1:nCenters
    centers(c).mean = ...
        sum(repmat(membership(:,c).*typicality(:,c),1,nDims).^m.*data,1)./...
        sum((membership(:,c).*typicality(:,c)).^m);
end


%------ GK-FCM -------

% distances: normed Mahalanobis distance
% center: struct with mean, covariance
% Krishnapuram: Eq. 11,12

%% GK-DISTANCE
function distances = gkDistanceFunction(data,centers,dummy)

% distance is normed mahalanobis distance
nCenters = length(centers);
[nData,nDims] = size(data);

distances = zeros(nData,nCenters);

for c = 1:nCenters
    % norm covariance matrix here, because there is no guarantee that the
    % input is normed
    distances(:,c) =  ...
        distMat2(data,centers(c).mean,...
        inv(centers(c).covariance)*det(centers(c).covariance)^(1/nDims)...
        ).^2;
end

%% GK-CENTER
function centers = gkCenterFunction(data,membership,typicality,m,dummy)

% find means via FCM-Mean
centers = fcmCenterFunction(data,membership,typicality,m);

% calculate covariance matrix
nCenters = length(centers);
[nData,nDims] = size(data);

for c = 1:nCenters
    centers(c).covariance = ...
        (repmat(membership(:,c).*typicality(:,c),1,nDims).^m.*...
        (data - repmat(centers(c).mean,nData,1)))' * ...
        (data - repmat(centers(c).mean,nData,1))./...
        sum((membership(:,c).*typicality(:,c)).^m);
end