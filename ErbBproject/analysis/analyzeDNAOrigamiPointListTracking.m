function success=analyzeDNAOrigamiPointListTracking(list,name,varargin)
%ANALYZEDNAOrigamiPointListTracking analyzes points list from Wyss DNA
%   paint experiments and applys Khulouds tracking software
%
%   Input:
%         required arguments
%           dataDirectory  -> directory of experimental data
%              wavelength  -> wavelength(s) [nm], must be in proper ordering
%
%         optional arguments
%            NA  ->  numerical aperture, default: 1.49
%           MAG  ->  magnification, default: 100x
%           rep  ->  number of frames in activation cycle, default: 10
%    createMask  ->  should a mask be created? default: false
%         doMMF  ->  use mixture model fitting default: false
%        GapLen  ->  maximum gap size to be closed, default: 4
%        Radius  ->  maxium search radius for linking defalt: 2
%  
%   Output:
%         success  ->  1/0: analysis was successful/unsuccessful
%
% Ulrich Schmidt, March 13, 2012
% Jeffrey Werbin, March 2013

success=1;

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('list',@isnumeric);
ip.addRequired('name',@ischar);

ip.addOptional('NA',1.49,@isscalar);
ip.addOptional('MAG',100,@isscalar);
ip.addOptional('rep',10,@isnumeric);
ip.addOptional('createMask',false,@islogical);
ip.addOptional('doMMF',false,@islogical);
ip.addOptional('GapLen',4,@isnumeric);
ip.addOptional('Radius',2,@isnumeric);

ip.parse(list,name,varargin{:});

NA=ip.Results.NA;
MAG=ip.Results.MAG;
rep=ip.Results.rep;
createMask=ip.Results.createMask;
doMMF=ip.Results.doMMF;
GapLen=ip.Results.GapLen;
Radius=ip.Results.Radius;

%% Creates parameter structures for tracking
%% general gap closing parameters
gapCloseParam.timeWindow = GapLen; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatStationaryLink';

%parameters
parameters.searchRadius = Radius;

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatStationaryCloseGaps';

%parameters
parameters.searchRadius = Radius; 
parameters.gapPenalty = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions = [];

%% additional input
saveResults = 0;

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% Reformat data and run tracking

    sigma = 0.3; %approximately 12 nm assuming
    amp = 10;
    amp_sig = 1;

    
    ind = hist(list(:,1),1:max(list(:,1)));
    ind = cumsum(ind);
    ind = [0,ind];
    

    
    %shift=10;
    numframe = numel(ind)-1;
    
    fprintf(1,'list being analyzed: %s\n',name);
    for i=1:numframe
        try
        temp = list(ind(i)+1:ind(i+1),2:3);
        catch
            i
        end
       if numel(temp) > 0
        movieInfo(i).xCoord = [temp(:,1),sigma*ones(size(temp(:,1)))];
        movieInfo(i).yCoord = [temp(:,2),sigma*ones(size(temp(:,2)))];
        movieInfo(i).amp = [amp*ones(size(temp(:,2))),amp_sig*ones(size(temp(:,2)))];
       else
        movieInfo(i).xCoord = [[],[]];
        movieInfo(i).yCoord = [[],[]];
        movieInfo(i).amp = [[],[]];
       end
           
        
        
        progressText(i/numframe,'All work and no play makes Jack a dull boy');
    end
    
    fprintf(1,'\n');
    
    %Apply tracking and Gap closing to localized data
    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
    
    
    % save results
    d = cd()
    save([d filesep name],'features','tracksFinal','GapLen','Radius');
        
end