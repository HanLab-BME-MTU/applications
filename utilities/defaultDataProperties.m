function dataProperties = defaultDataProperties(varargin);
%DEFAULTDATAPROPERTIES loads default data properties and updates them
%
% SYNOPSIS dataProperties = defaultDataProperties(structure)
%          dataProperties = defaultDataProperties('propertyName',propertyValue)
%
% INPUT    structure : dataProperties-substructure or movieHeader structure
%          'propertyName'/propertyValue : name/value pairs for structure
%
% For defaults, the program will first look for a file called
% 'default_dataProperties' in the current directory
%
% OUTPUT   dataProperties-structure
%
% c: 8/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================
% TEST INPUT
%========================

% inputType: 0 - none, 1 - dataProperties, 2 - movieHeader, 3 - pairs
inputType = 0;
% structure or pairs?
if nargin == 1 && isstruct(varargin{1})
    % structure
    inputStruct = varargin{1};

    % dataProperties or movieHeader?
    if isfield(inputStruct,'pixelX') && isfield(inputStruct,'numCols')
        inputType = 2;
    else
        inputType = 1;
    end

elseif isEven(nargin) && nargin > 0

    % we assume that there are name/value pairs. They will be sorted out
    % later
    inputType = 3;

else

    % inputType stays 0;
end

%=========================


%=========================
% READ DEFAULTS
%=========================

% 1. try to read from file in current path
file = dir('default_dataProperties');
if ~isempty(file)
    % load default_DataProperties
    load(file.name)
else
    dataProperties = loadDefault;
end

%==========================


%==========================
% ADJUST PARAMETERS
%==========================

switch inputType
    case 0
        % do nothing
    case 1
        % dataProperties substructure
        fn = fieldnames(inputStruct);
        for i = 1:length(fn)
            dataProperties.(fn{i}) = inputStruct.(fn{i});
        end

    case 2
        % movieHeader
        dataProperties.PIXELSIZE_XY = inputStruct.pixelX;
        if inputStruct.pixelX ~= inputStruct.pixelY
            error('unequal pixelsize x/y!')
        end
        dataProperties.PIXELSIZE_Z = inputStruct.pixelZ;
        dataProperties.LENSID = inputStruct.lensID;
        dataProperties.NA = naFromLensID(dataProperties.LENSID,1);
        dataProperties.movieSize(2) = inputStruct.numCols;
        dataProperties.movieSize(1) = inputStruct.numRows;
        dataProperties.movieSize(3) = inputStruct.numZSlices;
        dataProperties.movieSize(4) = inputStruct.numTimepoints;
        dataProperties.WVL = inputStruct.wvl;
        if isfield(inputStruct, 'Time')
            dataProperties.frameTime = reshape(inputStruct.Time,...
                dataProperties.movieSize(3:4));
        end
        % time could also be per frame only (metamorph)
        if isfield(inputStruct,'frameTime')
            dataProperties.frameTime = inputStruct.frameTime;
        end
        if isfield(inputStruct,'expTime')
            dataProperties.expTime = inputStruct.expTime;
        end
        if isfield(inputStruct,'ndFilter')
            dataProperties.NDfilter = inputStruct.ndFilter;
        end
        if isfield(inputStruct,'cropInfo')
            dataProperties.crop = inputStruct.cropInfo
        end

    case 3
        % name/value pairs
        for i = 1:2:nargin
            propertyName = varargin{i};
            if ~isempty(propertyName) && ischar(propertyName)
                dataProperties.(propertyName) = varargin{i+1};
            end
        end

end

% calculate filterparms
if isempty(dataProperties.NA)
    warning('DefaultDataProperties:UsingDefaultNA',...
        'Using default NA for calculating filter parameters!')
end
[FT_XY, FT_Z] = calcFilterParms(...
    dataProperties.WVL,dataProperties.NA,1.51,'gauss',[], ...
    [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
dataProperties.FT_SIGMA = [FT_XY,FT_XY,FT_Z];

% calculate frameTime
if ~isfield(dataProperties,'frameTime')
    movieLength = dataProperties.movieSize(4);
    stackSize = dataProperties.movieSize(1:3);
    meanTime = [1:movieLength]';
    %center acquisition timepoints around the mean, allow 0.03 sec per slice
    stackTime = ((stackSize(3)-1)/2-[stackSize(3)-1:-1:0])*0.03;
    frameTime = repmat(meanTime,1,stackSize(3))+repmat(stackTime,movieLength,1);
    %frameTime starts with 0
    dataProperties.frameTime = frameTime-frameTime(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataProperties = loadDefault
dataProperties.PIXELSIZE_XY = 0.05;
dataProperties.PIXELSIZE_Z = 0.2;
dataProperties.movieSize = [64 64 16 10];


dataProperties.LENSID=12003;
dataProperties.NA=1.4000;
dataProperties.WVL=0.5250;
dataProperties.timeLapse=1;
dataProperties.expTime=NaN;
dataProperties.NDfilter=NaN;
dataProperties.cellCycle=NaN;
dataProperties.strains=NaN;
dataProperties.drugs=NaN;
dataProperties.temperature={'Nan'};
dataProperties.crop=[];
dataProperties.maxSize=200*2e6;
dataProperties.F_TEST_PROB=0.9000;
dataProperties.IDopt= [];
dataProperties.PATCHSIZE=7;
dataProperties.CH_MAXNUMINTERV=1000;
dataProperties.OVERLPSIZE=[15 15 15];
dataProperties.sigmaCorrection=[1 1];
dataProperties.split=[];
dataProperties.MAXSPOTS=5;
dataProperties.T_TEST_PROB=0.0500;
dataProperties.help=['synthMovie'];
dataProperties.maxSize = 200 * 2^20; % 200 Mb
dataProperties.amplitudeCutoff = 0; % undefined
dataProperties.fitNPlusOne = 1; % super-resolution fitting

% linker properties
dataProperties.linker_relativeMaxDistance = -1;
dataProperties.linker_absoluteMaxDistance=-1;
dataProperties.linker_relAmpWeight=1/1.5;
dataProperties.linker_useCOM = 1;







