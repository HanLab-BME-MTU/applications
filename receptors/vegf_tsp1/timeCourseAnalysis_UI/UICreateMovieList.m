function [ML,param] = UICreateMovieList(varargin)
%Creates ML and MD if it doesn't already exist of movies selected by user
%Prompts the user to select movie (.tif files) that correspond to MD.
%This takes advantage of the fact that uTrackPackageGUI will always create
%the track analysis file in same folder as the movie itself.
%Tae H Kim, July 2015

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.StructExpand = true;
ip.addParameter('pixelSize_', 90, @isnumeric);
ip.addParameter('timeInterval_', 0.1, @isnumeric);
ip.addParameter('numAperature_', 1.49, @isnumeric);
ip.addParameter('emissionWavelength_', [], @isnumeric);%590nm for Rhod Red X, 525nm for GFP
ip.addParameter('exposureTime_', 20, @isnumeric);
ip.addParameter('imageType_', {'TIRF'}, @(x) ischar(x) || iscellstr(x));
% Parameters obtained by user input
ip.addParameter('fileNameML',[],@ischar);
ip.addParameter('filePathML',[],@ischar);
ip.addParameter('fileName',[],@iscellstr);
ip.addParameter('filePath',[],@ischar);
ip.addParameter('relTimeZero',[]);
ip.addParameter('relTimeZero2',[]);
ip.addParameter('zeroSelect',[]);
ip.addParameter('zeroSelect2',[]);
ip.addParameter('addTZ1',[]);
ip.addParameter('addTZ2',[]);
ip.parse(varargin{:});
param = ip.Results;
% emissionWL = ip.Results.emissionWavelength_;
ML = [];

if(~iscellstr(param.imageType_))
    param.imageType_ = {param.imageType_};
end

%% Ask for emission wavelength
if isempty(param.emissionWavelength_)
    canceled = false;
    while(~canceled)
        [param.emissionWavelength_(end+1),emissionStr,canceled] = timeCourseAnalysis.UI.waveLengthPrompt;
        if(~canceled && isempty(param.emissionWavelength_(end)))
            param.imageType_{length(param.emissionWavelength_)} = '';
        end
        if(canceled)
            param.emissionWavelength_ = param.emissionWavelength_(1:end-1);
        else
            disp([ emissionStr ' selected. Select another wavelength for multichannel movies or cancel to continue.']);
        end
    end
end

%% MovieList creation part 1
if(isempty(param.fileNameML) || isempty(param.filePathML))
    [param.fileNameML, param.filePathML] = uiputfile('*.mat', ... 
        'Find a place to save your movie list');
end
% remove filesep (robust enough?)
outputDir = param.filePathML(1:end-1);
%outputDir = uigetdir(filePathML, 'Select the directory to store the list analysis output');

%% MovieData creation
%User file selection for MD creation
if(isempty(param.fileName) || isempty(param.filePath))
    [param.fileName, param.filePath] = uigetfile( ...
        {'*.tif','*.ome.tiff','*488.ome.tiff','*561.ome.tiff','*640.ome.tiff'}', ...
        'Select Movies', 'MultiSelect', 'on');
end
if(~iscellstr(param.fileName))
    param.fileName = {param.fileName};
end

%Moved this prompts to beginning to keep all user input together
%% Relative time zero selection
if(isempty(param.relTimeZero))
    [param.relTimeZero, param.addTZ1, param.zeroSelect] = timeCourseAnalysis.relativeTimeZeroSelection('Enter start time',[],param);
end

%% Relative time zero selection2
if(isempty(param.relTimeZero2))
    [param.relTimeZero2, param.addTZ2, param.zeroSelect2] = timeCourseAnalysis.relativeTimeZeroSelection('Enter VEGF addition time',{'none'},param);
end

try
    %User file selection
%     MD(length(param.fileName)) = MovieData;
    MD = cell(1,length(param.fileName));
    movies = cell(1,length(MD));
    for iMD = 1:length(MD)
        printLength = fprintf('MD %g/%g\n', iMD, length(MD));

        evalc('[MD{iMD}, movies{iMD}] = timeCourseAnalysis.configureMovie(param.fileName{iMD},param.filePath,param)');       

        fprintf(repmat('\b',1,printLength));
    end

    %% MovieList creation part 2
    ML = MovieList(MD, outputDir, 'movieListFileName_', param.fileNameML, 'movieListPath_', param.filePathML(1:end-1), 'createTime_', clock());
    %%% evalc('ML.sanityCheck();');
    ML.sanityCheck();

    %% Obtain relTimeZero{2} if needed
    if( param.zeroSelect  < 0)
        param.relTimeZero  = MD{-param.zeroSelect }.acquisitionDate_;
    end
    if( param.zeroSelect2 < 0)
        param.relTimeZero2 = MD{-param.zeroSelect2}.acquisitionDate_;
    end

    %% Save ML
    ML.addProcess(TimePoints(ML));
    if param.addTZ1
        ML.processes_{1}.addTimePoint(param.relTimeZero, 'start');
    end
    if param.addTZ2
        ML.processes_{1}.addTimePoint(param.relTimeZero2, 'VEGF_added');
    end
    ML.save;
catch err
    disp(getReport(err));
    warning('UICreateMovieList failed.');
end

end
