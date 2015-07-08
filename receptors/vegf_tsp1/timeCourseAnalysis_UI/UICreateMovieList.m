function [ML] = UICreateMovieList(varargin)
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
ip.addOptional('pixelSize_', 90, @isnumeric);
ip.addOptional('timeInterval_', 0.1, @isnumeric);
ip.addOptional('numAperature_', 1.49, @isnumeric);
ip.addOptional('emissionWavelength_', [], @isnumeric);%590nm for Rhod Red X, 525nm for GFP
ip.addOptional('exposureTime_', 20, @isnumeric);
ip.addOptional('imageType_', 'TIRF', @isstr);
ip.parse(varargin{:});
emissionWL = ip.Results.emissionWavelength_;
%Ask for emission wavelength
if isempty(emissionWL)
    stringListWL = {'590nm : Rhod Red X', '525nm : GFP'};
    userChoiceWL = listdlg('PromptString','Select wavelength:', 'SelectionMode','single', 'ListString', stringListWL);
    if userChoiceWL == 1
        emissionWL = 590;
    end
    if userChoiceWL == 2
        emissionWL = 525;
    end
end
%% MovieList creation part 1
[fileNameML, filePathML] = uiputfile('*.mat', 'Find a place to save your movie list');
outputDir = uigetdir(filePathML, 'Select the directory to store the list analysis output');
%% MovieData creation part 1
%User file selection for MD creation
[fileName, filePath] = uigetfile('*.tif', 'Select Movies', 'MultiSelect', 'on');
%% Relative time zero selection
%This is for timecourse analysis
%Obtain relative time zero from user
%User can enter 6 element array [yr month day hr min sec]
% or scalar of the index number of the MD to be used as the relative time zero
% or 'select' to bring up another dialogue box with list dialogue box later
% or 'min' to use MD with earliest acquisition / observation time as time zero
userInputStr = inputdlg('6 element time array -or- index number -or- ''select''', 'Enter start time', 1, {'2015 7 2 14 44 18.9'});
userInputNum = str2num(userInputStr{1});
%if time array
if numel(userInputNum) == 6
    relTimeZero = userInputNum;
    zeroSelect = 0; %done no need to do again
end
%if index number
if isscalar(userInputNum)
    if userInputNum > numel(fileName)
        error('Input index for relative time zero out of bounds');
    end
    zeroSelect = 1; %need to get relTime0
end
%if not a number
if isempty(userInputNum)
    %if select
    if strcmpi(userInputStr, 'select')
        zeroSelect = 2; %need to get relTime0
    end
    %if min
    %if strcmpi(userInputStr, 'min')
    %    zeroSelect = 3; %need to get relTime0
    %end
end
%% Relative time zero selection2
%'same' means same as time above
userInputStr2 = inputdlg('6 element time array -or- index number -or- ''select'' -or- ''no VEGF''', 'Enter VEGF addition time', 1, {'no VEGF'});
userInputNum2 = str2num(userInputStr2{1});
addTZ2 = true;
%if time array
if numel(userInputNum2) == 6
    relTimeZero2 = userInputNum2;
    zeroSelect2 = 0; %done no need to do again
end
%if index number
if isscalar(userInputNum2)
    if userInputNum2 > numel(fileName)
        error('Input index for relative time zero out of bounds');
    end
    zeroSelect2 = 1; %need to get relTime0
end
%if not a number
if isempty(userInputNum2)
    %if select
    if strcmpi(userInputStr2, 'select')
        zeroSelect2 = 2; %need to get relTime0
    end
    if strcmpi(userInputStr2, 'no VEGF')
        addTZ2 = false;
        zeroSelect2 = 0; %done no need to do again
    end
    %if min
    %if strcmpi(userInputStr, 'min')
    %    zeroSelect = 3; %need to get relTime0
    %end
end
%% MovieData creation part2
%User file selection
if iscellstr(fileName)
    nMD = length(fileName);
    iMD = 1;
    for MCounter = nMD:-1:1
        printLength = fprintf('MD %g/%g\n', iMD, nMD);
        name = fileName{MCounter}(1:end-4);
        movies{MCounter} = [filePath name filesep name '.mat'];
        %create new MD if one doesn't exist already
        if exist(movies{MCounter}, 'file') == 0
            evalc('MD = MovieData([filePath fileName{MCounter}])');
            MD.pixelSize_ = ip.Results.pixelSize_;
            MD.timeInterval_ = ip.Results.timeInterval_;
            MD.numAperture_ = ip.Results.numAperature_;
            evalc('MD.sanityCheck;');
            MD.channels_.emissionWavelength_ = emissionWL;
            MD.channels_.exposureTime_ = ip.Results.exposureTime_;
            MD.channels_.imageType_ = ip.Results.imageType_;
            evalc('MD.channels_.sanityCheck;');
            MD.save;
        end
        fprintf(repmat('\b',1,printLength));
        iMD = iMD + 1;
    end
else
    name = fileName(1:end-4);
    movies = [filePath name filesep name '.mat'];
end
%% MovieList creation part 2
ML = MovieList(movies, outputDir, 'movieListFileName_', fileNameML, 'movieListPath_', filePathML(1:end-1), 'createTime_', clock());
evalc('ML.sanityCheck();');
%% Relative time zero set
%if index number
if zeroSelect == 1
    evalc('MD = MovieData.load([filePath fileName{userInputNum}])');
    relTimeZero = MD.acquisitionDate_;
end
%if select
if zeroSelect == 2
    userChoiceMD = listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', fileName);
    evalc('MD = MovieData.load([filePath fileName{userChoiceMD}])');
    relTimeZero = MD.acquisitionDate_;
end
%% Relative time zero 2 set
%If one is adding clathrin 
if zeroSelect2 == 1
    evalc('MD = MovieData.load([filePath fileName{userInputNum2}])');
    relTimeZero2 = MD.acquisitionDate_;
end
%if select
if zeroSelect2 == 2
    userChoiceMD2 = listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', fileName);
    evalc('MD = MovieData.load([filePath fileName{userChoiceMD2}])');
    relTimeZero2 = MD.acquisitionDate_;
end
%% Save relative time
ML.addProcess(TimePoints(ML));
ML.processes_{1}.addTimePoint(relTimeZero, 'start');
if addTZ2
    ML.processes_{1}.addTimePoint(relTimeZero2, 'VEGF_added');
end
ML.save;
end