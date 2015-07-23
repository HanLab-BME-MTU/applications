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
ip.addParameter('pixelSize_', 90, @isnumeric);
ip.addParameter('timeInterval_', 0.1, @isnumeric);
ip.addParameter('numAperature_', 1.49, @isnumeric);
ip.addParameter('emissionWavelength_', [], @isnumeric);%590nm for Rhod Red X, 525nm for GFP
ip.addParameter('exposureTime_', 20, @isnumeric);
ip.addParameter('imageType_', 'TIRF', @isstr);
ip.parse(varargin{:});
param = ip.Results;
emissionWL = ip.Results.emissionWavelength_;
%Ask for emission wavelength
if isempty(emissionWL)
    stringListWL = {'525nm : Alexa 488', '530nm : GFP', '590nm : Rhod Red X', '668nm : Alexa 640', '669nm : Atto 547N', 'Brightfield'};
    userChoiceWL = listdlg('PromptString','Select wavelength:', 'SelectionMode','single', 'ListString', stringListWL);
    if userChoiceWL == 1
        emissionWL = 525;
    end
    if userChoiceWL == 2
        emissionWL = 530;
    end
    if userChoiceWL == 3
        emissionWL = 590;
    end
    if userChoiceWL == 4
        emissionWL = 668;
    end
    if userChoiceWL == 5
        emissionWL = 669;
    end
    if userChoiceWL == 6
        emissionWL = [];
        param.imageType_ = 'Brightfield';
    end
end
%% MovieList creation part 1
[fileNameML, filePathML] = uiputfile('*.mat', 'Find a place to save your movie list');
outputDir = filePathML(1:end-1);
%outputDir = uigetdir(filePathML, 'Select the directory to store the list analysis output');
%% MovieData creation
%User file selection for MD creation
[fileName, filePath] = uigetfile('*.tif', 'Select Movies', 'MultiSelect', 'on');
%User file selection
if iscellstr(fileName)
    nMD = length(fileName);
    iMD = 1;
    for MCounter = nMD:-1:1
        printLength = fprintf('MD %g/%g\n', iMD, nMD);
        name = fileName{MCounter}(1:end-4);
        movies{MCounter} = [filePath name filesep name '.mat'];
        %create new MD if one doesn't exist already
        %use evalc to silence the output
        if exist(movies{MCounter}, 'file') == 0
            evalc('MD(MCounter) = MovieData([filePath fileName{MCounter}])');
            MD(MCounter).pixelSize_ = param.pixelSize_;
            MD(MCounter).timeInterval_ = param.timeInterval_;
            MD(MCounter).numAperture_ = param.numAperature_;
            evalc('MD(MCounter).sanityCheck;');
            MD(MCounter).channels_.emissionWavelength_ = emissionWL;
            MD(MCounter).channels_.exposureTime_ = param.exposureTime_;
            MD(MCounter).channels_.imageType_ = param.imageType_;
            evalc('MD(MCounter).channels_.sanityCheck;');
            MD(MCounter).save;
        else
            evalc('MD(MCounter) = MovieData.load(movies{MCounter})');
            evalc('MD(MCounter).channels_.sanityCheck;');
        end
        fprintf(repmat('\b',1,printLength));
        iMD = iMD + 1;
    end
else
    name = fileName(1:end-4);
    movies = {[filePath name filesep name '.mat']};
end
%% MovieList creation part 2
ML = MovieList(movies, outputDir, 'movieListFileName_', fileNameML, 'movieListPath_', filePathML(1:end-1), 'createTime_', clock());
evalc('ML.sanityCheck();');
%% Relative time zero selection
%This is for timecourse analysis
%Obtain relative time zero from user
%User can enter 6 element array [yr month day hr min sec]
% or scalar of the index number of the MD to be used as the relative time zero
% or 'select' to bring up another dialogue box with list dialogue box later
% or 'min' to use MD with earliest acquisition / observation time as time zero
userInputStr = inputdlg('6 element time array -or- index number -or- ''select'' -or- ''no start''', 'Enter start time', 1, {'2015 7 2 14 44 18.9'});
userInputNum = str2num(userInputStr{1}); %#ok<ST2NM>
addTZ1 = true;
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
    if strcmpi(userInputStr, 'no start')
        addTZ1 = false;
        zeroSelect = 0; %done no need to do again
    end
    %if min
    %if strcmpi(userInputStr, 'min')
    %    zeroSelect = 3; %need to get relTime0
    %end
end
%% Relative time zero set
%if index number
if zeroSelect == 1
    %evalc('MD = MovieData.load(movies{userInputNum})');
    relTimeZero = MD(userInputNum).acquisitionDate_;
end
%if select
if zeroSelect == 2
    userChoiceMD = listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', fileName);
    %evalc('MD = MovieData.load(movies{userChoiceMD})');
    relTimeZero = MD(userChoiceMD).acquisitionDate_;
end
%% Relative time zero selection2
%'same' means same as time above
userInputStr2 = inputdlg('6 element time array -or- index number -or- ''select'' -or- ''no VEGF''', 'Enter VEGF addition time', 1, {'no VEGF'});
userInputNum2 = str2num(userInputStr2{1}); %#ok<ST2NM>
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
%% Relative time zero 2 set
%If one is adding clathrin 
if zeroSelect2 == 1
    %evalc('MD = MovieData.load(movies{userInputNum2})');
    relTimeZero2 = MD(userInputNum2).acquisitionDate_;
end
%if select
if zeroSelect2 == 2
    userChoiceMD2 = listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', fileName);
    %evalc('MD = MovieData.load(movies{userChoiceMD2})');
    relTimeZero2 = MD(userChoiceMD2).acquisitionDate_;
end
%% Save ML
ML.addProcess(TimePoints(ML));
if addTZ1
    ML.processes_{1}.addTimePoint(relTimeZero, 'start');
end
if addTZ2
    ML.processes_{1}.addTimePoint(relTimeZero2, 'VEGF_added');
end
ML.save;
end
