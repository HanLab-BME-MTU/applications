function [] = UITimeCourseAnalysis()
%TimeCourseAnalysis of user selected CombinedMovieList objects
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%
%SYNOPSIS [] = UITimeCourseAnalysis()
%
%Tae H Kim, July 2015

%% Prompt user
%prompt user to select a folder where all figures and data will be stored
outputDir = uigetdir('', 'Select output folder');
%prompt user to select Combined Movie Data objects
%until they press cancel.
CML_FullPath = {};
[fileName, filePath] = uigetfile('*.mat', 'Select CombinedMovieLists', 'MultiSelect', 'on');
while ~isnumeric(fileName)
    if iscell(fileName)
        CML_FullPath = [CML_FullPath cellfun(@(x) [filePath x], fileName, 'UniformOutput', false)]; %#ok<*AGROW>
    else
        CML_FullPath{end+1} = [filePath fileName];
    end
    [fileName, filePath] = uigetfile('*.mat', 'Select CombinedMovieLists', 'MultiSelect', 'on');
end
%prompts user for analysis parameters
%the function will call actual timeCourseAnalysis
parameterCheckGUI(outputDir, CML_FullPath);
end
%% Parameter Checkbox GUI
function parameterCheckGUI(outputDir, CML_FullPath)
% Create figure
h.f = figure('units','pixels','position',[400,400,500,100],...
             'toolbar','none','menu','none','name','Choose TimeCourseAnalysis to be done');
% Create parameter checkboxes
h.c(1) = uicontrol('style','checkbox','units','pixels',...
                'position',[10,40,450,20],'string','Partitioning Analysis');
h.c(2) = uicontrol('style','checkbox','units','pixels',...
                'position',[10,70,450,20],'string','Testing');    
% Create OK pushbutton   
h.p = uicontrol('style','pushbutton','units','pixels',...
                'position',[40,5,70,20],'string','OK',...
                'callback',@p_call);
    % Pushbutton callback
    function p_call(varargin)
        parameter = get(h.c,'Value');
        %closes the dialogue box
        close(h.f);
        %calls the function that does the timeCourseAnalysis
        timeCourseAnalysis(CML_FullPath, outputDir, 'doPartitionAnalysis', parameter{1});
    end
end

