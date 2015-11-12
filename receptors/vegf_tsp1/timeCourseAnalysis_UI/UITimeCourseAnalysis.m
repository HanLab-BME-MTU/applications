function [] = UITimeCourseAnalysis()
%TimeCourseAnalysis of user selected CombinedMovieList objects
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%
%SYNOPSIS [] = UITimeCourseAnalysis()
%
%Tae H Kim, July 2015

%% Initialize
%Progresstext
clear progressTextMultiple;

p = TimeCourseAnalysisConfig;

if(~isempty(p))
    disp(p);
    timeCourseAnalysis(p.CML_FullPath, p.outputDir, ...
        'doNewAnalysis', p.doNewAnalysis, ...
        'doPartitionAnalysis', p.partitioningAnalysis, ...
        'start2zero', p.start2zero, ...
        'shiftPlotPositive', p.shiftPlotPositive, ...
        'channelNames', p.channelTable(:,2), ...
        'channels',find([p.channelTable{:,1}]));
end


% %% Prompt user
% %prompt user to select a folder where all figures and data will be stored
% outputDir = uigetdir('', 'Select output folder');
% %prompt user to select Combined Movie Data objects
% %until they press cancel.
% CML_FullPath = {};
% [fileName, filePath] = uigetfile('*.mat', 'Select CombinedMovieLists', 'MultiSelect', 'on');
% while ~isnumeric(fileName)
%     if iscell(fileName)
%         CML_FullPath = [CML_FullPath cellfun(@(x) [filePath x], fileName, 'UniformOutput', false)]; %#ok<*AGROW>
%     else
%         CML_FullPath{end+1} = [filePath fileName];
%     end
%     [fileName, filePath] = uigetfile('*.mat', 'Select CombinedMovieLists', 'MultiSelect', 'on');
% end
%prompts user if new timeCourse analysis should be done
% userChoiceNTCA = listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', {'Use old timeCourseAnalysis if possible', 'Do new timeCourseAnalysis'}, 'ListSize', [300, 300]);
% if userChoiceNTCA == 1
%     doNewAnalysis = false;
% else
%     doNewAnalysis = true;
% end
%prompts user for analysis parameters
%the function will call actual timeCourseAnalysis
% parameterCheckGUI(outputDir, CML_FullPath);
end
% %% Parameter Checkbox GUI
% function parameterCheckGUI(outputDir, CML_FullPath)
% % % Create figure
% % h.f = figure('units','pixels','position',[400,400,500,100],...
% %              'toolbar','none','menu','none','name','Choose TimeCourseAnalysis to be done');
% % % Create parameter checkboxes
% % h.c(1) = uicontrol('style','checkbox','units','pixels',...
% %                 'position',[10,40,450,20],'string','Partitioning Analysis');
% % h.c(2) = uicontrol('style','checkbox','units','pixels',...
% %                 'position',[10,70,450,20],'string','Set Average Start Time to Zero');
% % % Create OK pushbutton   
% % h.p = uicontrol('style','pushbutton','units','pixels',...
% %                 'position',[40,5,70,20],'string','OK',...
% %                 'callback',@p_call);
%     h = openfig('TimeCourseAnalysisConfig.fig');
%     % Pushbutton callback
%     function p_call(varargin)
%         parameter = get(h.c,'Value');
%         channelNames = get(h.channelNames,'Value');
%         %closes the dialogue box
%         close(h.f);
%         clear progressTextMultiple;
%         pause(1);
%         %calls the function that does the timeCourseAnalysis
%         timeCourseAnalysis(CML_FullPath, outputDir, 'doNewAnalysis', doNewAnalysis, 'doPartitionAnalysis', parameter{1}, 'start2zero', parameter{2}, 'channelNames', channelNames);
%     end
% end

