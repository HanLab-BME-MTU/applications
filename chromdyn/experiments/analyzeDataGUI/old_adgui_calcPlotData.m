function [plotData,currentData,selectedTags] = adgui_calcPlotData(handles,selectedLabels,labelState)
%function to calculate the data to plot for analyzeDataGUI (plotData.xData,plotData.x_Label etc)

%get menu item number
currentData = get(handles.adgui_filename_PD,'Value')-1;

if currentData<1
    errordlg('You have to load data first!')
    return
end

anaDat = handles.data(currentData).anaDat;
dataProperties = handles.data(currentData).dataProperties;

labelCase = cat(1,handles.PD_data{labelState,4});

%init xyz-string
xyzList = 'xyz';

%init legend text
legendText = {};

%init tag list
selectedTags = {};

if size(selectedLabels,1)>1
    selectedLabels = selectedLabels';
end

for i = selectedLabels
    
    switch labelCase(i)
        case 1 %NONE
            %don't do anything
            
        case 2 %time (of frame)
            %calc data
            eval(['plotData.',xyzList(i),'Data = cat(1,anaDat.time);']);
            eval(['plotData.',xyzList(i),'_Label = ''Time [sec]'';']);
            
            stats.sigma = cat(1,anaDat.sigmaTime);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
            selectedTags{i} = [];
            
            legendText = [legendText;{'Time [sec]'}];
            
        case 3 %distance between tags
            
            [data,label,legendText,selectedTags{i},stats] = adgui_calcPlotData_distance(handles,anaDat,xyzList(i),legendText,dataProperties);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
        case 4 %rayleigh (distance)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_rayleigh(handles,anaDat,xyzList(i),legendText,dataProperties);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 5 %delta distance
            
            [data,label,legendText,selectedTags{i},stats] = adgui_calcPlotData_deltaDistance(handles,anaDat,xyzList(i),legendText);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
        case 6 %velocity delta distance
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_deltaDistanceVelocity(handles,anaDat,xyzList(i),legendText);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
           
            
        case 7 %angle between vectors
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_angleVectors(handles,anaDat,xyzList(i),legendText);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 8 %projection (parallel)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_parallelProjection(handles,anaDat,xyzList(i),legendText);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 9 %projection (perpendicular)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_perpendicularProjection(handles,anaDat,xyzList(i),legendText);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            
        case 10 %time between frames
            %calc data
            time = cat(1,anaDat.time);
            timeBF = (time(1:end-1)+time(2:end))/2;
            
            eval(['plotData.',xyzList(i),'Data = timeBF;']);
            eval(['plotData.',xyzList(i),'_Label = ''Time [sec]'';']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            legendText = [legendText;{'Time [sec]'}];
            
            selectedTags{i} = [];
            
        case 11 %displacement (uncorrected)
            
            [data,label,legendText,selectedTags{i},stats] = adgui_calcPlotData_tagDisplacementCorrected(handles,anaDat,xyzList(i),legendText,dataProperties,0);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
        case 12 %velocity (uncorrected)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagVelocityUncorrected(handles,anaDat,xyzList(i),legendText);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            
        case 13 %diffusion (uncorrected)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyzList(i),legendText,0);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 14 %displacement (centroid)
            
            [data,label,legendText,selectedTags{i},stats] = adgui_calcPlotData_tagDisplacementCorrected(handles,anaDat,xyzList(i),legendText,dataProperties,1);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
        case 15 %velocity (centroid corrected)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagVelocityCorrected(handles,anaDat,xyzList(i),legendText,1);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 16 %diffusion (centroid corrected)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyzList(i),legendText,1);
            
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            
        case 17 %timepoint (per frame)
            %get data
            timePoints = cat(1,anaDat.timePoint);
            
            eval(['plotData.',xyzList(i),'Data = timePoints;']);
            eval(['plotData.',xyzList(i),'_Label = ''Timepoints'';']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            legendText = [legendText;{'Timepoints'}];
           
            selectedTags{i} = [];
            
        case 18 %timepoint (between frames)
            %calc data
            timePoints = cat(1,anaDat.timePoint);
            timePointsBF = (timePoints(1:end-1)+timePoints(2:end))/2;
            
            eval(['plotData.',xyzList(i),'Data = timePointsBF;']);
            eval(['plotData.',xyzList(i),'_Label = ''Timepoints'';']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            legendText = [legendText;{'Timepoints'}];
            
            selectedTags{i} = [];
            
        case 19 %timesteps for diffusion analysis
            %get data
            timePoints = cat(1,anaDat.timePoint);
            
            %change into seconds
            time = timePoints*dataProperties.timeLapse;
            
            eval(['plotData.',xyzList(i),'Data = time(1:floor(length(time)/2));']);
            eval(['plotData.',xyzList(i),'_Label = ''Time [sec]'';']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
            selectedTags{i} = [];
            
            legendText = [legendText;{'delta Time [sec]'}];
            
        case 20 %displacement (tag centered)
            
            [data,label,legendText,selectedTags{i},stats] = adgui_calcPlotData_tagDisplacementCorrected(handles,anaDat,xyzList(i),legendText,dataProperties,2);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = stats;']);
            
       case 21 %velocity (tag centered)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagVelocityCorrected(handles,anaDat,xyzList(i),legendText,2);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
        case 22 %diffusion (tag centered)
            
            [data,label,legendText,selectedTags{i}] = adgui_calcPlotData_tagDiffusion(handles,anaDat,xyzList(i),legendText,2);
            
            eval(['plotData.',xyzList(i),'Data = data;']);
            eval(['plotData.',xyzList(i),'_Label = label;']);
            eval(['plotData.',xyzList(i),'Stats = [];']);
            
    end %switch
    
end %for xyz

plotData.legendText = char(legendText);