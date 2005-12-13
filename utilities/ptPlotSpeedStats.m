function ptPlotSpeedStats (radioButtons, imageName, SaveDir, xAxis, avgVelocityStats, windowSize, drugTimepoint, handles, finalPersistence)
% ptPlotSpeedStats generates the plots for the velocity statistics
%
% SYNOPSIS       ptPlotSpeedStats (radioButtons, imageName, SaveDir, xAxis, velocityStats)
%
% INPUT          radioButtons : struct containing the value of all gui radiobuttons
%                imageName : name of the first image in the movie (used as title)
%                SaveDir : name of the directory where the plots are saved in files
%                xAxis : matrix with frame numbers (this should have the same length as
%                        (all the other matrices that follow)
%                velocityStats : struct with the following fields:
%
%                   avgVelocity : vector with avg velocity all cells
%                   avgVelocitySquared : vector with squared avg velocity all cells
%                   avgSingleVelocity : vector with avg single cell vel.
%                   avgClusteredVelocity : vector with avg clustered cell vel.
%                windowSize : size of the averaging window 
%                drugTimepoint : frame number where the drug/EGF is applied
%
% OUTPUT         None (plots are directly shown on the screen and written to disk) 
%
% DEPENDENCIES   ptPlotsSpeedStats  uses { nothing }
%                                  
%                ptPlotSpeedStats is used by { Polytrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Rewrite of all plotting functions

%get the global data set in ptCalculateSpeedValues
global finalPersistence
global jrSaveDir
global jrimageName
jrSaveDir=SaveDir;
jrimageName=imageName;

% Get data from struct
avgVelocity = avgVelocityStats.avgVelocity;
avgVelocitySquared = avgVelocityStats.avgVelocitySquared;
avgSingleVelocity = avgVelocityStats.avgSingleVelocity;
avgClusteredVelocity = avgVelocityStats.avgClusteredVelocity;
avgSinglePersistence=avgVelocityStats.avgSinglePersistence;
avgPersistence=avgVelocityStats.avgPersistence;

% Calculate average values if we have to
if radioButtons.runningaverage    
    raAvgVelocity = movingAverage (avgVelocity, windowSize, 'median');
    raAvgVelocitySquared = movingAverage (avgVelocitySquared, windowSize, 'median');
    raAvgSingleVelocity = movingAverage (avgSingleVelocity, windowSize, 'median');
    raAvgClusteredVelocity = movingAverage (avgClusteredVelocity, windowSize, 'median');
    raAvgPersistence=movingAverage (avgPersistence, windowSize, 'median');
    raAvgSinglePersistence=movingAverage (avgSinglePersistence, windowSize, 'median');
end

% Here is where all the plotting starts
if radioButtons.speedplot_2

    if ~radioButtons.donotshowplots
        
        % Generate the avg velocity plot (all cells)
        h_fig = figure('Name', imageName);

        % Draw a plot showing average velocity of all cells
        ymax = max (avgVelocity) + 1;
        subplot (2,1,1); plot (xAxis, avgVelocity); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgVelocity, 1, [2 1 1], drugTimepoint);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveAllVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesAllVelocity.csv'], est);
           end
        end
        
        title ('Avg Velocity All Cells');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        ymax = max (avgVelocitySquared) + 1;
        subplot (2,1,2); plot (xAxis, avgVelocitySquared); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgVelocitySquared, 'r'); hold off;
        end
                
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgVelocitySquared, 1, [2 1 2], drugTimepoint);
           hold off;
        end        
        
        title ('Avg Squared Velocity All Cells');
        xlabel ('Frames');
        ylabel ('Velocity^2 (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end

        % Save the figures in fig, eps and tif format        
        hgsave (h_fig,[SaveDir filesep [imageName '_avgVelocityAllCells.fig']]);
        % print (h_fig, [SaveDir filesep [imageName '_avgVelocityAllCells.eps']],'-depsc2','-tiff');
        % print (h_fig, [SaveDir filesep [imageName '_avgVelocityAllCells.tif']],'-dtiff');      

        % Generate the figure and title
        h_fig2 = figure('Name', imageName);

        % Draw a subplot showing the avg velocity of a single cell    
        ymax = max (avgSingleVelocity) + 1;
        subplot (2,1,1); plot (xAxis, avgSingleVelocity);
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgSingleVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgSingleVelocity, 1, [2 1 1], drugTimepoint);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveSingleVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesSingleVelocity.csv'], est);
           end
        end  
        
        title ('Avg Single Cell Velocity');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end   

        % Draw a subplot showing the avg velocity of a cluster
        ymax = max (avgClusteredVelocity) + 1;
        subplot (2,1,2); plot (xAxis, avgClusteredVelocity); 

        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgClusteredVelocity, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgClusteredVelocity, 1, [2 1 2], drugTimepoint);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveClusteredVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesClusteredVelocity.csv'], est);
           end
        end  
        
        title ('Avg Clustered Cell Velocity');
        xlabel ('Frames');
        ylabel ('Velocity (um/min)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
                
        % Save the figures in fig, eps and tif format     
        hgsave (h_fig2,[SaveDir filesep [imageName '_avgSingleAndClusterVelocity.fig']]);
        % print (h_fig2, [SaveDir filesep [imageName '_avgSingleAndClusterVelocity.eps']],'-depsc2','-tiff');
        % print (h_fig2, [SaveDir filesep [imageName '_avgSingleAndClusterVelocity.tif']],'-dtiff');
        
        % Generate the figure and title
        h_fig3 = figure('Name', imageName);
        
        % Draw a subplot showing the avg velocity of a single cell    
        ymax = max (avgPersistence) + 1;
        subplot (3,1,1); plot (xAxis, avgPersistence);
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgPersistence, 'r'); hold off;
        end
        
        if radioButtons.plotestimate
           hold on;
           [yPlot, est] = ptPlotEstimate (xAxis, avgPersistence, 1, [2 1 1], drugTimepoint);
           hold off;
           
           if ~radioButtons.donotshowplots 
               csvwrite ([SaveDir filesep imageName '_fittedCurveAllVelocity.csv'], [xAxis ; yPlot]);
               csvwrite ([SaveDir filesep imageName '_curveEstimatesAllVelocity.csv'], est);
           end
        end
        
        title ('Avg Persistency All Cells');
        xlabel ('Frames');
        ylabel ('Persistency)');
        if length (xAxis) > 1
           axis ([xAxis(1) xAxis(end) 0 ymax]);
        else
           axis ([xAxis(1) xAxis(1)+1 0 ymax]);
        end
        
        ymax = max (avgPersistence) + 1;
        subplot (3,1,2); plot (xAxis, avgPersistence); 
        
        if radioButtons.runningaverage
            hold on; plot (xAxis, raAvgPersistence, 'r'); hold off;
        end
        
        % Save the figures in fig, eps and tif format        
        hgsave (h_fig3,[SaveDir filesep [imageName '_avgPersistencyAllCells.fig']]);
        
        
        % JvR: create buttons for minimum and max velocity
        [dum, timer]=size(finalPersistence);
        %initialize
        MaxSpeed=[];
        for i=1 : 1 : timer
           MaxSpeed(i,1)= max(finalPersistence{1,i}.Speed);
           dum=isnan(finalPersistence{1,i}.Speed);
           dum= dum-1;
           dum=dum.^2;
           dumNum=find(dum==1);
           avgSpeed(i)=sum(finalPersistence{1,i}.Speed(dumNum))/sum(dum);
        end
        maxVelocity=max(MaxSpeed);
        avgVelo=sum(avgSpeed)/timer;
        %make sliders and textboxen and a OK button for post-analysis of
        %the persistency
        SliderMin=uicontrol('Tag','SliderMin','Style', 'slider', 'Max', maxVelocity, 'Min', 0.1, 'Value', 0.1, 'Units', 'normalized', 'Position', [0.05 0.02 0.13 0.035], 'Callback',{@Slider1_callback});
        Text1=uicontrol('Tag','TxtMin','Style', 'text', 'String', num2str(0.1), 'Units', 'normalized','Position', [0.18 0.02 0.173 0.035]);
        SliderMax=uicontrol('Tag','SliderMax','Style', 'slider', 'Max', maxVelocity, 'Min', 0, 'Value', avgVelo, 'Units', 'normalized', 'Position', [0.38 0.02 0.13 0.035], 'Callback',{@Slider2_callback});
        Text2=uicontrol('Tag','TxtMax','Style', 'text', 'String', num2str(avgVelo), 'Units', 'normalized','Position', [0.58 0.02 0.173 0.035]);
        cmdGo=uicontrol('Tag','cmdGo','Style', 'pushbutton', 'String', 'Go', 'Units', 'normalized','Position', [0.78 0.02 0.173 0.035],'Callback',{@cmdGo_callback});
        
    end  % if ~radioButtons.donotshowplots
        
    % Save MAT files for avg all, single and clustered cell velocity
    cd (SaveDir);
    save ([imageName '_avgCellVelocity.mat'],'avgVelocity');
    save ([imageName '_avgSingleCellVelocity.mat'],'avgSingleVelocity');
    save ([imageName '_avgClusteredCellVelocity.mat'],'avgClusteredVelocity');
    save ([imageName '_avgPersistence.mat'],'avgPersistence');
    save ([imageName '_avgSingleCellVelocity.mat'],'avgSinglePersistence');

    % Save CSV files for avg all, single and clustered cell velocity
    csvwrite ([imageName '_avgCellVelocity.csv'], [xAxis ; avgVelocity]);
    csvwrite ([imageName '_avgCellPersistence.csv'], [xAxis ; avgPersistence]);
    csvwrite ([imageName '_avgSingleCellVelocity.csv'], [xAxis ; avgSingleVelocity]);
    csvwrite ([imageName '_avgSingleCellPersistence.csv'], [xAxis ; avgSinglePersistence]);
    csvwrite ([imageName '_avgClusteredCellVelocity.csv'], [xAxis ; avgClusteredVelocity]);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
     % This is a function when the slider is pressed for SliderMin
     function Slider1_callback(src, eventdata) 
    % JvR Oct 2005
    % save data in handles. Because we are in a figure and not in a GUI we use guihandles instead of guidata 
        handles = guihandles(src);
        % Read the slider
        SliderValue=get(handles.SliderMin, 'Value');
        SliderValue=num2str(SliderValue);
        %Fill in the txtbox
        set(handles.TxtMin, 'String',SliderValue)
%--------------------------------------------------------------------------

     % This is a function when the slider is pressed for SliderMin
     function Slider2_callback(src, eventdata) 
    % JvR Oct 2005
    % save data in handles. Because we are in a figure and not in a GUI we use guihandles instead of guidata 
        handles = guihandles(src);
        % Read the slider
        SliderValue=get(handles.SliderMax, 'Value');
        SliderValue=num2str(SliderValue);
        %Fill in the txtbox
        set(handles.TxtMax, 'String',SliderValue)
%--------------------------------------------------------------------------

     % This is a function when the slider is pressed for SliderMin
     function cmdGo_callback(src, eventdata) 
     % JvR Oct 2005
     % save data in handles. Because we are in a figure and not in a GUI we use guihandles instead of guidata 
     
     %get the global data set in ptCalculateSpeedValues
     global finalPersistence
     global xAxis
     global jrSaveDir
     global jrimageName
     
     % get guihandles
     handles = guihandles(src);
     SliderValueMin=get(handles.SliderMin, 'Value');
     SliderValueMax=get(handles.SliderMax, 'Value');
     [dum, time]=size(finalPersistence); 
     for i=1 : 1 : time
       %find the NaN (division by zeros and leave them out of the
       %calculation
       dumA=isnan(finalPersistence{1,i}.Speed);
       dumA=(dumA-1);
       dumA=dumA.^2;
       dumSpeed=finalPersistence{1,i}.Speed.*dumA;
       dum=find(dumSpeed< SliderValueMax & dumSpeed>SliderValueMin);
       avgPersistenceWithCriterea(i,1)= sum((finalPersistence{1,i}.Persis(dum))) / length(finalPersistence{1,i}.Persis(dum)); 
     end


     % Draw a subplot showing the avg velocity of a single cell
     ymax = max (avgPersistenceWithCriterea) + 1;
     subplot (3,1,3); plot ( avgPersistenceWithCriterea);
     cd (jrSaveDir);
    % Save the figures in fig, eps and tif format        
  %  hgsave (h_fig3,[jrSaveDir filesep [jrimageName '_avgPersistencyAllCells.fig']]);
    
    % Save MAT files for avg all, single and clustered cell velocity
    
  %  save ([imageName '_avgPersistence.mat'],'avgPersistence');
   
    % Save CSV files for avg all, single and clustered cell velocity
    csvwrite ([jrSaveDir '\New_avgSingleCellPersistence.csv'], [xAxis ;  avgPersistenceWithCriterea']);
