function dataSet = runEdgeVelocityAnalysis

clear
clc

label.condition = {'Before Treatment','1-10 min post treatment','10-20 min post treatment'};
label.variable  = {'Protrusion Persistence Time [sec]'
                   'Protrusion Maximum Velocity [nm/sec]'
                   'Protrusion Minimum Velocity [nm/sec]'
                   'Protrusion Average Velocity [nm/sec]'
                   'Protrusion Median Velocity [nm/sec]'
                   'Retraction Persistence Time [sec]'
                   'Retraction Maximum Velocity [nm/sec]'
                   'Retraction Minimum Velocity [nm/sec]'
                   'Retraction Average Velocity [nm/sec]'
                   'Retraction Median Velocity [nm/sec]'};
                 
label.dataSet  = {'Control','Treated'};

label.controlCells = {'Cont1','Cont3','Expt1','Expt2','Expt3'};
label.treatedCells = {'AZD1','AZD2b','AZD3b','AZD3c'};



path = '/home/mv89/files/LCCB/fsm/Michelle/LifeActEDGE/goodControlMovieList.mat';
ML   = MovieList.load(path);

excludeC   = cell(1);
dataSet(1) = callSteps(ML,excludeC,0,label.controlCells);


path = '/home/mv89/files/LCCB/fsm/Michelle/LifeActEDGE/goodTreatedMovieList.mat';
ML   = MovieList.load(path);

excludeT{2} = [1:35 95:108];
dataSet(2)  = callSteps(ML,excludeT,0,label.treatedCells);


plotDataSetsResult(dataSet,label)

end



function dataSet = callSteps(ML,exclude,plotYes,labels)

minLen   = 179;
outLevel = 8;
scale    = true;
missObs  = 10;

formatEdgeVelocity(ML,'minLength',minLen,'outLevel',outLevel,'scale',scale,'missingObs',missObs);

%exclude    = cell(1);


excludeWindowsFromAnalysis(ML,'excBorder',3,'excludeW',exclude);

interval = {1:59,60:100,101:179};
%interval = {1:179};
[cellData,dataSet] = edgeVelocityQuantification(ML,'interval',interval);

if plotYes
    nCell = numel(ML.movies_);
    for iCell = 1:nCell
        
        figure
        set(gcf,'Name',labels{iCell})
        set(gcf,'number','off')
        plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,cellData(iCell).protrusionAnalysis.total.percentage,'LineWidth',2)
        hold on
        plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,cellData(iCell).retractionAnalysis.total.percentage,'r','LineWidth',2)
        plot((interval{1}*ML.movies_{iCell}.timeInterval_)/60,ones(179,1) - (cellData(iCell).protrusionAnalysis.total.percentage + cellData(iCell).retractionAnalysis.total.percentage),'g','LineWidth',2)
        xlabel('Time [min]')
        ylabel('Percentage of the cell Edge')
        title(labels{iCell})
        legend({'Protrusion','Retraction','Quiescent'})
    end
end

end