function [ output_args ] = GCARunMovieList(steps2Run,runMode,projList)
%Quick Personal MovieList Wrapper 20141126 for re-running with the newest code:
% currently the input is just a cell array of project files eventually put in
% movieList format... 

%INPUT: 
% steps2Run: a cell array of strings dictating the steps to perform 
% runMode: Select a predefined selection to run (ignore the first input)
% projList : cell array of directory strings
% createNew: if 1 create new Analysis folder first round if zero do not . 
% Calls GCARun_MakeItSo- which creates or adds to the movie data file
if strcmpi(runMode ,'startAnal');  
   clear steps2Run
   steps2Run{1} = 'thresh1'; 
   steps2Run{2} = 'neuriteOrient';
   steps2Run{3} = 'fixNeuriteOrient'; 
   createNew = 1; 
end 

if strcmpi(runMode,'veilStem')
    clear stepsToRun
    steps2Run{1} = 'veilStem'; 
    createNew = 0; 
end 

if strcmpi(runMode,'filopodia'); 
    clear steps2Run
    steps2Run{1} = 'bodySwitch'; 
    steps2Run{2} = 'prot'; 
    steps2Run{3} = 'reconstruct';
    steps2Run{4} = 'fitFilo'; 
  %  steps2Run{4} = 'fit'; 
    createNew = 0; 
end 

if strcmpi(runMode,'wind')
    clear steps2Run
    steps2Run{1} = 'wind'; 
    createNew = 0;
end 


    
 
    
    



for iProj = 1:numel(projList) 

    dataFolder = projList{iProj}; 
    cd(dataFolder); 
    
    display(['Running Analysis for ' dataFolder '...']); 
    if createNew ==1
        if exist([projList{iProj} filesep 'ANALYSIS'],'dir')==7;
            % rename
            newName  = ([projList{iProj} filesep 'ANALYSIS_Before20150321']);
            movefile([projList{iProj} filesep 'ANALYSIS'], newName);
            display('Old ANALYSIS file detected- Renaming and Copying');
        end
    end
    
    for iStep = 1:numel(steps2Run)
        display(['Step ' steps2Run{iStep}]);
        if (iStep == 1 && createNew ==1 )
            createNewIn =1;
            dataIn = dataFolder;
        else
           dataIn = [dataFolder filesep 'ANALYSIS' ];
            createNewIn = 0;
        end
        % make a new movieData if necessary or load the old movie data and
        % run the analysis...
        GCARunMovie(dataIn,createNewIn,steps2Run{iStep})      
    end
end



end

