function varargout = calcMTDynamics(inputData,saveDirResults,verbose,saveDirInput,saveDirInputType,analysisType)
%calculate MT dynamics: wrapper to load data and specify jobs for calcMTDynamics_main
%
%input(future): inputData(opt): structure with a field anaDat for every anaDat
%               to be analyzed. allowing job-files as input will follow
%               saveDirResults(opt): directory where the ouptut-data should be saved 
%                             (0 if it should not be saved)
%               verbose(opt): whether to display anything (1,0)
%               saveDirInput(opt): directory where the file-list is saved
%                                   that is used for generating the input
%                                   ({0} if no saving, 1 if it should ask)
%                                  fileList: textfile containing the full
%                                  path of all the files: CMTD-files-*.txt
%                                  dataFiles: CMTD-data-*.mat
%                                  
%               saveDirInputType(opt): 1 if fileList, 2 if data, 3 if both
%               analysisType: 1 for normal, 2 for convergence, 3 for both
%
%
%mtdDat(1:end-1) : data for each movie analyzed
%
%fields:  .data  .distance [distance, sigma] in mum
%                .time [time, sigma] in seconds
%
%               stateList contains all the data for analysis (point might not 
%                         necessarily be needed outside the program)
%                         NOTE: all slopes are given in microns/minute!
%
%         .stateList   .point [#of point, state (1/-1/0), slope, slopeSigma]
%                      .single [#of group (binary), state, startData#, endData#, deltaT, deltaD, slope, slopeSigma, significance]
%             (currently not used) .group = [counter, state(1/2/-1/-2/0), startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2...
%                                           avg. time in undetermined state, sum of single counters, significance1, significance2]
%
%               plotData contains all the data necessary to generate the plots
%               (currently not saved)
%
%         .plotData  fields .point, .single, .group1, .group2 with subfields
%                       xData1,xData2,xData3, yData1,yData2,yData3
%
%         .info .dataProperties
%
%         .statistics   cell array, structured {name}{mean value}{std}
%                       distance;
%                       catFreq, resFreq, growthSpeed, shrinkageSpeed,
%                       transition rates, avg distance travelled %time
%                       undetermined
%                       for group only
%
%mtdDat(end) .statistics: same statistics, but for all the analyzed movies
%
%
% 06-03-jonas


%-------------check input
errorMsg = [];

if nargin == 0 | isempty(inputData)
    askLoadData = 1;
else
    askLoadData = 0;
    if ~isfield(inputData,'anaDat')
        errorMsg = [errorMsg, '   there is no field ''anaDat'' in the input structure!',char(10)];
    end
end

if nargin < 2 | isempty(saveDirResults)
    asksaveDirResults = 1;
else
    asksaveDirResults = 0;
    if saveDirResults~=0
        if ~isdir(saveDirResults)
            errorMsg = [errorMsg, '   save directory does not exist',char(10)];
        end
    end
end

if nargin < 3 | isempty(verbose)
    verbose = 1;
else
    if verbose~=1 & verbose~=0
        errorMsg = [errorMsg, '   verbose has to be either 1 or 0!',char(10)];
    end
end

if nargin < 4 | isempty(saveDirInput)
    asksaveDirInput = 0; 
    saveDirInput = 0;
else %if either 0 or 1, asksaveDirInput is adjusted, otherwise it has to be a directory
    if saveDirInput == 0
        asksaveDirInput = 0;
    elseif saveDirInput == 1
        asksaveDirInput = 1;
    elseif ~isdir(saveDirInput)
        errorMsg = [errorMsg, '   save directory does not exist',char(10)];
    end
end

if nargin < 5 | isempty(saveDirInputType)
    saveDirInputType = 1;
elseif ~any(saveDirInputType == [1 2 3])
    errorMsg = [errormsg, '   saveDirInputType has to be 1 to 3 or empty!',char(10)];
end

if nargin < 6 | isempty(analysisType)
    analysisType = 1;
elseif ~any(analysisType == [1 2 3])
    errorMsg = [errormsg, '   analysisType has to be 1 to 3 or empty!',char(10)];
end

if ~isempty(errorMsg)
    error(['error: bad input for calcMTDynamics:',char(10),errorMsg]);
end
%-------------end test input



%-------------define constants

TESTPROBLOW = 0.8;
TESTPROBHIGH = 0.95;

%-------------end define constants



%-------------load data if necessary

if askLoadData
    
    done = 0;
    ct = 1;
    inputData = [];
    
    
    while ~done
        
        %check if default biodata-dir exists and cd if exist
        mainDir = cdBiodata(2);
        
        %get project data file or saved data. filterIndex returns the
        %choice or 0 if cancel
        [fileName,pathName,filterIdx] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files';...
                'CMTD-file-*.txt','project data fileList';...
                'anaDat-*.mat','saved anaDat-structures';...
                '*.mte;*.mts;*.mtx','results files'},'select input file: project data, fileList or saved anaDat');
        
        switch filterIdx 
            %0: nothing
            %1: projectData-File
            %2: fileList
            %3: anaDat
            %4: results-file
            
            case 0
                
                %user does not want to load any further data
                done = 1;
                break; %don't go further in the loop with no data!
                
                
                
            case 1 %standard way
                
                %calculate anaDat
                anaDat = calcAnaDatFromProjectData(pathName,fileName);
                %store anaDat & fileInfo
                inputData(ct).anaDat = anaDat;
                
                %for storage: find identifier and rest of path
                envHome = getenv('HOME');
                envBiodata = getenv('BIODATA');
                envSimdata = getenv('SIMDATA');
                %go through the list until we find a match, else there is a
                %NONE-identifier
                
                if ~isempty(envHome)&strmatch(envHome,pathName)
                    
                    identifier = 'HOME';
                    secondPath = [pathName(length(envHome)+1:end),fileName];
                    
                elseif ~isempty(envBiodata)&strmatch(envBiodata,pathName)
                    
                    identifier = 'BIODATA';
                    secondPath = [pathName(length(envBiodata)+1:end),fileName];
                    
                elseif ~isempty(envSimdata)&strmatch(envSimdata,pathName)
                    
                    identifier = 'SIMDATA';
                    secondPath = [pathName(length(envSimdata)+1:end),fileName];
                    
                else
                    
                    identifier = 'NONE';
                    secondPath = [pathName,fileName];
                    
                end
                
                inputData(ct).fileInfo.identifier = identifier;
                inputData(ct).fileInfo.secondPath = secondPath;
                disp([anaDat(1).info.name,' loaded']);
                
                %count
                ct = ct + 1;
                
                
                
            case {2,4} %fileList from file (pure fileList or results-File)
                %files are stored 'HOME \matlab\chromtrack\shared_work\calcMTDynamics'
                
                %read text with textread (this could surely be nicer, but I don't fully understand textread)
                list = {[]};
                n = 1;
                while isempty(findstr(list{end},'***'))
                list = textread([pathName,fileName],'%s',n,'commentstyle','matlab');
                n = n+1;
                %make sure we're not in an infinite loop
                if n == 10000000
                    error(['no string ''***'' in file to mark end of fileList in ',fileName])
                end
            end
                
                %the list is now a cell array of strings, containing empty cells at the beginning, then an indicator
                %(HOME,BIODATA,SIMDATA) followed by the relative path starting with the filesep
                %now loop through the list and recover the files
                done = 0;
                
                while ~done
                    %drop empty entries
                    if isempty(list{1})
                        list(1) = [];
                    else %read identifier
                        identifier = list{1};
                        list(1) = [];
                        firstPath = getenv(identifier);
                        %check whether identifier exists
                        if isempty(identifier)
                            disp([identifier,' not found - ',list{1},' not read']);
                        else
                            %read second part of path
                            secondPath = list{1};
                            list(1) = [];
                            
                            %find all 'fileseps'
                            filesepList = strfind(secondPath,secondPath(1));
                            
                            %check first filesep: if it is not the correct one,
                            %replace all fileseps
                            if secondPath(1)~=filesep
                                secondPath(filesepList) = filesep;
                            end
                            
                            pathName = [firstPath,secondPath(1:filesepList(end)-1)];
                            fileName = [secondPath(filesepList(end)+1:end)];
                            %calculate anaDat
                            anaDat = calcAnaDatFromProjectData(pathName,fileName);
                            
                            %store anaDat & file info
                            inputData(ct).anaDat = anaDat;
                            inputData(ct).fileInfo.identifier = identifier;
                            inputData(ct).fileInfo.secondPath = secondPath;
                            disp([anaDat(1).info.name,' loaded']);
                            
                            %count
                            ct = ct + 1;
                        end
                    end
                    
                    %stop if there is no list left
                    if isempty(list)|strmatch(list{1},'***')
                        done = 1;
                    end
                end %while ~done
                
                
                
            case 3 %anaDat-file
                
                anaDatFile = load([pathName,fileName]);
                if ~isfield(anaDatFile,'anaDat')
                    error('bad anaDatFile - does not contain field anaDat')
                end
                
                for i = 1:length(anaDatFile)
                    %store anaDat
                    inputData(ct).anaDat = anaDatFile(i).anaDat;
                    if isfield(anaDatFile,'fileInfo')
                        inputData(ct).fileInfo.identifier = anaDatFile(i).fileInfo.identifier;
                        inputData(ct).fileInfo.fileName = anaDatFile(i).fileInfo.secondPath;
                    else
                        
                        inputData(ct).fileInfo.identifier = 'ANADAT';
                        inputData(ct).fileInfo.secondPath = anaDatFile(i).anaDat(1).info.name;
                    end
                    disp([anaDat(1).info.name,' loaded']);
                    
                    %count
                    ct = ct + 1;
                end
                
                
        end %switch filterIdx    
        
        
    end %while ~done
end %if askLoadData
%-----------------end load data

%make sure there is data at all
if isempty(inputData)
    error('no file loaded')
end

%ask for saveDir etc beforehand, so that you can walk away after all the
%input has been given

%-----------------ask for saveDirResults (save mtdDat only!)
if asksaveDirResults & any(bsum2bvec(analysisType)==1)
    %to do: write saveDir/saveName for all results
    %get directory to save results. If cancel, saveDirResults will be numeric 0
    [saveResultsName,saveResultsDir,expOrSim] = uiputfile({'*.mte','experimental MT data';'*.mts','simulation MT data';'*.mtx','any MT data'},'please select filename, directory and filetype to save results file');
elseif saveDirResults ~= 0
    if isdir(saveDirResults)
        %I admire how creative I am late at night figuring out variable names
        %(saveDirResults, saveResultsDir)
        saveResultsName = ['mtDynamics-',nowString];
        saveResultsDir = saveDirResults;
        expOrSim = 3;
    end
else
    saveResultsDir = 0;
    saveResultsName = 0;
    expOrSim = 0;
end
%-----------------end ask for saveDirResults


%-----------------save fileList/data

%first figure out whether we want to save at all
if saveDirInput == 0
    [saveInputFileListName,saveInputFileListDir] = deal(0);
    [saveInputDataName,saveInputDataDir] = deal(0);
else
    saveTypeList = bsum2bvec(saveDirInputType);
    if any(saveTypeList == 1)
        if saveDirInput == 1
            %get directory to save input files. If cancel, saveDir will be numeric 0
            [saveInputFileListName,saveInputFileListDir] = uiputfile('CMTD-file-*','please select directory and filename to save input file list');
        else
            saveInputFileListDir = saveDirInput;
            saveInputFileListName = ['CMTD-file-',nowString];
        end
    else
        [saveInputFileListName,saveInputFileListDir] = deal(0);
    end
    if any(saveTypeList == 2)
        if saveDirInput == 1
            %get directory to save input files. If cancel, saveDir will be numeric 0
            [saveInputDataName,saveInputDataDir] = uiputfile('anaDat-*','please select directory and filename to save input data');
        else
            saveInputDataDir = saveDirInput;
            saveInputDataName = ['anaDat-',nowString];
        end
    else
        [saveInputDataName,saveInputDataDir] = deal(0);
    end
end
%-----------------end save fileList/data

%------write fileList file (so that it is still there if there was an error in the calculations)
if saveInputFileListDir ~= 0
    %go to right file
    cd(saveInputFileListDir);
    
    %createFile with results (mte, mts, mtx if experiment, simulation, unknown)
    fname = [saveInputFileListName,'.txt'];
    fidData = fopen([fname],'w');
    
    %write first a general comment, then all the fileNames (including path)
    
    numData = length(inputData);
    
    %write introduction
    fprintf(fidData,'%s\n%s\n%s\n%s\n','%  MICROTUBULE DYNAMICS ANALYSIS - list of filenames',...
        '%    this file contains a list of filenames that can be used for MT dynamic analysis',...
        '%    the list is made up as: Identifier  restOfPathIncludingFileName, where Identifier is the environment variable for the particular path',...
        '%    Please do not uncomment this header or delete the ''***'' that mark the end of the filenames. Fileseps can be windows or linux type.');
    
    %write file info
    for i = 1:numData
        fprintf(fidData,'%s    %s\n', inputData(i).fileInfo.identifier, inputData(i).fileInfo.secondPath);
    end
    
    %write end-of-fileList marker
    fprintf(fidData,'\n%s\n','***');
    
    %close file
    fprintf(fidData,'file created %s',nowString);
    fclose(fidData);
    
end

%------write anaDat file
if saveInputDataDir ~= 0
    cd(saveInputDataDir);
    save(saveInputDataName,inputData);
end


%============================================

%-----start main program

varargout = calcMTDynamics_main(inputData,verbose,TESTPROBLOW,TESTPROBHIGH,analysisType);

%-----end main program

%===========================================

%-----------write everything to files------------------------------


%-----results file
if expOrSim ~=0
    
    %go to right file
    cd(saveResultsDir);
    
    %save output data as mat-file
    save(['calcMTData-',nowString],'varargout');
    
    %save mtdDat only if it has been calculated
    if any(bsum2bvec(analysisType)==1)
        
        mtdDat = varargout{1};
        
        %createFile with results (mte, mts, mtx if experiment, simulation, unknown)
        switch expOrSim
            case 1
                fname = [saveResultsName,'.mte'];
            case 2
                fname = [saveResultsName,'.mts'];
            case 3
                fname = [saveResultsName,'.mtx'];
        end
        
        fidData = fopen([fname],'w');
        
        
        
        %write first a general comment, then all the fileNames (including path),
        %then the overall statistics, 
        %then the individual ones (general, point, single, group)
        statsAllLength = size(mtdDat(end).statisticsCell,1);
        statsIndLength = size(mtdDat(1).statisticsCell,1);
        numData = length(mtdDat)-1;
        
        %write introduction
        fprintf(fidData,'%s\n%s\n%s\n%s\n','%  MICROTUBULE DYNAMICS ANALYSIS',...
            '%    this file contains first the filenames used to generate the file (can be reused directly!)',...
            '%    next are overall statistics followed by the individual statistics',...
            '%    Please do not uncomment this header or delete the ''***'' that mark the end of the filenames');
        
        %write file info
        for i = 1:numData
            fprintf(fidData,'%s    %s\n', inputData(i).fileInfo.identifier, inputData(i).fileInfo.secondPath);
        end
        
        %write end-of-fileList marker
        fprintf(fidData,'\n%s\n','***');
        
        %write overall statistics
        fprintf(fidData,'\n\n---OVERALL STATISTICS---\n\n');
        
        for i = 1:statsAllLength
            %(I so love these formatted strings)
            fprintf(fidData,'%20s\t%15.6f\t%15.6f',mtdDat(end).statisticsCell{i,:});
            fprintf(fidData,'\n');
        end
        
        %loop through the individual trajectories, print their overall statistics
        fprintf(fidData,'\n\n---individual data I---');
        for i = 1:numData
            fprintf(fidData,'\n\noverall statistics for %s\n',mtdDat(i).info.name);
            for j = 1:statsIndLength
                fprintf(fidData,'%20s\t%15.6f\t%15.6f',mtdDat(i).statisticsCell{j,:});
                fprintf(fidData,'\n');
            end
        end
        
        %loop through individual trajectories, print measurements
        fprintf(fidData,'\n\n---individual data II---');
        for i = 1:numData
            fprintf(fidData,'\n\ndetailed statistics for %s\n',mtdDat(i).info.name);
            
            %point
            numEntries = size(mtdDat(i).stateList.point,1);
            fprintf(fidData,'point data [#of point, state (1/-1/0), slope, slopeSigma]\n');
            for j = 1:numEntries
                fprintf(fidData,'%g\t%g\t%g\t%g\n',mtdDat(i).stateList.point(j,:));
            end
            
            %single
            numEntries = size(mtdDat(i).stateList.single,1);
            fprintf(fidData,'single data [#of group (binary), state, startData#, endData#, deltaT, deltaD, slope, slopeSigma, significance]\n');
            for j = 1:numEntries
                fprintf(fidData,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n*',mtdDat(i).stateList.single(j,:));
            end
            
            %group
%             numEntries = size(mtdDat(i).stateList.group,1);
%             fprintf(fidData,'group data [counter, state(1/2/-1/-2/0), startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2, avg. time in undetermined state, sum of single counters, significance1, significance2]\n');
%             for j = 1:numEntries
%                 fprintf(fidData,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n*',mtdDat(i).stateList.group(j,:));
%             end
        end
        
        %close file
        fprintf(fidData,'file created %s',nowString);
        fclose(fidData);
        
    end %if any(bsum2bvec(analysisType)==1)
end %if expOrSim ~=0
%------end write results file



