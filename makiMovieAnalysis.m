function makiMovieAnalysis(job)
%MAKIMOVIEANALYSIS is the main function to process mammalian kinetochore movies
%
% SYNOPSIS: makiMovieAnalysis
%
% INPUT job: structure containing information about movies to be
%                analyzed. Best set up via a GUI
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn, modified by: kjaqaman
% DATE: 27-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===============
%% DEFAULTS
%===============

%testMovieName = 'cell2.DV';
testMoviePath = 'D:\makiTestData\210607_cell2';
testDataFile = '210607_cell2-makiData-28-Jun-2007-14-38-13.mat';
testMode = true;
if nargin > 0 && ~isempty(job)
    testMode = false;
end


%===============
%% MAIN LOOP
%===============

% status: rows
% 1: crop
% 2: dataFile
% 3: initCoord
% 4: plane fit
% 5: track
% 6: identify sisters
% 7: slist

if testMode

    % settings - use as defaults for batch definition later
    % sigmaCorrection = 1.5

    % create job
    job = struct('dataStruct',makiLoadDataFile(...
        fullfile(testMoviePath,testDataFile)));

    % set status
    % MMF
    job(1).dataStruct.status(5) = -1;

    job(1).jobPath = 'D:\makiTestData';
    job(1).jobName = sprintf('testJob-%s.mat',nowString);
else
    % revert job paths
    job = makiMakeJobPlatformIndependent(job);
end

% collect status of all the jobs
status = catStruct(2,'job.dataStruct.status');
nJobs = length(job);

% set up logfile - maybe add field logPath
logFileName = fullfile(job(1).jobPath,sprintf('%s.log',job(1).jobName(1:end-4)));
generalLog = fopen(logFileName,'a+');

% loop though all. Allow for multiple passes
done = false;
individualLog = [];

while ~done

    for iJob = 1:nJobs

        try

            % open individual logfile
            logFileName = fullfile(job(iJob).dataStruct.dataFilePath,[...
                job(iJob).dataStruct.projectName,'_analysis.log']);
            individualLog = fopen(logFileName,'a+');
            fprintf(individualLog,'\n\n-----------------------\n');

            % get the numbers of the jobs to do
            jobs2do = find(status(:,iJob) == -1);


            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %--------------- initial coords -----------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if any(jobs2do == 3)

                % write current job to log files
                fprintf(1,'%s : find initial coords for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : find initial coords\n',nowString);
                fprintf(individualLog,'%s : find initial coords\n', nowString);

                % pass and retrieve dataStruct
                job(iJob).dataStruct = makiInitCoord(job(iJob).dataStruct);


                % save job
                job(iJob).dataStruct.status(3) = 1;
                job(iJob).dataStruct.statusHelp{3,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                % save dataStruct. Do not overwrite older initCoord
                makiSaveDataFile(job(iJob).dataStruct,'initCoord');

            end
            
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %----------------- congression ------------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if any(jobs2do == 4)
                % write current job to log files
                fprintf(1,'%s : metaphase plate fit %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : metaphase plate fit\n',nowString);
                fprintf(individualLog,'%s : metaphase plate fit\n', nowString);

                % pass and retrieve dataStruct
                job(iJob).dataStruct = makiFitPlane(job(iJob).dataStruct,0);

                % save job
                job(iJob).dataStruct.status(4) = 1;
                job(iJob).dataStruct.statusHelp{4,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                % save dataStruct. Do not overwrite older plane fits
                makiSaveDataFile(job(iJob).dataStruct,'planeFit');
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %------------------- tracking -------------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if any(jobs2do == 5)
                
                %write current job to log files
                fprintf(1,'%s : generate tracks for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : generate tracks\n',nowString);
                fprintf(individualLog,'%s : generate tracks\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiGenerateTracks(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(5) = 1;
                job(iJob).dataStruct.statusHelp{5,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older tracks
                makiSaveDataFile(job(iJob).dataStruct,'tracks');
                
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %---------------- group sisters -----------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if any(jobs2do == 6)
                
                %write current job to log files
                fprintf(1,'%s : group sisters for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : group sisters\n',nowString);
                fprintf(individualLog,'%s : group sisters\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiGroupSisters(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(6) = 1;
                job(iJob).dataStruct.statusHelp{6,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older sisterLists
                makiSaveDataFile(job(iJob).dataStruct,'sisterList');
                
            end

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %--------------- mixture model fitting ----------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if any(jobs2do == 7)

                % write current job to log files
                fprintf(1,'%s : mixture model fitting for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : mixture model fitting\n',nowString);
                fprintf(individualLog,'%s : mixture model fitting\n', nowString);

                % read and prepare data
                rawMovieName = fullfile(job(iJob).dataStruct.rawMoviePath,...
                    job(iJob).dataStruct.rawMovieName);

                nTimepoints = job(iJob).dataStruct.dataProperties.movieSize(end);

                % create a 'cord' structure (yes, I know it's a typo)
                cordStruct = makiCoord2Cord(initCoord);

                % loop with movie-chunks
                loopDone = 0;

                % preassign slist
                slist(1:nTimepoints) = ...
                    struct('sp',[],...
                    'statistics',[],...
                    'parms',[],...
                    'COM',[]);

                %!!!!!!!!!! update this !!!!!!!!!!!!
                % for testing, set amplitudeCutoff to 6
                % activate lines in makiInitCoord to get good A.C.
                job(iJob).dataStruct.dataProperties.amplitudeCutoff = 6;


                % load first part
                [rawMovie, movieHeader, loadStruct] = ...
                    cdLoadMovie({rawMovieName,'raw'},[],job(iJob).dataStruct.dataProperties);

                while ~loopDone

                    lf = loadStruct.loadedFrames;
                    fprintf(generalLog,sprintf('%s : MMF frames %i:%i\n',nowString,lf(1),lf(end)));
                    fprintf(individualLog,sprintf('%s : MMF frames %i:%i\n',nowString,lf(1),lf(end)));
                    slist(lf) = ...
                        detectSpots_MMF_main(rawMovie,cordStruct(lf),...
                        job(iJob).dataStruct.dataProperties,[],sprintf('MMF frames %i:%i/%i',lf(1),lf(end),...
                        nTimepoints)); %#ok<AGROW>

                    % load more
                    if ~isempty(loadStruct.frames2load)
                        [rawMovie, movieHeader, loadStruct] = ...
                            cdLoadMovie(loadStruct.movieType,[],loadStruct);
                    else
                        loopDone = 1;
                    end

                end % while loop


                % save job
                job(iJob).dataStruct.status(7) = 1;
                job(iJob).dataStruct.statusHelp{7,2} = date;
                job(iJob).dataStruct.slist = slist;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                % save dataStruct. Do not overwrite previous slist
                makiSaveDataFile(job(iJob).dataStruct,'slist');

            end

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %--------------- THE END (of job-loop) ------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        catch
            % display and log error
            err = lasterror;

            % print to screen
            fprintf(1,'%s : job %i - error: %s\n',nowString,iJob,err.message);
            for iErr = 1:length(err.stack)
                fprintf(1, 'in %s at %i\n',...
                    err.stack(iErr).name,err.stack(iErr).line);
            end

            % print to individual log
            if ~isempty(individualLog)
                fprintf(individualLog,...
                    '%s : job %i - error: %s\n',nowString,iJob,err.message);
                for iErr = 1:length(err.stack)
                    fprintf(individualLog, 'in %s at %i\n',...
                        err.stack(iErr).name,err.stack(iErr).line);
                end
            end

            % print to global log
            fprintf(generalLog,...
                '%s : job %i - error: %s\n',nowString,iJob,err.message);
            for iErr = 1:length(err.stack)
                fprintf(generalLog, 'in %s at %i\n',...
                    err.stack(iErr).name,err.stack(iErr).line);
            end

        end
        try
            fprintf(generalLog,'%s : evaluation finished\n',nowString);
            fprintf(individualLog,'%s : evaluation finished\n',nowString);
            fprintf(1,'%s : evaluation finished\n',nowString);
            fclose(individualLog);
        catch
        end
        individualLog = [];
    
        %save history in dataStruct
        
    end % loop jobs

    % there is no redo for now
    done = true;


end % while ~ done

% close general log
try
    fclose(all);
catch
end



% removed step(s)

% %++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %--------------- filter movie -------------------------
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++
% if any(jobs2do == 3)
% 
%     % write current job to log files
%     fprintf(1,'%s : filtering movie %s\n',nowString,...
%         job(iJob).dataStruct.rawMovieName);
%     fprintf(generalLog,'%s : filtering movie\n',nowString);
%     fprintf(individualLog,'%s : filtering movie\n', nowString);
% 
%     filteredMovieName = fullfile(job(iJob).dataStruct.dataFilePath,...
%         job(iJob).dataStruct.filteredMovieName);
%     rawMovieName = fullfile(job(iJob).dataStruct.rawMoviePath,...
%         job(iJob).dataStruct.rawMovieName);
%     dataProperties = job(iJob).dataStruct.dataProperties;
% 
%     %check wheter another movie already
%     %exists (within the loading loop, we
%     %want to append!)
%     if exist(filteredMovieName,'file')
%         fprintf(generalLog,'%s : delete(%s);\n',nowString, filteredMovieName);
%         fprintf(individualLog,'%s : delete(%s);\n',nowString, filteredMovieName);
%         delete(filteredMovieName);
%     end
% 
%     % loop with movie-chunks
%     loopDone = 0;
% 
%     % load first part
%     [rawMovie, movieHeader, loadStruct] = ...
%         cdLoadMovie({rawMovieName,'raw'},[],dataProperties);
% 
%     while ~loopDone
% 
%         %filter movie
%         lf = loadStruct.loadedFrames;
%         fprintf(generalLog,sprintf('%s : filtermovie frames %i:%i\n',nowString,lf(1),lf(end)));
%         fprintf(individualLog,sprintf('%s : filtermovie frames %i:%i\n',nowString,lf(1),lf(end)));
%         filteredMovie = filtermovie(...
%             rawMovie,dataProperties.FILTERPRM,...
%             sprintf('filter frames %i:%i/%i',lf(1),lf(end),...
%             dataProperties.movieSize(end)));
% 
%         %save. Writemat appends to an existing file
%         writemat(filteredMovieName,filteredMovie,1,5);
% 
%         % load more
%         if ~isempty(loadStruct.frames2load)
%             [rawMovie, movieHeader, loadStruct] = ...
%                 cdLoadMovie(loadStruct.movieType,[],loadStruct);
%         else
%             loopDone = 1;
%         end
% 
%     end % while loop
% 
%     clear('filteredMovie','rawMovie'); %to prevent memory problems
% 
%     % save job
%     job(iJob).dataStruct.status(3) = 1;
%     job(iJob).dataStruct.statusHelp{3,2} = date;
%     save(fullfile(job(1).jobPath,job(1).jobName),'job');
%     % save dataStruct
%     makiSaveDataFile(job(iJob).dataStruct);
% 
% end % ------------ filter ---------------


