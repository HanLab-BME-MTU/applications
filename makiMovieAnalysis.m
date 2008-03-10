function makiMovieAnalysis(serverType,job)
%MAKIMOVIEANALYSIS is the main function to process mammalian kinetochore movies
%
% SYNOPSIS: makiMovieAnalysis
%
% INPUT serverType: 'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW',
%                   'MCAINSH', or 'MADDOX'
%       job: structure containing information about movies to be
%            analyzed. Best set up via a GUI
%       
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn, kjaqaman, gdanuser
% DATE: 27-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFAULTS

%testMovieName = 'cell2.DV';
testMoviePath = 'D:\makiTestData\210607_cell2';
testDataFile = '210607_cell2-makiData-28-Jun-2007-14-38-13.mat';
testMode = true;
if nargin > 1 && ~isempty(job)
    testMode = false;
end


%% MAIN LOOP

% status: rows
% 1: crop
% 2: dataFile
% 3: initCoord
% 4: slist (mmf)
% 5: plane fit
% 6: track
% 7: sister identification
% 8: update frame and kinetochore classification
% 9: frame alignment

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
    job = makiMakeJobPlatformIndependent(job,serverType);
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

            % get the numbers of the tasks to do
            jobs2do = find(status(:,iJob) == -1);
            
            if isfield(job(iJob).dataStruct,'history')
                
                % increase the numRuns index in the history
                job(iJob).dataStruct.history.numRuns = job(iJob).dataStruct.history.numRuns+1;
                
                % set the current version of the data properties for
                % this specific run
                %
                % ASSUMPTION: none of the jobs to be executed
                % below writes new dataProperties -- gd/July-20/2007
                job(iJob).dataStruct.history.dataProperties(job(iJob).dataStruct.history.numRuns) = ...
                    getVersion(job(iJob).dataStruct.dataPropertiesName);
                
            else
                
                % exception handling for old data structures with no
                % history
                error('Job uses old data file format with no history. Run makiUpdateDataFile on this movie');
                
            end
            
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
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'initCoord');
                
                % set the current version of the initCoord
                if success
                    job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns) = ...
                        getVersion(job(iJob).dataStruct.initCoordName);
                else
                    error('unable to secure-save initial coordinates');
                end
                
            else
                
                % write temporarily a zero to history.initCoord
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns)= 0;
                
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %--------------- mixture model fitting ----------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            %I recently moved mmf from task # 7 to task #4, where it should
            %be. Keep an eye on things. -KJ (02-Aug-2007)
            %             if any(jobs2do == 4)
            %
            %                 % write current job to log files
            %                 fprintf(1,'%s : mixture model fitting for %s\n',...
            %                     nowString,job(iJob).dataStruct.projectName);
            %                 fprintf(generalLog,'%s : mixture model fitting\n',nowString);
            %                 fprintf(individualLog,'%s : mixture model fitting\n', nowString);
            %
            %                 % read and prepare data
            %                 rawMovieName = fullfile(job(iJob).dataStruct.rawMoviePath,...
            %                     job(iJob).dataStruct.rawMovieName);
            %
            %                 nTimepoints = job(iJob).dataStruct.dataProperties.movieSize(end);
            %
            %                 % create a 'cord' structure (yes, I know it's a typo)
            %                 cordStruct = makiCoord2Cord(initCoord);
            %
            %                 % loop with movie-chunks
            %                 loopDone = 0;
            %
            %                 % preassign slist
            %                 slist(1:nTimepoints) = ...
            %                     struct('sp',[],...
            %                     'statistics',[],...
            %                     'parms',[],...
            %                     'COM',[]);
            %
            %                 %!!!!!!!!!! update this !!!!!!!!!!!!
            %                 % for testing, set amplitudeCutoff to 6
            %                 % activate lines in makiInitCoord to get good A.C.
            %                 job(iJob).dataStruct.dataProperties.amplitudeCutoff = 6;
            %
            %
            %                 % load first part
            %                 [rawMovie, movieHeader, loadStruct] = ...
            %                     cdLoadMovie({rawMovieName,'raw'},[],job(iJob).dataStruct.dataProperties);
            %
            %                 while ~loopDone
            %
            %                     lf = loadStruct.loadedFrames;
            %                     fprintf(generalLog,sprintf('%s : MMF frames %i:%i\n',nowString,lf(1),lf(end)));
            %                     fprintf(individualLog,sprintf('%s : MMF frames %i:%i\n',nowString,lf(1),lf(end)));
            %                     slist(lf) = ...
            %                         detectSpots_MMF_main(rawMovie,cordStruct(lf),...
            %                         job(iJob).dataStruct.dataProperties,[],sprintf('MMF frames %i:%i/%i',lf(1),lf(end),...
            %                         nTimepoints)); %#ok<AGROW>
            %
            %                     % load more
            %                     if ~isempty(loadStruct.frames2load)
            %                         [rawMovie, movieHeader, loadStruct] = ...
            %                             cdLoadMovie(loadStruct.movieType,[],loadStruct);
            %                     else
            %                         loopDone = 1;
            %                     end
            %
            %                 end % while loop
            %
            %
            %                 % save job
            %                 job(iJob).dataStruct.status(4) = 1;
            %                 job(iJob).dataStruct.statusHelp{4,2} = date;
            %                 job(iJob).dataStruct.slist = slist;
            %                 save(fullfile(job(1).jobPath,job(1).jobName),'job');
            %
            %                 % save dataStruct. Do not overwrite previous slist
            %                 [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'slist');
            %
            %                 % set the current version of the slist
            %                 if success
            %                     job(iJob).dataStruct.history.slist(job(iJob).dataStruct.history.numRuns)= ...
            %                         getVersion(job(iJob).dataStruct.slistName);
            %                     % mixture model fitting depends on initial coordinates -- thus,
            %                     % fill in the history.initCoord field as well
            %                     job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns) = ...
            %                         getVersion(job(iJob).dataStruct.initCoordName);
            %                 else
            %                     error('unable to secure-save mixture model coordinates');
            %                 end
            %
            %             else
            %
            % write temporarily a zero to history.slist
            % Dependent on the jobs to follow, this field might be
            % overwritten again
            job(iJob).dataStruct.history.slist(job(iJob).dataStruct.history.numRuns)= 0;

            %             end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %----------------- plane fitting ----------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            if any(jobs2do == 5)
                
                % write current job to log files
                fprintf(1,'%s : plane fit for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : plane fit\n',nowString);
                fprintf(individualLog,'%s : plane fit\n', nowString);

                % pass and retrieve dataStruct
                job(iJob).dataStruct = makiFitPlane(job(iJob).dataStruct,0);

                % save job
                job(iJob).dataStruct.status(5) = 1;
                job(iJob).dataStruct.statusHelp{5,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                % save dataStruct. Do not overwrite older plane fits
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'planeFit');
                
                % set the current version of the plane fit
                if success
                    job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns) = ...
                        getVersion(job(iJob).dataStruct.planeFitName);
                    % plane fit depends on initial coordinates -- thus,
                    % fill in the history.initCoord field as well
                    job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns) = ...
                        getVersion(job(iJob).dataStruct.initCoordName);
                else
                    error('unable to secure-save plane fits');
                end
                
            else
                
                % write temporarily a zero to history.planeFit
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns)= 0;
                
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %------------------- tracking -------------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            if any(jobs2do == 6)
                
                %write current job to log files
                fprintf(1,'%s : generate tracks for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : generate tracks\n',nowString);
                fprintf(individualLog,'%s : generate tracks\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiGenerateTracks(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(6) = 1;
                job(iJob).dataStruct.statusHelp{6,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older tracks
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'tracks');
                
                % set the current version of the tracks
                if success
                    job(iJob).dataStruct.history.tracks(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.tracksName);
                    % tracking depends on initial coordinates -- thus,
                    % fill in the history.initCoord field as well
                    job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.initCoordName);
                    % if the tracking involves pre-rotation of the
                    % coordinate system into the fitted plane, there is
                    % dependency on planeFit as well -- thus,
                    % fill in the history.planeFit field as well
                    if job(iJob).dataStruct.dataProperties.tracksParam.rotate
                        job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns)= ...
                            getVersion(job(iJob).dataStruct.planeFitName);
                    end;
                else
                    error('unable to secure-save tracks');
                end
                
            else
                
                % write temporarily a zero to history.tracks
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.tracks(job(iJob).dataStruct.history.numRuns) = 0;
                
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %---------------- group sisters -----------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            if any(jobs2do == 7)
                
                %write current job to log files
                fprintf(1,'%s : group sisters for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : group sisters\n',nowString);
                fprintf(individualLog,'%s : group sisters\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiGroupSisters(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(7) = 1;
                job(iJob).dataStruct.statusHelp{7,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older sisterLists
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'sisterList');
                
                % set the current version of the sister list
                if success
                    job(iJob).dataStruct.history.sisterList(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.sisterListName);
                    % sister identification depends on plane fit and tracks -- thus,
                    % fill in the history.tracks and history.planeFit fields as well
                    job(iJob).dataStruct.history.tracks(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.tracksName);
                    job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.planeFitName);
                else
                    error('unable to secure-save sister list');
                end
                
            else
                
                % write temporarily a zero to history.tracks
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.sisterList(job(iJob).dataStruct.history.numRuns) = 0;
                
            end

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %------------- update classification ------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            if any(jobs2do == 8)
                
                %write current job to log files
                fprintf(1,'%s : update classification for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : update classification\n',nowString);
                fprintf(individualLog,'%s : update classification\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiUpdateClass(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(8) = 1;
                job(iJob).dataStruct.statusHelp{8,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older classification
                %updates
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'updatedClass');
                
                % set the current version of the classification update
                if success
                    job(iJob).dataStruct.history.updatedClass(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.updatedClassName);
                    % classification update depends on plane fit, tracks 
                    % and sisterList -- thus fill in those fields as well
                    job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.planeFitName);
                    job(iJob).dataStruct.history.tracks(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.tracksName);
                    job(iJob).dataStruct.history.sisterList(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.sisterListName);
                else
                    error('unable to secure-save updated classification');
                end
                
            else

                % write temporarily a zero to history.updatedClass
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.updatedClass(job(iJob).dataStruct.history.numRuns)= 0;
                
            end

            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %---------------- align frames ------------------------
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            if any(jobs2do == 9)
                
                %write current job to log files
                fprintf(1,'%s : align frames for %s\n',...
                    nowString,job(iJob).dataStruct.projectName);
                fprintf(generalLog,'%s : align frames\n',nowString);
                fprintf(individualLog,'%s : align frames\n', nowString);

                %pass and retrieve dataStruct
                job(iJob).dataStruct = makiAlignFrames(job(iJob).dataStruct);

                %save job
                job(iJob).dataStruct.status(9) = 1;
                job(iJob).dataStruct.statusHelp{9,2} = date;
                save(fullfile(job(1).jobPath,job(1).jobName),'job');
                
                %save dataStruct. Do not overwrite older frame alignments
                [success, job(iJob).dataStruct ] = makiSaveDataFile(serverType,job(iJob).dataStruct,'frameAlignment');
                
                % set the current version of the frame alignment
                if success
                    job(iJob).dataStruct.history.frameAlignment(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.frameAlignmentName);
                    % frame alignment depends on initCoord, plane fit, tracks
                    % and updatedClass -- thus fill in those fields as well
                    job(iJob).dataStruct.history.initCoord(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.initCoordName);
                    job(iJob).dataStruct.history.planeFit(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.planeFitName);
                    job(iJob).dataStruct.history.tracks(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.tracksName);
                    job(iJob).dataStruct.history.updatedClass(job(iJob).dataStruct.history.numRuns)= ...
                        getVersion(job(iJob).dataStruct.updatedClassName);
                else
                    error('unable to secure-save frame alignment');
                end
                
            else

                % write temporarily a zero to history.frameAlignment
                % Dependent on the jobs to follow, this field might be
                % overwritten again
                job(iJob).dataStruct.history.frameAlignment(job(iJob).dataStruct.history.numRuns)= 0;
                
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
    
        %save history in dataStruct to the dataFile a last time
        makiSaveDataFile(serverType,job(iJob).dataStruct);
        
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
%     makiSaveDataFile(serverType,job(iJob).dataStruct);
% 
% end % ------------ filter ---------------


%% LOCAL FUNCTIONS 

function v = getVersion(fName)
% Service function to retrieve the current version from a filename in the
% job data structure
%
% The function assumes that the version index is the last number before
% '.mat' and is preceeded by an '_'; 
%
% for example: 'gaga_cpi11_anythingelse_23.mat' 
% will generate a numerical value 23
v = str2num(fName(regexp(fName,'_\d+\.mat')+1:regexp(fName,'.mat')-1));



