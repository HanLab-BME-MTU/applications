function experiment = yeastAnalyzeMovies(forcedStatus,resetCutoff,stopEval)
%YEASTANALYZEMOVIES is a batch processing tool for analyzing yeast movies
%
% SYNOPSIS: yeastAnalyzeMovies
%
% INPUT none
%
% OUTPUT none
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: jdorn
% DATE: 24-Feb-2006
%
%

% forcedStatus forces the status to a certain level
if nargin < 1 || isempty(forcedStatus)
    forcedStatus = 0;
end
% resetCutoff sets cutoff to 0, which forces new selection of amplitude
% cutoff
if nargin < 2 || isempty(resetCutoff)
    resetCutoff = 0;
end

% load directory
done = 0;
nExperiments = 0;
while ~done
    selectedDir = uigetdir(cdBiodata);
    if selectedDir == 0
        done = 1;
    else
        nExperiments = nExperiments + 1;
        dirName{nExperiments,1} = selectedDir;
    end
end

% check whether anything was loaded
if nExperiments == 0
    experiment = [];
    return
end

experiment(1:nExperiments) = struct('rawMovieName',[],'projectName',[],...
    'movieDir',[],'filteredMovieName',[],'dataPropertiesName',[],...
    'dataProperties',[],'movieHeaderName',[],'movieHeader',[],...
    'slistName',[],'slist',[],'idlistName',[],'idlist',[],...
    'idlistLName',[],'idlist_L',[],'status',0,'testRatiosName',[],...
    'testRatios',[]);

for iExperiment = 1:nExperiments

    % look for movie
    rawMovieName = searchFiles('r3d$|3D.dv$','DIC',dirName{iExperiment},0);

    % check that we have found two movies
    if size(rawMovieName,1) == 1
        % good
    else
        error('not 1 movie found');
    end

    experiment(iExperiment).rawMovieName = rawMovieName{1};
    experiment(iExperiment).movieDir = rawMovieName{2};


    % find the movieName (without _R3D, _D3D etc)
    endOfNameIdx = regexpi(experiment(iExperiment).rawMovieName,'(_.3d)*\.')-1;

    experiment(iExperiment).projectName = ...
        experiment(iExperiment).rawMovieName(1:endOfNameIdx);

    % filtered movie
    filteredMovieName = searchFiles('fim','',experiment(iExperiment).movieDir);
    movieHeaderName = searchFiles('ovieHeader','',experiment(iExperiment).movieDir);



    % create names for
    % - dataProperties
    % - movieHeader
    % - filteredMovie
    % - slist
    % - idlist
    experiment(iExperiment).dataPropertiesName = ...
        ['dataProperties_',experiment(iExperiment).projectName,'.mat'];
    experiment(iExperiment).movieHeaderName = ...
        movieHeaderName{1};
    experiment(iExperiment).filteredMovieName = ...
        filteredMovieName{1};
    experiment(iExperiment).slistName = ...
        ['slist_',experiment(iExperiment).projectName,'.mat'];
    experiment(iExperiment).idlistName = ...
        ['idlist_',experiment(iExperiment).projectName,'.mat'];
    experiment(iExperiment).idlist_LName = ...
        ['idlist_L_',experiment(iExperiment).projectName,'.mat'];
    experiment(iExperiment).testRatioName = ...
        ['testRatios_',experiment(iExperiment).projectName,'.mat'];



    % load files, determine status
    % 0: nothing done yet
    % 1: dataProperties, movieHeader have been saved
    % 2: filtered movie
    % 3: slist
    % 4: idlist
    % 5: idlist_L

    % change directory to movieDir
    oldDir = cd(experiment(iExperiment).movieDir);

    % check for dataProperties. If they exist, load dataFile. Otherwise, go
    % the standard route.
    if ~exist(experiment(iExperiment).dataPropertiesName,'file')
        % load from dataFile
        [dummy,dataProperties] = ...
            loadProjectData([],pwd,'last');
        if ~isfield(dataProperties,'maxSize')
            dataProperties.maxSize = 200000000;
        end
        if ~isfield(dataProperties,'amplitudeCutoff')
            dataProperties.amplitudeCutoff = 0;
        end
        if dataProperties.MAXSPOTS == 5
            answer = inputdlg(...
                sprintf('enter maxSpots for %s',...
                experiment(iExperiment).projectName),...
                'bla',1,{'2'});
            dataProperties.MAXSPOTS = str2double(answer{1});
        end
        experiment(iExperiment).dataProperties = dataProperties;
        experiment(iExperiment).status = 2;

        % write out movieHeader, dataProperties
        movieHeader = readr3dheader(experiment(iExperiment).rawMovieName);
        % correct for the mistake of reading (and saving) the header before
        % which removed all the correction info. In the future, just load
        % movieHeader, save with name "movieHeader"
        if exist('correctionData.mat','file')
            load correctionData
            movieHeader.correctInfo = correctionData.info;
            movieHeader.numTimepoints = ...
                movieHeader.numTimepoints - ...
                sum(correctionData.info.correctFrames);
        end
        save(experiment(iExperiment).movieHeaderName,'movieHeader');
        save(experiment(iExperiment).dataPropertiesName,'dataProperties')
        experiment(iExperiment).movieHeader = movieHeader;

    else


        % status 1
        if exist(experiment(iExperiment).dataPropertiesName,'file')
            load(experiment(iExperiment).dataPropertiesName)
            experiment(iExperiment).dataProperties = dataProperties;
            load(experiment(iExperiment).movieHeaderName)
            if exist('r3dMovieHeader','var')
                experiment(iExperiment).movieHeader = r3dMovieHeader;
            else
                experiment(iExperiment).movieHeader = movieHeader;
            end
            experiment(iExperiment).status = 1;

            % status 2
            if exist(experiment(iExperiment).filteredMovieName,'file')
                % don't load filtered movie!
                experiment(iExperiment).status = 2;

                % status 3
                if exist(experiment(iExperiment).slistName,'file')
                    load(experiment(iExperiment).slistName)
                    experiment(iExperiment).slist = slist;
                    experiment(iExperiment).status = 3;

                    % status 4
                    if exist(experiment(iExperiment).idlistName,'file')
                        load(experiment(iExperiment).idlistName)
                        experiment(iExperiment).idlist = idlist;
                        experiment(iExperiment).status = 4;

                        % status 5
                        if exist(experiment(iExperiment).idlist_LName,'file')
                            load(experiment(iExperiment).idlist_LName)
                            experiment(iExperiment).idlist_L = idlist_L;
                            experiment(iExperiment).status = 5;
                        end
                    end
                end
            end
        else

            error('undefined loading from raw movie')
            %         % create movie header, dataProperties
            %
            %         % read movieHeader
            %         movieHeader = readr3dheader(...
            %             fullfile(experiment(iExperiment).movieDir,experiment(iExperiment).rawMovieName));
            %
            %         % create dataProperties
            %         dataProperties = defaultDataProperties(movieHeader);
            %         % set maxSpots
            %         dataProperties.MAXSPOTS = movieProperties.MAXSPOTS;
            %         dataProperties.NA = movieProperties.NA;
            %         dataProperties.maxSize = movieProperties.maxSize;
            %         dataProperties.name = ...
            %             experiment(iExperiment).projectName;
            %         dataProperties.amplitudeCutoff = movieProperties.amplitudeCutoff;
            %
            %         % find filter parameters
            %         [FT_XY, FT_Z] = calcFilterParms(...
            %             dataProperties.WVL,dataProperties.NA,1.51,'gauss',...
            %             [1 1], ...
            %             [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
            %         patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
            %         dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
            %         dataProperties.FT_SIGMA = [FT_XY,FT_XY,FT_Z];
            %
            %         % store dataProperties, movieHeader
            %         experiment(iExperiment).dataProperties = dataProperties;
            %         experiment(iExperiment).movieHeader = movieHeader;
            %
            %         % save to disk
            %         save(fullfile(experiment(iExperiment).movieDir,...
            %             experiment(iExperiment).dataPropertiesName),...
            %             'dataProperties');
            %         save(fullfile(experiment(iExperiment).movieDir,...
            %             experiment(iExperiment).movieHeaderName),...
            %             'movieHeader');
            %
            %         % update status
            %         experiment(iExperiment).status = 1;
            %
            %         % close imaris - don't clear movieProperties!
            %         clear imarisHandle dataProperties movieHeader
            %
            %

        end % if exist...



    end % if loaded from projectData
    
    
    % force status
    if forcedStatus
        experiment(iExperiment).status = min(forcedStatus,...
            experiment(iExperiment).status);
    end
    
    % reset cutoff
    if resetCutoff
        % reset cutoff
        experiment(iExperiment).dataProperties.amplitudeCutoff = 0;
        % force status
        experiment(iExperiment).status = min(2,...
            experiment(iExperiment).status);
    end
    
    
    % go back to oldDir
    cd(oldDir)
end % loop experiment

% loop to get all status to the level that we need to launch labelgui
for iExperiment = 1:nExperiments
    while experiment(iExperiment).status < 4

        try

            switch experiment(iExperiment).status

                case 1

                    %---------------
                    % filter movie
                    %---------------

                    % talk to user
                    evaluationStartTime = clock;
                    disp(sprintf(...
                        '\n\nexperiment %i: filtering movie ''%s''. Start at %s',...
                        iExperiment,...
                        experiment(iExperiment).projectName,...
                        datestr(evaluationStartTime)));
                    
                    % correct background here!
                    %
                    % - only possible with cdEditDatafile

                    % load movie
                    [rawMovie,dummy,loadStruct] = ...
                        cdLoadMovie('corr/raw',experiment(iExperiment).movieDir,...
                        experiment(iExperiment).dataProperties);

                    % loop movieChunks
                    % filter the movie
                    done = 0;
                    while ~done
                        % filter movie
                        filteredMovie = filtermovie(...
                            rawMovie,...
                            experiment(iExperiment).dataProperties.FILTERPRM);
                        % write movie to file
                        writemat(fullfile(...
                            experiment(iExperiment).movieDir,...
                            experiment(iExperiment).filteredMovieName),...
                            filteredMovie,1,5);
                        % load next bit of movie
                        if isempty(loadStruct.frames2load)
                            rawMovie = [];
                            done = 1;
                        else
                            rawMovie = [];
                            [rawMovie,dummy,loadStruct] = ...
                                cdLoadMovie('raw',experiment(iExperiment).movieDir,...
                                loadStruct);
                        end
                    end % while-loop

                    % update status
                    experiment(iExperiment).status = 2;

                    % display elapsed time
                    disp(sprintf('time elapsed: %i s',...
                        round(etime(clock,evaluationStartTime))))



                case 2

                    %----------------
                    % detect spots
                    %----------------

                    % talk to user
                    evaluationStartTime = clock;
                    disp(sprintf(...
                        '\n\nexperiment %i: detecting spots for ''%s''. Start at %s',...
                        iExperiment,...
                        experiment(iExperiment).projectName,...
                        datestr(evaluationStartTime)));


                    % pass movieNames to detector
                    [slist, dataProperties, testRatios, debugData] = detectSpots(...
                        fullfile(...
                        experiment(iExperiment).movieDir,...
                        experiment(iExperiment).rawMovieName), ...
                        fullfile(...
                        experiment(iExperiment).movieDir,...
                        experiment(iExperiment).filteredMovieName), ...
                        experiment(iExperiment).dataProperties,2,...
                        struct('debug',[1,2]));

                    % store debugData
%                     experiment(iExperiment).detectorFStats = ...
%                         debugData.fStats;
%                     experiment(iExperiment).residualImages = ...
%                         debugData.residualImages;

                    % store slist and dataProperties
                    experiment(iExperiment).slist = slist;
                    experiment(iExperiment).dataProperties = dataProperties;


                    % write slist to file
                    save(fullfile(experiment(iExperiment).movieDir,...
                        experiment(iExperiment).slistName),'slist');
                    save(fullfile(experiment(iExperiment).movieDir,...
                        experiment(iExperiment).dataPropertiesName),'dataProperties')
                    save(fullfile(experiment(iExperiment).movieDir,...
                        experiment(iExperiment).testRatioName),'testRatios');

                    % update status
                    experiment(iExperiment).status = 3;

                    % display elapsed time
                    disp(sprintf('time elapsed: %i s',...
                        round(etime(clock,evaluationStartTime))))

                case 3

                    %-----------
                    % link tags
                    %-----------

                    % talk to user
                    evaluationStartTime = clock;
                    disp(sprintf(...
                        '\n\nexperiment %i: linking tags for ''%s''. Start at %s',...
                        iExperiment,...
                        experiment(iExperiment).projectName,...
                        datestr(evaluationStartTime)));

                    idlist = linker(...
                        experiment(iExperiment).slist,...
                        experiment(iExperiment).dataProperties,1);

                    % label idlist
                    nTags = length(idlist(1).stats.labelcolor);
                    labels = cellstr(num2str([1:nTags]'));
                    idlist(1).stats.labelcolor = labels;
                    idlist(1).stats.labellist = labels;

                    % store idlist
                    experiment(iExperiment).idlist = idlist;

                    % write idlist to file
                    save(fullfile(experiment(iExperiment).movieDir,...
                        experiment(iExperiment).idlistName),'idlist');

                    % update status
                    experiment(iExperiment).status = 4;

                    % display elapsed time
                    disp(sprintf('time elapsed: %i s',...
                        round(etime(clock,evaluationStartTime))))
                case 99
                    % failed movie

                otherwise
                    error('unknown status %i',...
                        experiment(iExperiment).status)
            end % switch status
        catch
            err = lasterror;
            experiment(iExperiment) = ...
                handleError(err,...
                experiment(iExperiment),...
                evaluationStartTime);
        end
    end % loop till movie is analyzed
end % loop till all experiments are analyzed



%====================
% SHOW RESULT
%====================
labelguiH = [];
for iExperiment = 1:nExperiments
    if experiment(iExperiment).status == 99  ||...
            experiment(iExperiment).status > 4
        % don't go further
    else
        try
            % talk to user
            evaluationStartTime = clock;
            disp(sprintf(...
                '\n\nexperiment %i: launch labelgui for ''%s''. Start at %s',...
                iExperiment,...
                experiment(iExperiment).projectName,...
                datestr(evaluationStartTime)));
            % load filtered movie
            [filteredMovie, dummy, loadStruct] = ...
                cdLoadMovie({fullfile(experiment(iExperiment).movieDir,...
                experiment(iExperiment).filteredMovieName),'filtered'}, [], ...
                experiment(iExperiment).dataProperties);

            labelguiH = LG_loadAllFromOutside(...
                filteredMovie,...
                experiment(iExperiment).movieDir,...
                loadStruct,...
                experiment(iExperiment).dataProperties,...
                experiment(iExperiment).idlist,'idlist');
            uiwait(labelguiH);
            % read idlist_L
            idlist_L = LG_readIdlistFromOutside;

            if isempty(idlist_L)
                error('checking idlist aborted')
            end

            % store and save
            experiment(iExperiment).idlist_L = idlist_L;
            save(fullfile(experiment(iExperiment).movieDir,...
                experiment(iExperiment).idlist_LName),'idlist_L');

            % update status
            experiment(iExperiment).status = 5;
            % display elapsed time
            disp(sprintf('time elapsed: %i s',...
                round(etime(clock,evaluationStartTime))))

        catch
            err = lasterror;
            experiment(iExperiment) = ...
                handleError(err,...
                experiment(iExperiment),...
                evaluationStartTime);
        end
    end
end % loop experiments
if ishandle(labelguiH)
    close(labelguiH);
end


% ---- track here ------
for iExperiment = 1:nExperiments
    if experiment(iExperiment).status == 99  ||...
            experiment(iExperiment).status > 6
        % don't track
    else
        try
            % talk to user
            evaluationStartTime = clock;
            disp(sprintf(...
                '\n\nexperiment %i: track tags for ''%s''. Start at %s',...
                iExperiment,...
                experiment(iExperiment).projectName,...
                datestr(evaluationStartTime)));

            % track tags
            [idlisttrack,debugData] = ...
                tagTracker({ fullfile(...
                experiment(iExperiment).movieDir,...
                experiment(iExperiment).rawMovieName),...
                'corr/raw'},experiment(iExperiment).idlist_L,...
                experiment(iExperiment).dataProperties,...
                1,struct('fStats',1));


            % store and save
            experiment(iExperiment).trackerFStats = debugData.fStats;
            experiment(iExperiment).idlisttrack = idlisttrack;
            save(fullfile(experiment(iExperiment).movieDir,...
                experiment(iExperiment).idlist_LName),'idlist_L');

            % update status
            experiment(iExperiment).status = 6;
            % display elapsed time
            disp(sprintf('time elapsed: %i s',...
                round(etime(clock,evaluationStartTime))))

        catch
            err = lasterror;
            experiment(iExperiment) = ...
                handleError(err,...
                experiment(iExperiment),...
                evaluationStartTime);
        end
    end
end % loop experiments


%====================
% SHOW RESULT
%====================
labelguiH = [];
for iExperiment = 1:nExperiments
    if experiment(iExperiment).status == 99  ||...
            experiment(iExperiment).status > 7
        % don't go further
    else
        try
            % talk to user
            evaluationStartTime = clock;
            disp(sprintf(...
                '\n\nexperiment %i: launch labelgui for ''%s''. Start at %s',...
                iExperiment,...
                experiment(iExperiment).projectName,...
                datestr(evaluationStartTime)));
            % load filtered movie
            [filteredMovie, dummy, loadStruct] = ...
                cdLoadMovie({fullfile(experiment(iExperiment).movieDir,...
                experiment(iExperiment).filteredMovieName),'filtered'}, [], ...
                experiment(iExperiment).dataProperties);

            labelguiH = LG_loadAllFromOutside(...
                filteredMovie,...
                experiment(iExperiment).movieDir,...
                loadStruct,...
                experiment(iExperiment).dataProperties,...
                experiment(iExperiment).idlisttrack,'idlisttrack');
            uiwait(labelguiH);
            % read idlist_L
            idlisttrack_L = LG_readIdlistFromOutside;

            if isempty(idlisttrack_L)
                error('checking idlist aborted')
            end

            % store and save
            experiment(iExperiment).idlisttrack_L = idlisttrack_L;
            save(fullfile(experiment(iExperiment).movieDir,...
                experiment(iExperiment).idlisttrack_LName),'idlisttrack_L');

            % update status
            experiment(iExperiment).status = 7;
            % display elapsed time
            disp(sprintf('time elapsed: %i s',...
                round(etime(clock,evaluationStartTime))))

        catch
            err = lasterror;
            experiment(iExperiment) = ...
                handleError(err,...
                experiment(iExperiment),...
                evaluationStartTime);
        end
    end
end % loop experiments
if ishandle(labelguiH)
    close(labelguiH);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
function emi = handleError(err,emi,evaluationStartTime)

% handleError displays meaningful error messages
disp(sprintf('movie ''%s'' failed: %s',...
    emi.projectName,...
    err.message))
for iErr = 1:length(err.stack)
    disp(sprintf( 'in %s at %i',...
        err.stack(iErr).name,err.stack(iErr).line))
end
% set status to failed
emi.status = 99;
% display elapsed time
disp(sprintf('time elapsed: %i s',...
    round(etime(clock,evaluationStartTime))))
