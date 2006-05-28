function varargout = testChromdynExperiment(test,subTest,localDir)
%TESTCHROMDYNEXPERIMENT tests the chromdyn/experiments algorithms
%
% INPUT test: test/list of tests to run (1-by-n)
%       subTest: if a test is split into parts: which parts to run (1-by-n)
%       localDir: local directory for temporary storage of data

% case 1 in-plane straight approach
% case 2 two frames, to stationary tags
% case 3 one tag moving inside a voxel (detector)
% case 4 two tags at different orientations and fractions of resolution
%           limit (detector)
% case 5 one tag moving inside a voxel (tracker on top of detector)
% case 6 two tags at different orientations and fractions of resolution
%           limit (tracker on top of detector)
% case 7 variant 1 of 5, but with transposed movie
% case 8 variant 2 of 5, but with randomized order of points
% case 9 tracker-convergence-test
%
% case 10 in-plane straight crossing (-200%PSF:10%:200%PSF) to allow perfect
%         initial guesses. The movie will be analyzed once with Gaussian,
%         and once with Posson noise and return movie(1:2).detectorFstats,
%         movie(1:2).trackerFstats
%


if nargout
    [varargout{1:nargout}] = deal([]);
end

if nargin < 1 || isempty(test)
    help(testChromdynExperiment)
    return
else
    test = test(:)';
end
if nargin < 2 || isempty(subTest)
    subTest = [];
else
    subTest = subTest(:)';
end
if nargin < 3
    localDir = [];
elseif ~isdir(localDir)
    error('cannot find %s',localDir)
end



% loop through the tests
for currentTest = test

    % switch according to current test
    switch currentTest
        case 1
            %===load or generate raw movie

            % project number 1: in-plane straight approach, equal tag
            % intensity, 13 frames
            projectNumber = 100;

            % load in-plane approach
            positions = loadPositions(projectNumber);

            % set amplitudes
            positions(:,4,:) = 1;
            maxAmplitude = 1;
            minAmplitude = 1;
            % projectNumber = 101;
            % positions(1:2:end,4,:) = 1.1;
            % positions(2:2:end,4,:) = 0.9;
            %             maxAmplitude = 1.1;
            %             minAmplitude = 0.7;

            % load or generate data
            inputDataProperties.movieSize = [64,64,16,13];
            [rawMovie, slist, idlist, dataProperties] = ...
                generateProject(projectNumber, positions,inputDataProperties);

            % snr=10
            snr = 10;
            % poisson noise
            [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude);


            % perturb idlist. All frames are separated, but shifted around
            % by 0.1 pixel
            %             nTimepoints = length(idlist);
            %             fuse = [];
            %fuse = [[4:nTimepoints]',repmat([2,1],nTimepoints-3,1)];
            %             perturbation = 0; % in pixels - 5/20 nm;
            %             idlist_P = perturbIdlist(idlist,fuse,perturbation, dataProperties);

            %             labelguiH = labelgui;
            %             label_loadmovieCB(snrMovie,dataProperties.name,pwd);
            %             labelPanelH = GetUserData(labelguiH,'currentWindow');
            %             SetUserData(labelPanelH,dataProperties,1);
            %             label_loadslistCB(idlist_P,'idlist',slist);
            %             uiwait(labelPanelH)

            % instead of faking something: make slist, idlist
            filteredMovie = filtermovie(snrMovie,dataProperties.FILTERPRM);

            % run spotDetection
            dataProperties.MAXSPOTS = 2;
            [slist] = detectSpots(snrMovie, filteredMovie, dataProperties,2);

            idlist = linker(slist,dataProperties,0);
            idlist(1).info.randomState = randomState;

            labelguiH = LG_loadAllFromOutside(snrMovie,[],[],dataProperties,idlist,'idlist');
            uiwait(labelguiH);
            % theoretically we could read the idlist now
            idlist = LG_readIdlistFromOutside;

            if isempty(idlist)
                return
            end
            %close(labelguiH);

            % run tagTracker
            [idlisttrack,debugData] = tagTracker(snrMovie,idlist,dataProperties,1,struct('fStats',[],'objectiveFunction',1));
            idlisttrack(1).info.randomState = randomState;

            labelguiH = LG_loadAllFromOutside(snrMovie,[],[],dataProperties,idlisttrack,'idlisttrack');

            % assign output
            varargout{1} = idlisttrack;
            varargout{2} = idlist;
            %varargout{3} = idlist_P;

        case 2
            % project number 2: two isolated, stationary tags, two frames
            projectNumber = 101;

            % load two-frame positions
            positions = loadPositions(projectNumber);

            % set amplitudes
            positions(:,4,:) = 1;
            maxAmplitude = 1;
            minAmplitude = 1;

            % load or generate data
            [rawMovie, slist, idlist, dataProperties] = ...
                generateProject(projectNumber, positions);

            % we just do SNR 100
            snr = 100;
            [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude);


            % perturb idlist. All frames are separated, but shifted around
            % by 0.1 pixel
            nTimepoints = length(idlist);
            fuse = [];
            perturbation = 1; % in pixels
            idlist_P = perturbIdlist(idlist,fuse,perturbation, dataProperties);
            idlist_P(1) = idlist(1);

            %             labelguiH = labelgui;
            %             label_loadmovieCB(snrMovie,dataProperties.name,pwd);
            %             labelPanelH = GetUserData(labelguiH,'currentWindow');
            %                         SetUserData(labelPanelH,dataProperties,1);
            %                         label_loadslistCB(idlist_P,'idlist',slist);
            %                         uiwait(labelPanelH)


            % run tagTracker
            [idlisttrack] = tagTracker(snrMovie,idlist_P,dataProperties);

            % compare idlists
            %[A,B,C] = compareIdlists(idlisttrack,idlist,idlist_P);

            % assign output
            varargout{1} = idlisttrack;
            varargout{2} = idlist;
            varargout{3} = idlist_P;

        case 3
            % test for detector (might later be used for tracker, too).
            % One tag in a mini-movie (20x20x10-frame). Tag is placed in
            % one eigth of a voxel on a 0.1-pixel grid to check the
            % influence of position on the accuracy of detection.
            % (the results should be symmetrical along the pixel-diagonal)

            projectNumber = 900;

            % load 216 positions
            positions = loadPositions(projectNumber);

            % leave amplitudes at 1
            maxAmplitude = 1;

            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [20,20,10,216];
            [rawMovie, slist, idlist, dataProperties, dataPath] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            cd(dataPath);
            % run the filter/detector 100 times for every single SNR in the
            % list. Get positions and uncertainties from the slists.
            snrList = [15:-0.5:1.5];
            %snrList = [50,5];
            nSNR = length(snrList);
            nRepeats = 25;
            %nRepeats = 2;
            [deltaCoord,sigmaPos] = deal(zeros(18,36,length(snrList),nRepeats));
            nTotal = nSNR * nRepeats;
            tic
            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats



                    % display status
                    time = toc;
                    timeStr = datestr(datenum(num2str(time),'SS'),13);
                    nCurrent = nDone + k;
                    disp(sprintf(['snr:%2.1f,'...
                        ' #:%i, %i/%i, time elapsed at start %s'],...
                        snr,k,nCurrent,nTotal,timeStr));

                    [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude);

                    % filter the movie
                    filteredMovie = filtermovie(snrMovie,dataProperties.FILTERPRM,0);

                    % run spotDetection
                    cord = spotfind(filteredMovie,dataProperties,0);
                    slist=findoverlap(snrMovie,cord,dataProperties,0);

                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[216,3,1]));
                    for t=1:216
                        if ~isempty(slist(t).sp)
                            deltaTmp(t,[2,1,3],1) = slist(t).sp(1).cord;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,1:3,1);
                            d = diag(slist(t).statistics.Q);
                            qTmp(t,[2,1,3],1) = sqrt(d*slist(t).statistics.chi(1));
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        [reshape(deltaTmp(:,1),[6,36]);...
                        reshape(deltaTmp(:,2),[6,36]);...
                        reshape(deltaTmp(:,3),[6,36])];
                    sigmaPos(:,:,iSNR,k) = ...
                        [reshape(qTmp(:,1),[6,36]);...
                        reshape(qTmp(:,2),[6,36]);...
                        reshape(qTmp(:,3),[6,36])];

                    % save slist
                    slistName = ['slist_S',...
                        num2str(snr),'_i',num2str(k),'.mat'];
                    slist(1).randomState = randomState;
                    save(slistName,'slist');
                end % loop 25x
            end % loop SNR

            deltaOverSigma = deltaCoord./sigmaPos;
            pDelta = normcdf(deltaOverSigma);

            % save everything
            save(['test_3_' nowString],...
                'deltaCoord', 'sigmaPos', 'deltaOverSigma','pDelta')

            % return variables
            varargout = {deltaCoord,sigmaPos,deltaOverSigma,pDelta};

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));



        case 4
            % second test for detector, and later also tracker
            % Two tags in a special-size movie (35x35x24), one fixed, one
            % moving

            projectNumber = 901;

            % load 225 positions
            positions = loadPositions(projectNumber);

            % define other properties. We test tags of different amplitude,
            % at different SNR (respective to fainter tag). Repeat every
            % test 25 times.
            ampList = 1./[1,1.25,1.5,1.75,2,2.4,2.8,3.2,3.6,4];
            snrList = [15:-0.5:1.5];
            %snrList = [50,3];
            nRepeats = 25;
            %nRepeats = 2;
            nAmp = length(subTest);
            nSNR = length(snrList);
            nTotal = nAmp*nSNR*nRepeats;


            % subTest defines which amplitude to take (subTest should be
            % 1-by-n)
            tic % start counting time
            for iAmp = subTest

                % get amplitude to create movie
                minAmp = ampList(iAmp);
                maxAmp = 1;
                positions(:,4,2) = minAmp;

                % make movieNumber from project number: 901_A'1/amp'
                movieNumber = [num2str(projectNumber), '_A',...
                    num2str(1/minAmp)];

                % generate/load movie. Don't forget to supply movieSize
                inputDataProperties.movieSize = [35,35,24,225];
                [rawMovie, slist, idlist, dataProperties, dataPath] = ...
                    generateProject(movieNumber, positions,inputDataProperties);

                % goto local dir
                if ~isempty(localDir)
                    cd(localDir)
                else
                    cd(dataPath);
                end

                % run the filter/detector 25 times for every single SNR in the
                % list. Get positions and uncertainties from the slists.

                %                 oldResultsFileName = [];
                %                 [deltaCoord,sigmaPos] = deal(zeros(45,15,2,1,nRepeats));
                %                 nsp = zeros(225,nSNR,nRepeats);


                for iSNR = 1:nSNR
                    % read snr
                    snr = snrList(iSNR);
                    % count how many iterations we've done so far
                    nDone = (find(iAmp==subTest)-1) * nSNR * nRepeats + ...
                        (iSNR-1) * nRepeats;

                    for k = 1:nRepeats

                        time = toc;
                        timeStr = datestr(datenum(num2str(time),'SS'),13);
                        nCurrent = nDone + k;
                        disp(sprintf(['amp:%1.2f, snr:%2.1f,'...
                            ' #:%i, %i/%i, time elapsed at start %s'],...
                            minAmp,snr,k,nCurrent,nTotal,timeStr));


                        [snrMovie,randomState] = addNoise(rawMovie,snr,minAmp);

                        % filter the movie
                        filteredMovie = ...
                            filtermovie(snrMovie,dataProperties.FILTERPRM,0);

                        % run spotDetection
                        cord = spotfind(filteredMovie,dataProperties,0);
                        slist=findoverlap(snrMovie,cord,dataProperties,0);




                        % save slist. Create name indicative of amplitude
                        % etc, then write random state into slist before
                        % saving
                        slistName = ['slist_A',num2str(1/minAmp),'_S',...
                            num2str(snr),'_i',num2str(k),'.mat'];
                        slist(1).randomState = randomState;
                        save(slistName,'slist');

                    end % loop 25x


                end % loop SNR

                % move everything to dataPath
                if ~isempty(localDir)
                    movefile([localDir,filesep,'*'],dataPath);
                end
            end % loop Amp


            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %=============================================================

        case 5
            % save as case 3: single spot within 1/8th of a voxel, but this
            % time with tracker

            % load project stuff (900)
            % - positions
            % - movie (generateProject)
            tic
            projectNumber = 900;

            % load 216 positions
            positions = loadPositions(projectNumber);

            % leave amplitudes at 1
            maxAmplitude = 1;

            % define SNR, repeats, nSources
            snrList = [15:-0.5:1.5];
            if nargin > 1 && ~isempty(subTest) && all(subTest) > 0
                snrList = snrList(subTest);
            end
            nSNR = length(snrList);
            nRepeats = 25;

            %sourceList = [2,3,4,5,6];
            %nSources = 5;
            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [20,20,10,216];
            [rawMovie, dummy, dummy, dataProperties, dataPath, data2Name] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            % load and sort idlists
            [idlistList, sortNumbers] = ...
                readSynthXlist('idl', dataPath, '_S([\d.]+)_i([\d]+)', [-1,2]);
            idlistList = struct2cell(idlistList);

            goodSnrIdx = ismember(sortNumbers(:,1),snrList);
            idlistList = idlistList(goodSnrIdx);

            % rows: snr, cols: repeats
            idlistList = reshape(idlistList,nRepeats,nSNR)';
            % goto local dir
            if ~isempty(localDir)
                cd(localDir)
                copyfile([dataPath, filesep, data2Name]);
            else
                cd(dataPath);
            end


            micron2pixel = 1./[dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];


            % loop through noise/repeats: run tracker
            [deltaCoord,sigmaPos] = deal(zeros(18,36,length(snrList),nRepeats));
            fitList(1:nSNR,1:nRepeats) = struct('fittingMatrices',[]);

            % deltaStartEnd: nTp, nDim+2, nFits, 2 (=start/end), nSnr, nRepeats
            deltaStartEnd = zeros(216,3+2,5,2,length(snrList),nRepeats);
            % exitflag: nTp, nFits, [exitflag, isSuccess]
            exitFlag = zeros(216,5,2,length(snrList),nRepeats);


            nTotal = nSNR * nRepeats;
            disp(sprintf('time since start: %f',toc));
            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats

                    % display status
                    time = toc;
                    if strmatch('6.5',version)
                        nCurrent = nDone + k;
                        disp(sprintf(['snr:%2.1f,'...
                            ' #:%i, %i/%i, time elapsed at start %f'],...
                            snr,k,nCurrent,nTotal,toc));
                    else
                        timeStr = datestr(datenum(num2str(time),'SS'),13);
                        nCurrent = nDone + k;
                        disp(sprintf(['snr:%2.1f,'...
                            ' #:%i, %i/%i, time elapsed at start %s'],...
                            snr,k,nCurrent,nTotal,timeStr));
                    end

                    % add noise according to idlist(1).info.randomState
                    idlist = idlistList{iSNR,k};
                    randomState = idlist(1).info.randomState;
                    [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude,randomState);

                    % track tags. Use default nSources
                    [idlisttrack,fitList(iSNR,k).fittingMatrices, dse, xfl]...
                        = tagTracker(snrMovie,idlist,dataProperties,-1,1);
                    deltaStartEnd(:,:,:,:,iSNR,k) = dse;
                    exitFlag(:,:,:,iSNR,k) = xfl;


                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[216,3,1]));
                    for t=1:216
                        if ~isempty(idlisttrack(t).linklist)
                            deltaTmp(t,[2,1,3],1) = idlisttrack(t).linklist(1,9:11).*micron2pixel;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,1:3,1);
                            d = diag(idlisttrack(t).info.totalQ_Pix);
                            %qTmp(t,[2,1,3],1) = sqrt(d*idlisttrack(t).linklist(1,12));
                            qTmp(t,[2,1,3],1) = sqrt(d);
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        [reshape(deltaTmp(:,1),[6,36]);...
                        reshape(deltaTmp(:,2),[6,36]);...
                        reshape(deltaTmp(:,3),[6,36])];
                    sigmaPos(:,:,iSNR,k) = ...
                        [reshape(qTmp(:,1),[6,36]);...
                        reshape(qTmp(:,2),[6,36]);...
                        reshape(qTmp(:,3),[6,36])];
                    % save idlisttrack in data2-List as
                    % idlisttrack_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlisttrackName = sprintf('idlisttrack_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlisttrack(1).info.randomState = randomState;
                    saveData2(data2Name,idlisttrackName,idlisttrack);

                end % loop 25x
                save(['t5_intermediate',nowString],...
                    'deltaCoord', 'sigmaPos', 'fitList','deltaStartEnd','exitFlag');
            end % loop SNR


            % save everything
            save(['test_5_' nowString],...
                'deltaCoord', 'sigmaPos', 'fitList','deltaStartEnd','exitFlag');

            % move everything to dataPath
            if ~isempty(localDir)
                movefile([localDir,filesep,'*'],dataPath,'f');
            end

            % return variables
            varargout = {deltaCoord,sigmaPos,fitList,deltaStartEnd,exitFlag};

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %==========================================================

        case 6 % run tracker on slists from 2-spot-test
            % Two tags in a special-size movie (35x35x24), one fixed, one
            % moving

            projectNumber = 901;

            % load 225 positions
            positions = loadPositions(projectNumber);

            % define other properties. We test tags of different amplitude,
            % at different SNR (respective to fainter tag). Repeat every
            % test 25 times.
            ampList = 1./[1,1.25,1.5,1.75,2,2.4,2.8,3.2,3.6,4];
            snrList = [15:-0.5:1.5];
            %snrList = [50,3];
            nRepeats = 25;
            %nRepeats = 2;
            nAmp = length(subTest);
            nSNR = length(snrList);
            nTotal = nAmp*nSNR*nRepeats;


            % subTest defines which amplitude to take (subTest should be
            % 1-by-n)
            disp(sprintf('Start calculation: %s', nowString));
            tic % start counting time
            for iAmp = subTest

                % get amplitude to create movie
                minAmp = ampList(iAmp);
                maxAmp = 1;
                positions(:,4,2) = minAmp;

                % make movieNumber from project number: 901_A'1/amp'
                movieNumber = [num2str(projectNumber), '_A',...
                    num2str(1/minAmp)];

                % generate/load movie. Don't forget to supply movieSize
                inputDataProperties.movieSize = [35,35,24,225];
                [rawMovie, slist, idlist, dataProperties, dataPath, data2Name] = ...
                    generateProject(movieNumber, positions,inputDataProperties);

                % goto dataPath (we have the right amp already!) to load idlists

                [slistList, sortNumbers] = ...
                    readSynthXlist('slist', [dataPath],...
                    '_A([\d.]+)_S([\d.]+)_i([\d]+)',[1,-2,3]);

                % reshape idlist so that every col is another SNR
                slistList = struct2cell(slistList);
                slistList = reshape(slistList,nRepeats,nSNR);

                % goto local dir
                if ~isempty(localDir)
                    cd(localDir)
                    % copy data2-file
                    copyfile([dataPath,filesep,data2Name],localDir);
                else
                    cd(dataPath);
                end

                % run the filter/detector 25 times for every single SNR in the
                % list. Get positions and uncertainties from the slists.
                fitList(1:nSNR,1:nRepeats) = struct('fittingMatrices',[]);
                for iSNR = 1:nSNR
                    % read snr
                    snr = snrList(iSNR);
                    % count how many iterations we've done so far
                    nDone = (find(iAmp==subTest)-1) * nSNR * nRepeats + ...
                        (iSNR-1) * nRepeats;

                    for k = 1:nRepeats

                        time = toc;
                        timeStr = datestr(datenum(num2str(time),'SS'),13);
                        nCurrent = nDone + k;
                        disp(sprintf(['amp:%1.2f, snr:%2.1f,'...
                            ' #:%i, %i/%i, time elapsed at start %s'],...
                            minAmp,snr,k,nCurrent,nTotal,timeStr));
                        % read current idlist
                        currentSlist = slistList{k,iSNR};
                        currentIdlist = slist2idlist(currentSlist,projectNumber,dataProperties);

                        % add the correct noise
                        [snrMovie,randomState] = addNoise(rawMovie,snr,minAmp,currentIdlist(1).info.randomState);

                        % track tags. Use default nSources. No verbosity
                        [idlisttrack,fitList(iSNR,k).fittingMatrices]...
                            = tagTracker(snrMovie,currentIdlist,dataProperties,-1);

                        %LG_loadAllFromOutside(snrMovie,dataPath,[],dataProperties,idlisttrack,'idlisttrack')



                        % save idlisttrack. Create name indicative of amplitude
                        % etc, then write random state into slist before
                        % saving
                        idlisttrackName = sprintf(...
                            'idlisttrack_NS%i_A%1.2f_S%1.2f_i%i.mat',...
                            5,1/minAmp, snr, k);
                        idlisttrack(1).info.randomState = randomState;
                        saveData2(data2Name,idlisttrackName,idlisttrack);

                    end % loop 25x


                end % loop SNR

                % move everything to dataPath
                if ~isempty(localDir)
                    movefile([localDir,filesep,'*'],dataPath);
                end
            end % loop Amp


            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %==========================================================================

        case 7 % like 5, but transposed movie/data
            % save as case 3: single spot within 1/8th of a voxel, but this
            % time with tracker

            % load project stuff (900)
            % - positions
            % - movie (generateProject)
            projectNumber = 900;

            % load 216 positions
            positions = loadPositions(projectNumber);

            % make variant 1
            projectNumber = 900.1;

            % leave amplitudes at 1
            maxAmplitude = 1;

            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [20,20,10,216];
            [rawMovie, dummy, dummy, dataProperties, dataPath, data2Name] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            cd(dataPath);
            micron2pixel = 1./[dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];

            % define SNR, repeats, nSources
            snrList = [15:-0.5:1.5];
            if nargin > 1 && ~isempty(subTest) && subTest > 0
                snrList = snrList(subTest);
            end
            nSNR = length(snrList);
            nRepeats = 25;

            %sourceList = [2,3,4,5,6];
            %nSources = 5;


            % load and sort idlists
            [idlistList, sortNumbers] = ...
                readSynthXlist('idl', pwd, '_S([\d.]+)_i([\d]+)', [-1,2]);
            idlistList = struct2cell(idlistList);
            idlistList = idlistList(end-(nRepeats*nSNR-1):end);
            % rows: snr, cols: repeats
            idlistList = reshape(idlistList,nRepeats,nSNR)';
            %idlistList = idlistList(1);

            % goto local dir
            if ~isempty(localDir)
                cd(localDir)
                copyfile([dataPath, filesep, data2Name]);
            else
                cd(dataPath);
            end


            % loop through noise/repeats: run tracker
            [deltaCoord,sigmaPos] = deal(zeros(18,36,length(snrList),nRepeats));

            nTotal = nSNR * nRepeats;
            tic
            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats

                    % display status
                    time = toc;
                    timeStr = datestr(datenum(num2str(time),'SS'),13);
                    nCurrent = nDone + k;
                    disp(sprintf(['snr:%2.1f,'...
                        ' #:%i, %i/%i, time elapsed at start %s'],...
                        snr,k,nCurrent,nTotal,timeStr));

                    % add noise according to idlist(1).info.randomState
                    idlist = idlistList{iSNR,k};
                    randomState = idlist(1).info.randomState;
                    [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude,randomState);

                    % transpose movie and idlist
                    snrMovie = permute(snrMovie, [2,1,3,4,5]);
                    for t=1:216
                        if ~isempty(idlist(t).linklist)
                            idlist(t).linklist(:,9:10) = ...
                                idlist(t).linklist(:,[10,9]);
                        end
                    end

                    % track tags. Use default nSources
                    [idlisttrack] = tagTracker(snrMovie,idlist,dataProperties,-1);

                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[216,3,1]));
                    % switch also x,y of positions
                    for t=1:216
                        if ~isempty(idlisttrack(t).linklist)
                            deltaTmp(t,[2,1,3],1) = idlisttrack(t).linklist(1,9:11).*micron2pixel;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,[2,1,3],1);
                            d = diag(idlisttrack(t).info.totalQ_Pix);
                            qTmp(t,[2,1,3],1) = sqrt(d*idlisttrack(t).linklist(1,12));
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        [reshape(deltaTmp(:,1),[6,36]);...
                        reshape(deltaTmp(:,2),[6,36]);...
                        reshape(deltaTmp(:,3),[6,36])];
                    sigmaPos(:,:,iSNR,k) = ...
                        [reshape(qTmp(:,1),[6,36]);...
                        reshape(qTmp(:,2),[6,36]);...
                        reshape(qTmp(:,3),[6,36])];
                    % save idlisttrack in data2-List as
                    % idlisttrack_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlisttrackName = sprintf('idlisttrack_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlisttrack(1).info.randomState = randomState;
                    saveData2(data2Name,idlisttrackName,idlisttrack);

                end % loop 25x
            end % loop SNR



            % save everything
            save(['test_5.1_' nowString,'.mat'],...
                'deltaCoord', 'sigmaPos')

            % move everything to dataPath
            if ~isempty(localDir)
                movefile([localDir,filesep,'*'],dataPath);
            end

            % return variables
            varargout = {deltaCoord,sigmaPos};

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %==============================================================

        case 8 % like 5, but with randomized order
            % save as case 3: single spot within 1/8th of a voxel, but this
            % time with tracker

            % load project stuff (900)
            % - positions
            % - movie (generateProject)
            projectNumber = 900;

            % load 216 positions
            positions = loadPositions(projectNumber);

            % set variant 2
            projectNumber = 900.2;

            % leave amplitudes at 1
            maxAmplitude = 1;

            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [20,20,10,216];
            [rawMovie, dummy, dummy, dataProperties, dataPath, data2Name] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            cd(dataPath);
            micron2pixel = 1./[dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];

            % define SNR, repeats, nSources
            snrList = [15:-0.5:1.5];
            if nargin > 1 && ~isempty(subTest) && subTest > 0
                snrList = snrList(subTest);
            end
            nSNR = length(snrList);
            nRepeats = 25;

            %sourceList = [2,3,4,5,6];
            %nSources = 5;


            % load and sort idlists
            [idlistList, sortNumbers] = ...
                readSynthXlist('idl', pwd, '_S([\d.]+)_i([\d]+)', [-1,2]);
            idlistList = struct2cell(idlistList);
            idlistList = idlistList(end-(nRepeats*nSNR-1):end);
            % rows: snr, cols: repeats
            idlistList = reshape(idlistList,nRepeats,nSNR)';
            %idlistList = idlistList(1);

            % goto local dir
            if ~isempty(localDir)
                cd(localDir)
                copyfile([dataPath, filesep, data2Name]);
            else
                cd(dataPath);
            end

            % set random order - keep 1 as the first frame because of
            % idlist(1).xx
            randomOrder = [1,1+randperm(215)];
            % update positions
            positions = positions(randomOrder,:,:);

            % loop through noise/repeats: run tracker
            [deltaCoord,sigmaPos] = deal(zeros(18,36,length(snrList),nRepeats));

            nTotal = nSNR * nRepeats;
            tic
            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats

                    % display status
                    time = toc;
                    timeStr = datestr(datenum(num2str(time),'SS'),13);
                    nCurrent = nDone + k;
                    disp(sprintf(['snr:%2.1f,'...
                        ' #:%i, %i/%i, time elapsed at start %s'],...
                        snr,k,nCurrent,nTotal,timeStr));

                    % add noise according to idlist(1).info.randomState
                    idlist = idlistList{iSNR,k};
                    randomState = idlist(1).info.randomState;
                    [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude,randomState);

                    % make random order
                    snrMovie = snrMovie(:,:,:,:,randomOrder);
                    idlist = idlist(randomOrder);

                    % track tags. Use default nSources
                    [idlisttrack] = tagTracker(snrMovie,idlist,dataProperties,-1);

                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[216,3,1]));
                    for t=1:216
                        if ~isempty(idlisttrack(t).linklist)
                            deltaTmp(t,[2,1,3],1) = idlisttrack(t).linklist(1,9:11).*micron2pixel;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,1:3,1);
                            d = diag(idlisttrack(t).info.totalQ_Pix);
                            qTmp(t,[2,1,3],1) = sqrt(d*idlisttrack(t).linklist(1,12));
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        [reshape(deltaTmp(:,1),[6,36]);...
                        reshape(deltaTmp(:,2),[6,36]);...
                        reshape(deltaTmp(:,3),[6,36])];
                    sigmaPos(:,:,iSNR,k) = ...
                        [reshape(qTmp(:,1),[6,36]);...
                        reshape(qTmp(:,2),[6,36]);...
                        reshape(qTmp(:,3),[6,36])];
                    % save idlisttrack in data2-List as
                    % idlisttrack_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlisttrackName = sprintf('idlisttrack_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlisttrack(1).info.randomState = randomState;
                    idlisttrack(1).info.randomOrder = randomOrder;
                    saveData2(data2Name,idlisttrackName,idlisttrack);

                end % loop 25x
            end % loop SNR


            % save everything
            save(['test_5.2_' nowString],...
                'deltaCoord', 'sigmaPos')

            % move everything to dataPath
            if ~isempty(localDir)
                movefile([localDir,filesep,'*'],dataPath,'f');
            end

            % return variables
            varargout = {deltaCoord,sigmaPos};

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %==========================================================================
        case 9 % test tracker convergence:
            % 6 frames, always same spot position, synthetically very large
            % uncertainty from detector to not influence the result. Run
            % tracker with random bad initial positions and check
            % convergence. Do at 5 SNRs: 50,15,10,5,2.5, repeat 25 times

            % load project stuff (102) - movie is 30x30x16x10
            % - positions
            % - movie (generateProject)
            projectNumber = 102;
            
            error('can''t use test9')

            % load 10 positions
            positions = loadPositions(projectNumber);

            % leave amplitudes at 1
            maxAmplitude = 1;

            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [30,30,15,6];
            [rawMovie, dummy, rawIdlist, dataProperties, dataPath, data2Name] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            cd(dataPath);
            micron2pixel = 1./[dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];

            % define SNR, repeats, nSources
            snrList = [50,15,10,5,2.5];
            if nargin > 1 && ~isempty(subTest) && subTest > 0
                snrList = snrList(subTest);
            end
            nSNR = length(snrList);
            nRepeats = 25;


            % prepare idlist: assign arbitrarily large sigma0 for detector
            for t=1:6
                rawIdlist(t).linklist(:,12) = 1000;
            end

            % goto local dir
            if ~isempty(localDir)
                cd(localDir)
                copyfile([dataPath, filesep, data2Name]);
            else
                cd(dataPath);
            end

            % loop through noise/repeats: run tracker
            [deltaCoord,sigmaPos] = deal(zeros(6,3,length(snrList),nRepeats));
            % deltaStartEnd: nTp, nDim, nFits, S/E, nSnr, nRepeats
            deltaStartEnd = zeros(6,3,6,2,length(snrList),nRepeats);
            % exitflag: nTp, nFits, [exitflag, isSuccess]
            exitFlag = zeros(6,6,2,length(snrList),nRepeats);

            nTotal = nSNR * nRepeats;
            tic

            fitList(1:nSNR,1:nRepeats) = struct('fittingMatrices',[]);

            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats

                    % display status
                    time = toc;
                    timeStr = datestr(datenum(num2str(time),'SS'),13);
                    nCurrent = nDone + k;
                    disp(sprintf(['snr:%2.1f,'...
                        ' #:%i, %i/%i, time elapsed at start %s'],...
                        snr,k,nCurrent,nTotal,timeStr));

                    % add noise
                    [snrMovie,randomState] = addNoise(rawMovie,snr,...
                        maxAmplitude);



                    % track tags. Use default nSources
                    dbOpt.deltaInitPos = [1,1,0.5]; % sigma in pix
                    [idlisttrack, fitList(iSNR,k).fittingMatrices, dse, xfl] =...
                        tagTracker(snrMovie,rawIdlist,...
                        dataProperties,0,dbOpt);
                    deltaStartEnd(:,:,:,:,iSNR,k) = dse;
                    exitFlag(:,:,:,iSNR,k) = xfl;

                    % read results. Loop because there could be no spots
                    % found. Also: adjust dimensions
                    % Only the first spot will be taken. There should be no
                    % others, anyway.
                    [deltaTmp,qTmp] = deal(repmat(NaN,[6,3,1]));
                    for t=1:6
                        if ~isempty(idlisttrack(t).linklist)
                            deltaTmp(t,[2,1,3],1) = idlisttrack(t).linklist(1,9:11).*micron2pixel;
                            deltaTmp(t,:,1) = deltaTmp(t,:,1) - positions(t,1:3,1);
                            d = diag(idlisttrack(t).info.totalQ_Pix);
                            qTmp(t,[2,1,3],1) = sqrt(d*idlisttrack(t).linklist(1,12));
                        end
                    end

                    deltaCoord(:,:,iSNR,k) = ...
                        deltaTmp;
                    sigmaPos(:,:,iSNR,k) = ...
                        qTmp;

                    % save idlisttrack in data2-List as
                    % idlisttrack_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlisttrackName = sprintf('idlisttrack_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlisttrack(1).info.randomState = randomState;
                    saveData2(data2Name,idlisttrackName,idlisttrack);

                end % loop 25x
            end % loop SNR

            % save everything
            save(['test_9_' nowString],...
                'deltaCoord', 'sigmaPos', 'deltaStartEnd','fitList','exitFlag')

            % move everything to dataPath
            if ~isempty(localDir)
                movefile([localDir,filesep,'*'],dataPath,'f');
            end

            % return variables
            varargout = {deltaCoord,sigmaPos,deltaStartEnd};

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

        case 10
            %===load or generate raw movie

            % project number 10: in-plane straight approach, equal tag
            % intensity, 41 frames
            projectNumber = 103;

            % load in-plane approach
            positions = loadPositions(projectNumber);

            % set amplitudes
            positions(:,4,:) = 1;
            maxAmplitude = 1;
            minAmplitude = 1;

            % load or generate data
            inputDataProperties.movieSize = [64,64,16,41];
            [rawMovie, dummy, dummy, dataProperties] = ...
                generateProject(projectNumber, positions,inputDataProperties);

            % ==========

            % check subTest


            % snr=10
            snr = 10;

            out(1:2) = struct('detectorFstats',[],'trackerFstats',[]);

            % loop noise
            %for iNoise = 1:2
            for iNoise = 1
                % 1: gaussian noise
                % 2: poisson noise
                [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude,[],iNoise-1);

                % instead of faking something: make slist, idlist
                filteredMovie = filtermovie(snrMovie,dataProperties.FILTERPRM);

                % run spotDetection - get back fStats
                dataProperties.MAXSPOTS = 2;
                [slist, dummy, dummy, debugData] = ...
                    detectSpots(snrMovie, filteredMovie, dataProperties,2,struct('debug',1));
                % remember fStats
                out(iNoise).detectorFstats = debugData.fStats;

                idlist = linker(slist,dataProperties,1);
                idlist(1).info.randomState = randomState;

                labelguiH = LG_loadAllFromOutside(snrMovie,[],[],dataProperties,idlist,'idlist');
                uiwait(labelguiH);
                % theoretically we could read the idlist now
                idlist = LG_readIdlistFromOutside;

                if isempty(idlist)
                    return
                end
                %close(labelguiH);

                % run tagTracker
                [idlisttrack] = tagTracker(snrMovie,idlist,dataProperties,1);

                %[idlisttrack,debugData] = tagTracker(snrMovie,idlist,dataProperties,1,struct('fStats',[],'objectiveFunction',1));
                idlisttrack(1).info.randomState = randomState;
                out(iNoise).trackerFstats = debugData.fStats;
                out(iNoise).idlisttrack = idlisttrack;
                out(iNoise).idlist = idlist;
                out(iNoise).snrMovie = snrMovie;
                out(iNoise).dataProperties = dataProperties;
                out(iNoise).filteredMovie = filteredMovie;

            end % for iNoise = 1:2

            % assign output
            varargout{1} = out;
            %varargout{3} = idlist_P;
            
        case 11 
            
            % 16 spots, all independent, moving within 1/8th of a voxel 
            % Do at 5 SNRs: 50,15,10,5,2.5

            % load project stuff (902) - movie is 80x80x10x216
            % - positions
            % - movie (generateProject)
            projectNumber = 902;
            
           
            % load positions
            positions = loadPositions(projectNumber);

            % leave amplitudes at 1
            maxAmplitude = 1;

            % generate/load movie. Don't forget to supply movieSize
            inputDataProperties.movieSize = [80,80,10,216];
            [rawMovie, dummy, dummy, dataProperties, dataPath, data2Name] = ...
                generateProject(projectNumber, positions,inputDataProperties);
            cd(dataPath);
            micron2pixel = 1./[dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_XY,...
                dataProperties.PIXELSIZE_Z];

            % define SNR, repeats, nSources
            snrList = [50,15,10,5,2.5];
            if nargin > 1 && ~isempty(subTest) && subTest > 0
                snrList = snrList(subTest);
            end
            nSNR = length(snrList);
            
            % goto local dir
            if ~isempty(localDir)
                cd(localDir)
                copyfile([dataPath, filesep, data2Name]);
            else
                cd(dataPath);
            end

            nRepeats = 1;
            nTotal = nSNR;
            tic
            
            
            for iSNR = 1:nSNR
                % calculate number of already done evaluations
                nDone = (iSNR-1) * nRepeats;
                snr = snrList(iSNR);
                for k = 1:nRepeats

                    % display status
                    time = toc;
                    timeStr = datestr(datenum(num2str(time),'SS'),13);
                    nCurrent = nDone + k;
                    disp(sprintf(['snr:%2.1f,'...
                        ' #:%i, %i/%i, time elapsed at start %s'],...
                        snr,k,nCurrent,nTotal,timeStr));

                    % add noise
                    [snrMovie,randomState] = addNoise(rawMovie,snr,...
                        maxAmplitude);
                    
                    % instead of faking something: make slist, idlist
                filteredMovie = filtermovie(snrMovie,dataProperties.FILTERPRM);

                % run spotDetection - get back fStats
                dataProperties.MAXSPOTS = 16;
                dataProperties.amplitudeCutoff = 6.5;
                [slist] = ...
                    detectSpots(snrMovie, filteredMovie, dataProperties,0);
                
                % set dataProperties for linker
                dataProperties.linker_relAmpWeight = 0.01;
                dataProperties.linker_useCOM = 0;
                
                idlist = linker(slist,dataProperties,0);
                idlist(1).info.randomState = randomState;
                
                % label tags
                idlist = testing_labelTags(idlist,positions,dataProperties);

                % for testing: labelgui
%                 labelguiH = LG_loadAllFromOutside(snrMovie,[],[],dataProperties,idlist,'idlist');
%                 uiwait(labelguiH);
%                 % theoretically we could read the idlist now
%                 idlist = LG_readIdlistFromOutside;
                
                % save idlist in data2-List as
                    % idlist_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlistName = sprintf('idlist_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlist(1).info.randomState = randomState;
                    saveData2(data2Name,idlistName,idlist);

                    % track tags. Use default nSources
                    dbOpt.fStats = [];
                    dbOpt.fitStats = [];
                    dbOpt.trackResults = [];
                    [idlisttrack, debugData] =...
                        tagTracker(snrMovie,idlist,...
                        dataProperties,0,dbOpt);
                    
                   

                    % save idlisttrack in data2-List as
                    % idlisttrack_NS#_S#_i#.mat (# sources, snr, iteration)
                    idlisttrackName = sprintf('idlisttrack_NS%i_S%1.2f_i%i.mat', 5, snr, k);
                    idlisttrack(1).info.randomState = randomState;
                    saveData2(data2Name,idlisttrackName,idlisttrack);
                    
                    
                    debugDataName = sprintf('debugData_P%i_S%1.2f_i%i.mat',projectNumber,snr,k);
                    save(debugDataName,debugData);

                end % loop 1x
            end % loop SNR

           
            % move everything to dataPath
            if ~isempty(localDir)
                movefile([localDir,filesep,'*'],dataPath,'f');
            end

            time = toc;
            timeStr = datestr(datenum(num2str(time),'SS'),13);
            disp(sprintf(['done. total time elapsed %s'],timeStr));

            %==========================================================================

        case 100 %Khuloud's
            
            %assign ...
            projectNumber = 1000; %project number
            movieLength = 5; %length of movie (in seconds)
            numZSlices = 16; %number of z-slices
            timePerSlice = 0.6/numZSlices; %exposure time per z-slice (in seconds)
            pixelSizeXY = 0.010; %pixel size in x and y (in micrometers)
            numPixelXY = 70; %number of pixels in x and y
            
            %determine number of time points in generated MT length trajectory
            numPointsIn1s = ceil(1/timePerSlice);
            numTimePoints = movieLength*numPointsIn1s;
            
            %assign simulation parameter values
            modelParam.growthSpeed = [4 0.3]; %micrometers/minute
            modelParam.shrinkageSpeed = [4 0.3]; %micrometers/minute
            modelParam.growthTime = [0.8 0.3]; %seconds
            modelParam.shrinkageTime = [0.8 0.3]; %seconds
            
            %generate MT length trajectory (in micrometers)
            mtLength0 = mtGammaTdSd(modelParam,0.3,movieLength);
            
            %sample trajectory every timePerSlice seconds
            mtLength = sampleTraj(mtLength0,timePerSlice);
            mtLength = mtLength(1:numTimePoints,2);
            
            %remove from trajectory time points that are not observed
            mtLength = reshape(mtLength,numPointsIn1s,movieLength);
            mtLength = mtLength(1:numZSlices,:);
            mtLength = reshape(mtLength,numZSlices*movieLength,1);
            numTimePoints = length(mtLength);

            %get SPB and CEN positions (in micrometers)
            posSpb = zeros(numTimePoints,3); %spindle pole body position
            posCen = posSpb + [mtLength.*cos(pi/4) ...
                mtLength.*sin(pi/4) zeros(numTimePoints,1)]; %centromere position

            %             %put SPB and CEN in imaging area (coordinates in pixels)
            %             scaleMat = [1/pixelSizeXY*ones(numTimePoints,2) ...
            %                 1/0.2*ones(numTimePoints,1)];
            %             translateMat = repmat(floor([numPixelXY/3 numPixelXY/3 numZSlices/3]),...
            %                 numTimePoints,1);
            %             posSpb = posSpb.*scaleMat + translateMat;
            %             posCen = posCen.*scaleMat + translateMat;

            %             %put positions and amplitudes in one array
            %             positions(:,:,1) = [posSpb ones(numTimePoints,1)];
            %             positions(:,:,2) = [posCen ones(numTimePoints,1)];
            %             maxAmplitude = 1;
            %             minAmplitude = 1;

            %put CEN in imaging area (coordinates in pixels)
            scaleMat = [1/pixelSizeXY*ones(numTimePoints,2) ...
                1/0.2*ones(numTimePoints,1)];
            translateMat = repmat(floor([numPixelXY/2 numPixelXY/2 ...
                numZSlices/2]),numTimePoints,1);
            posCen = (posCen-repmat(posCen(1,:),numTimePoints,1))...
                .*scaleMat + translateMat;

            %put positions and amplitudes in one array
            positions(:,:,1) = [posCen ones(numTimePoints,1)];
            maxAmplitude = 1;
            minAmplitude = 1;

            %generate raw movie
            inputDataProperties.movieSize = [numPixelXY,numPixelXY,...
                numZSlices,movieLength];
            inputDataProperties.PIXELSIZE_XY = pixelSizeXY;
            [rawMovie, slist, idlist, dataProperties] = ...
                generateProject(projectNumber,positions,inputDataProperties,1);

            snr = 100;
            % poisson noise
            [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude);
            
            % instead of faking something: make slist, idlist
            %             filteredMovie = filtermovie(snrMovie,dataProperties.FILTERPRM);
            filteredMovie = snrMovie;

            % run spotDetection
            dataProperties.MAXSPOTS = 2;
            [slist] = detectSpots(snrMovie, filteredMovie, dataProperties,2);

            idlist = linker(slist,dataProperties,0);
            idlist(1).info.randomState = randomState;

            labelguiH = LG_loadAllFromOutside(snrMovie,[],[],dataProperties,idlist,'idlist');
            uiwait(labelguiH);
            % theoretically we could read the idlist now
            idlist = LG_readIdlistFromOutside;

            varargout{1} = idlist;

        otherwise
            disp(sprintf('test %i is not defined',currentTest))
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=======================================================
%   ADD SNR
%=======================================================
function [snrMovie,randomState] = addNoise(rawMovie,snr,maxAmplitude,randomState,poisson);
% adds noise to rawMovie
% snr = maxAmplitude/sigmaNoise
% randomState: state of random generator (for reproducibility)
% poisson: if 1, noise will depend on intensity

if nargin < 4 || isempty(randomState)
    randomState = randn('state');
else
    randn('state',randomState);
end
if nargin < 5 || isempty(poisson)
    poisson = 0;
end

if poisson
    % Have sigma depend on amplitude.
    % var(noise) = var(Gauss) + beta*(I-I0)
    % sigma(noise) = maxAmplitude/SNR
    % In the speckle software, beta is 0.5.
    % Try to choose var(Gauss) so that SNRmovie is equal to SNR where the
    % amplitude is 0.5*maxAmplitude.
    % Thus:
    sigmaNoise = ...
        sqrt( ((maxAmplitude/snr)^2 - 0.25*maxAmplitude) + 0.5 * rawMovie);
else
    % use uniform gaussian noise
    sigmaNoise = maxAmplitude/snr;
end

% add noise, ensure that intensity is never below 0
snrMovie = randn(size(rawMovie)) .* sigmaNoise;
snrMovie = rawMovie + snrMovie;
snrMovie = snrMovie - min(snrMovie(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=======================================================
%   PERTURB IDLIST
%=======================================================
function [idlist_P] = perturbIdlist(idlist, fuse, perturbation, dataProperties)
% perturbs idlist.
% fuse is a list of fusions in the form [t, tag#, tag# of new position]
% perturbation (in pixels): standard deviation of the random perturbation
% of the positions. Scalar or 1-by-3 matrix

% calculate perturbation in microns and make correct matrix size
nTags = length(idlist(1).stats.labelcolor);
perturbation = perturbation .* ...
    [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];
perturbation = repmat(perturbation, nTags, 1);

lastgood = 0;
for t = 1:length(idlist)
    if ~isempty(idlist(t).linklist)
        % perturb
        idlist(t).linklist(:,9:11) = idlist(t).linklist(:,9:11) + ...
            randn(nTags,3) .* perturbation;

        % check whether we have to fuse a spot
        if ~isempty(fuse)
            fuseIdx = find(fuse(:,1)==t);
            for i = 1:length(fuseIdx)
                idlist(t).linklist(fuse(fuseIdx(i),2),[2,9:11]) = ...
                    idlist(t).linklist(fuse(fuseIdx(i),3),[2,9:11]);
            end

            % adjust linkup/linkdown
            if lastgood > 0
                idlist(lastgood).linklist(:,7) = idlist(t).linklist(:,2);
                idlist(t).linklist(:,6) = idlist(lastgood).linklist(:,2);
            end

            % remember last good frame
            lastgood = t;
        end % if ~isempty(fuse)
    end % if ~isempty(idlist(t).linklist)
end % for t = 1:length(idlist)

idlist_P = idlist;



