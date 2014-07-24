function varargout=testFuzzyGroupArma(test,varargin)
%TESTFUZZYGROUPARMA is a test routine for fuzzyGroupArma
%
% SYNOPSIS: out=testFuzzyGroupArma(test)
%
% INPUT test: index of the test
%       varargin: for test 2, 3
%                   1) trueArma parameters of model
%                   2) addNoise
%
% OUTPUT out: test data
%
% REMARKS
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: jdorn
% DATE: 25-Jan-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% this is a bit of a hack at the moment, until I get the structure of it
% all sorted out


switch test

    case 0
        % tests possiblilistic clustering

        % set up scene
        n1 = 500; %figure: 500/500
        n2 = 200; %figure: 500/200
        nu = 500; %figure: 500
        nn = 1500; %figure: 1500

        % set default cluster parameters
        eta = -1;
        m = 1.5;
        initialGuess = struct('mean',{[4,1],[2,5]});
        c1 = [2,2];
        c2 = [4,3.5];
        verbose = true;

        if nargin > 1 && ~isempty(varargin{1})
            testOpt = varargin{1};
            if isfield(testOpt,'c1')
                c1 = testOpt.c1;
                c2 = testOpt.c2;
            end
            % assign n2
            if isfield(testOpt,'n2')
                n2 = testOpt.n2;
            end
            % set state of random generators
            if isfield(testOpt,'randState')
                rand('state',testOpt.randState);
                randn('state',testOpt.randnState);
            end
            if isfield(testOpt,'initialGuess')
                initialGuess = testOpt.initialGuess;
            end
            if isfield(testOpt,'verbose')
                verbose = testOpt.verbose;
            end
            if isfield(testOpt,'nn')
                nn = testOpt.nn;
            end
            if isfield(testOpt,'nu');
                nu = testOpt.nu;
            end
        end
        if nargin > 2 && ~isempty(varargin{2})
            % assign pcOptions
            m = varargin{2};

        else
            % default m has already been assigned
        end


        % cluster1 is at 2,2
        cluster1 = randn(n1,2) + repmat(c1,n1,1);

        % cluster2 is at 3,2.5
        cluster2 = randn(n2,2) + repmat(c2,n2,1);

        % extend the figure from -1 to 7
        outsideIdx1 = find(cluster1>7 | cluster1<-1);
        cluster1(outsideIdx1) = rand(size(outsideIdx1)) * 10 - 2;
        outsideIdx2 = find(cluster2>7 | cluster2<-1);
        cluster2(outsideIdx2) = rand(size(outsideIdx2)) * 10 - 2;

        % add no uniform noise
        uniformNoise = rand(nu,2) * 8 - 1;


        % add nonuniform noise. Distribute so that more is outside. Make
        % Gaussians from the border
        nn2 = floor(nn/2);
        nonuniformNoise = randn(nn2,2);
        aboveZeroIdxL = nonuniformNoise>0;
        belowZeroIdxL = nonuniformNoise<0;
        nonuniformNoise(aboveZeroIdxL) = nonuniformNoise(aboveZeroIdxL) - 2;
        nonuniformNoise(belowZeroIdxL) = nonuniformNoise(belowZeroIdxL) + 8;

        nonuniformNoise = [nonuniformNoise(:,1),rand(nn2,1)*10-2;...
            rand(nn-nn2,1)*10-2,nonuniformNoise(:,2)];

        % collect coordinates
        coords = [cluster1;cluster2;uniformNoise;nonuniformNoise];



        centerFunction = 'fcmCenterFunction';
        distanceFunction = 'fcmDistanceFunction';


        % run fuzzy clustering
        [centers, dataWeight, distances, eta, typicality, membership,diagnostics] = ...
            possibilisticClustering(coords,initialGuess,centerFunction,...
            distanceFunction,eta,m);

        % plot results
        classification = blkdiag(ones(n1,1),ones(n2,1),...
            ones(nu+nn,1));
        if verbose
            possibilisticClusteringPlot(coords,typicality,centers,[],[],classification,cat(1,initialGuess.mean),[2,2;4,3.5]);
        end
        if nargout > 0
            varargout{1} = centers;
            varargout{2} = coords;
            varargout{3} = classification;
        end
        %% Test with gaussians
    case -1
        % test possibilisticClustering

        currentTest = varargin{1};
        % current test:
        % 1: equally-sized clusters,
        % 2: unequally sized clusters
        %       start at 16 different positions
        %       run 100x for each algorithm
        % 3: as 1, but random starting positions [0-6]
        % 4: as 2, but random starting positions [0-6]
        %       start at 100 different positions
        %       use 100 different data sets

        % match centers via lap
        % store end-positions, deltaPos

        % check for no noise (negative test number
        noNoise = sign(currentTest)/2+0.5;
        
        % check for smaller start area
        smallStart = currentTest - round(currentTest) ~= 0;
        if smallStart
            currentTest = round(currentTest);
        end

        % initialize
        switch abs(currentTest)
            case 1

                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 500;

                if noNoise
                    testOpt.nn = 0;
                    testOpt.nu = 0;
                end

                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                startGrid = 1:2:5;
                nStartGrid = length(startGrid);
                nStart = nStartGrid^2*(nStartGrid^2+1)/2;

                % start positions: x= 0:2:6, y = 0:2:6
                % make 136 unique coordinate pairs [xRed,yRed,xGreen,yGreen]
                [xx1,yy1,xx2,yy2]=ndgrid(startGrid,startGrid,...
                    startGrid,startGrid);
                startCoords = [xx1(:),yy1(:),xx2(:),yy2(:)];
                % make indices for coords
                [startIdx1,startIdx2] = ndgrid(1:nStartGrid^2,1:nStartGrid^2);
                % cut redundant indices, coords
                triIdx = find(tril(ones(nStartGrid^2)));
                startCoords = startCoords(triIdx,:);
                startIdx1 = startIdx1(triIdx);
                startIdx2 = startIdx2(triIdx);

            case 2
                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 200;

                if noNoise
                    testOpt.nn = 0;
                    testOpt.nu = 0;
                end

                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                startGrid = 1:2:5;
                nStartGrid = length(startGrid);
                nStart = nStartGrid^2*(nStartGrid^2+1)/2;

                % start positions: x= 0:2:6, y = 0:2:6
                % make 136 unique coordinate pairs [xRed,yRed,xGreen,yGreen]
                [xx1,yy1,xx2,yy2]=ndgrid(startGrid,startGrid,...
                    startGrid,startGrid);
                startCoords = [xx1(:),yy1(:),xx2(:),yy2(:)];
                % make indices for coords
                [startIdx1,startIdx2] = ndgrid(1:nStartGrid^2,1:nStartGrid^2);
                % cut redundant indices, coords
                triIdx = find(tril(ones(nStartGrid^2)));
                startCoords = startCoords(triIdx,:);
                startIdx1 = startIdx1(triIdx);
                startIdx2 = startIdx2(triIdx);

            case 3
                % random start

                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 500;

                if noNoise
                    testOpt.nn = 0;
                    testOpt.nu = 0;
                end

                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                nStart = 100;

                % startGrid : nStartx4 array of xs1,ys1, xs2, ys2
                if smallStart
                    startCoords = rand(nStart,4)*4+1; % from 1 to 5
                else
                    startCoords = rand(nStart,4)*6; % from 0 to 6
                end

            case 4
                % random start, unequal clusters

                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 200;

                if noNoise
                    testOpt.nn = 0;
                    testOpt.nu = 0;
                end

                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                nStart = 100;

                % startCoords : nStartx4 array of xs1,ys1, xs2, ys2
               if smallStart
                    startCoords = rand(nStart,4)*4+1; % from 1 to 5
                else
                    startCoords = rand(nStart,4)*6; % from 0 to 6
                end

            case 5
                % random start
                % as 3, but with less noise. No uniform, only 500
                % nonuniform noise spots

                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 500;


                testOpt.nn = 250;
                testOpt.nu = 0;


                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                nStart = 100;

                % startGrid : nStartx4 array of xs1,ys1, xs2, ys2
                if smallStart
                    startCoords = rand(nStart,4)*4+1; % from 1 to 5
                else
                    startCoords = rand(nStart,4)*6; % from 0 to 6
                end

            case 6
                % random start, unequal clusters
                % no uniform, only 500 nonuniform noise spots

                % true centers
                trueMeans = [2,2;...
                    4,3.5];

                testOpt.c1 = trueMeans(1,:);
                testOpt.c2 = trueMeans(2,:);
                testOpt.verbose = false;
                pcOpt.m = 1.5;
                testOpt.n2 = 200;


                testOpt.nn = 250;
                testOpt.nu = 0;


                % nIterations
                nRep = 100;
                algList = [1,2,3];
                nAlg = length(algList);
                nStart = 100;

                % startCoords : nStartx4 array of xs1,ys1, xs2, ys2
                if smallStart
                    startCoords = rand(nStart,4)*4+1; % from 1 to 5
                else
                    startCoords = rand(nStart,4)*6; % from 0 to 6
                end

            otherwise
                % no other test yet
        end


        % set up store array
        switch currentTest
            case {-2,-1,1,2}
                store = zeros(nRep,11,nStart,nAlg);
            case {-4,-3,3,4}
                store = zeros(nRep,13,nStart,nAlg);
        end

        disp(sprintf('\n%s',repmat(' ',1,51)));
        cumTime = 0;
        switch abs(currentTest)
            case {3,4,5,6}
                for rep = 1:nRep
                    % create scene. Draw some random numbers, then remember
                    % state
                    rand(1000);
                    randn(1000);
                    testOpt.randState = rand('state');
                    testOpt.randnState = randn('state');
                    for alg = algList
                        % set algorithm
                        pcOpt.algorithm = alg;
                        for start = 1:nStart
                            % set start values
                            testOpt.initialGuess = struct('mean',...
                                {startCoords(start,1:2),...
                                startCoords(start,3:4)});
                            tic
                            % cluster
                            [centers,coords,classification]=testFuzzyGroupArma(0,testOpt,pcOpt);
                            t=toc;
                            % assign centers
                            centerMeans = cat(1,centers.mean);

                            dm = distMat2(trueMeans,centerMeans);
                            sortIdx = lap(dm);
                            centerMeans = centerMeans(sortIdx,:);
                            sc = [startCoords(start,1:2);startCoords(start,3:4)];
                            sc = sc(sortIdx,:);

                            % store
                            store(rep,:,start,(alg==algList)) = ...
                                [centerMeans(1,:),...
                                centerMeans(1,:) - trueMeans(1,:),...
                                centerMeans(2,:),...
                                centerMeans(2,:) - trueMeans(2,:),...
                                sc(1,:),sc(2,:),...
                                t];
                            cumTime = cumTime + t;
                            disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\brep %4.0f/%4.0f alg %4.0f/%4.0f sta %4.0f/%4.0f t %5.0f',...
                                rep,nRep,find(alg==algList),nAlg,start,nStart,cumTime))
                        end % loop start
                    end % loop alg
                end % loop rep

            case {1,2}
                for rep = 1:nRep
                    % create scene. Draw some random numbers, then remember
                    % state
                    rand(1000);
                    randn(1000);
                    testOpt.randState = rand('state');
                    testOpt.randnState = randn('state');
                    for alg = algList
                        % set algorithm
                        pcOpt.algorithm = alg;
                        for start = 1:nStart
                            % set start values
                            testOpt.initialGuess = struct('mean',...
                                {startCoords(start,1:2),...
                                startCoords(start,3:4)});
                            tic
                            % cluster
                            [centers,coords,classification]=testFuzzyGroupArma(0,testOpt,pcOpt);
                            t=toc;
                            % assign centers
                            centerMeans = cat(1,centers.mean);

                            dm = distMat2(trueMeans,centerMeans);
                            sortIdx = lap(dm);
                            centerMeans = centerMeans(sortIdx,:);

                            % store
                            store(rep,:,start,(alg==algList)) = ...
                                [centerMeans(1,:),...
                                centerMeans(1,:) - trueMeans(1,:),...
                                centerMeans(2,:),...
                                centerMeans(2,:) - trueMeans(2,:),...
                                startIdx1(start),startIdx2(start),...
                                t];
                            cumTime = cumTime + t;
                            disp(sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\brep %4.0f/%4.0f alg %4.0f/%4.0f sta %4.0f/%4.0f t %5.0f',...
                                rep,nRep,find(alg==algList),nAlg,start,nStart,cumTime))
                        end % loop start
                    end % loop alg
                end % loop rep

            otherwise
                error('no such test')
        end

        % plot all. In order to not overload plot with too many handles,
        % make one green and one red plot
        switch abs(currentTest)
            case {1,2}
                vectorScale = 0.2;
                for alg = 1:nAlg
                    figure('Name',sprintf('Algorithm %i',algList(alg)));
                    center1 = NaN(3*nRep*nStart,2);
                    % center: [startXY;startXY+delta;NaN,NaN]*n
                    center1(1:3:end,:) = [reshape(store(:,9,:,alg),[],1),...
                        reshape(store(:,10,:,alg),[],1)];
                    center1(2:3:end,:) = center1(1:3:end,:) + ...
                        [reshape(store(:,3,:,alg),[],1),...
                        reshape(store(:,4,:,alg),[],1)] * vectorScale;

                    center2 = NaN(3*nRep*nStart,2);
                    % center: [startXY;startXY+delta;NaN,NaN]*n
                    center2(1:3:end,:) = center1(1:3:end,:);
                    center2(2:3:end,:) = center2(1:3:end,:) + ...
                        [reshape(store(:,7,:,alg),[],1),...
                        reshape(store(:,8,:,alg),[],1)] * vectorScale;

                    plot(center1(:,1),center1(:,2),'+-r',...
                        center2(:,1),center2(:,2),'.-g',startIdx1,startIdx2,'ok')
                    % plot key to starting positions
                    [xxt,yyt]=ndgrid(2:nStartGrid+1,nStartGrid^2-nStartGrid+1:nStartGrid^2);
                    for i=1:nStartGrid^2,text(xxt(i),yyt(i),sprintf('%2.0f',i));end
                    for i=0:nStartGrid-1,text(1,nStartGrid^2-nStartGrid+i+1,sprintf('%2.0f',startGrid(i+1)),'Color','r');end
                    for i=0:nStartGrid-1,text(2+i,nStartGrid^2-nStartGrid,sprintf('%2.0f',startGrid(i+1)),'Color','r');end
                    xlim([0,nStartGrid^2+1])
                    ylim([0,nStartGrid^2+1])
                end
                save(sprintf('testPCM_store_t%i_%s',currentTest,nowString),'store','nStartGrid','nRep','nStart','nAlg','startIdx1','startIdx2','startGrid')

            case {3,4,5,6}
                % make list as output

                %       C1      C2
                % A1 - deltaMean, sigmaMax, ... % conv, %coinc
                % A2
                % A3

                out = zeros(3,6);
                nAll = nRep * nStart;
                clusterIdx = [3,4;7,8];

                %                 fhc = figure('Name','Convergence');
                %                 ahc(1) = subplot(2,2,1);hold on
                %                 ahc(2) = subplot(2,2,2);hold on
                %                 ahc(3) = subplot(2,2,3);hold on

                fhr = figure('Name','Results');
                for alg = 1:nAlg
                    ahr(alg) = subplot(2,2,alg);hold on
                end

                colors = 'rg';
                %
                % loop algorithm
                for alg = 1:nAlg
                    % find coincident clusters
                    storeAlg = store(:,:,:,alg);

                    % reshape to have list of coords. Make so that we have
                    % different starting points before different data sets
                    storeAlg = permute(storeAlg,[2,3,1]);
                    storeAlg = reshape(storeAlg,13,[])';




                    % if one cluster diverges, the other
                    % might be incorrectly labeled, and it may be at the
                    % wrong place. I define divergence as being outside
                    % [0,6], and do not count the position of the other
                    % cluster in that case
                    divIdxL = any(storeAlg(:,[1,2,5,6])>6|storeAlg(:,[1,2,5,6])<0,2);

                    out(alg,5) = sum(divIdxL)/nAll;

                    % find coincident clusters (but not divergent)
                    coincIdxL = all(isApproxEqual(...
                        storeAlg(:,[1,2]),storeAlg(:,[5,6]),0.5,'absolute'),2)...
                        & ~divIdxL;

                    % store percentage
                    out(alg,6) = sum(coincIdxL)/nAll;


                    % loop clusters
                    for c = 1:2

                        clusterPos = storeAlg(:,clusterIdx(c,:));
                        %startPos = storeAlg(:,clusterIdx(c,:)+6);



                        % goodIdx is all nice ones
                        goodIdxL = ~divIdxL & ~coincIdxL;







                        %

                        % delta
                        meanPos = mean(clusterPos(goodIdxL,:),1);
                        out(alg,1+(c-1)*2) = norm(meanPos);
                        if any(goodIdxL)
                            % largest eingenValue
                            covariance = cov(clusterPos(goodIdxL,:));
                            varMax = max(eig(covariance));
                            out(alg,2+(c-1)*2) = sqrt(varMax);
                        end

                        % plot results
                        if c==1
                            % plot coincident
                            plot(ahr(alg),...
                                clusterPos(coincIdxL,1)+trueMeans(c,1),clusterPos(coincIdxL,2)+trueMeans(c,2),'.b');
                        end
                        plot(ahr(alg),...
                            clusterPos(goodIdxL,1)+trueMeans(c,1),clusterPos(goodIdxL,2)+trueMeans(c,2),['+',colors(c)],...
                            clusterPos(divIdxL,1)+trueMeans(c,1),clusterPos(divIdxL,2)+trueMeans(c,2),['.',colors(c)]);


                        % plot true positions
                        plot(ahr(alg),trueMeans(c,1),trueMeans(c,2),'ok')

                        % plot error ellipses twice: around meanPos, and
                        % around [0,0]
                        if any(goodIdxL)
                            plot(ahr(alg),meanPos(:,1),meanPos(:,2),['+',colors(c)])

                            eh=errorEllipse(ahr(alg),covariance,meanPos);
                            set(eh,'Color',colors(c))

                            plot(ahr(alg),0,0,'ok')
                            % change meanPos to true positions
                            meanPos = meanPos + trueMeans(c,:);
                            plot(ahr(alg),meanPos(:,1),meanPos(:,2),'+k')
                            eh=errorEllipse(ahr(alg),covariance,meanPos);
                            set(eh,'Color','k')
                        end


                    end
                end

                disp(out)

        end

        fh = possibilisticClusteringPlot(coords,[],[],[],[],classification,[],[2,2;4,3.5]);
        set(fh,'Name','Sample Data')



        save(sprintf('testPCM_store_t%i_%s',currentTest,nowString),'store','out')

        %%
    case 1

        % map arma models in the vicinity of experimental ones and find
        % the p-values between them
        ar = 0.79:0.025:0.99;
        nAr = length(ar);
        ma = -0.1:-0.05:-0.5;
        nMa = length(ma);
        nTimepoints = 2000;

        arma(nAr,nMa) = struct('traj',[],'fitResults',[]);

        for a = 1:nAr
            for m = 1:nMa

                % trajectory - wnv = 0.05^2
                arma(a,m).traj = simSetarma([],[],ar(a),ma(m),0.05,nTimepoints);

                arma(a,m).traj = perturbTrajectory(arma(a,m).traj,0.01,0.05);

                % fit - first get initial guess, then fit
                initialGuess.arParamP0 = inverseLevDurbExpoAR(ar(a));
                initialGuess.maParamP0 = inverseLevDurbExpoMA(ma(m));

                arma(a,m).fitResults = armaxFitKalman(arma(a,m).traj,[],initialGuess,'tl');

            end
        end


        % loop to get p-values
        %         totalN = nAr*nMa;
        %         pVals = zeros(totalN);
        %         for idx = 1:totalN
        %             for jdx = 1:totalN
        %                 pVals(idx,jdx)=armaxModelComp(arma(idx).fitResults,...
        %                     arma(jdx).fitResults);
        %             end
        %         end

        % pValues: to find the distribution for point (ar1,ma1):
        % pVals(:,:,1,1)
        pVals = zeros(nAr,nMa,nAr,nMa);
        for a1 = 1:nAr
            for m1 = 1:nMa
                for a2 = 1:nAr
                    for m2 = 1:nMa
                        pVals(a1,m1,a2,m2)=armaxModelComp(arma(a1,m1).fitResults,...
                            arma(a2,m2).fitResults);
                    end
                end
            end
        end

        % plot individual models
        figure
        hold on
        for a = 1:nAr
            for m = 1:nMa
                plot(ar(a),ma(m),'+k');
                plot(arma(a,m).fitResults.arParamK(1),arma(a,m).fitResults.maParamK(1),'or');
                line([ar(a),arma(a,m).fitResults.arParamK(1)],[ma(m),arma(a,m).fitResults.maParamK(1)])
                errorEllipse(arma(a,m).fitResults.varCovMatF,...
                    [arma(a,m).fitResults.arParamK(1),arma(a,m).fitResults.maParamK(1)],...
                    'conf',-1);

            end
        end

        varargout = {arma,pVals};

    case 2
        % test clustering

        def_trueArma = {[0.9,-0.3] [0.8,-0.4]};
        def_addNoise = true;
        def_nTimepoints = 2000;
        def_trajLength = 100;
        if nargin == 1 || isempty(varargin{1})
            trueArma = def_trueArma;
        else
            trueArma = varargin{1};
        end
        if nargin < 3 || isempty(varargin{2})
            addNoise = def_addNoise;
        else
            addNoise = varargin{2};
        end
        if nargin < 4 || isempty(varargin{3})
            nTimepoints = def_nTimepoints;
        else
            nTimepoints = varargin{3};
        end
        if nargin< 5 || isempty(varargin{4})
            trajLength=def_trajLength;
        else
            trajLength=varargin{4};
        end

        % generate model - it's fast
        %trueArma = {[0.9,-0.3] [0.8,-0.4]};
        wnSigma = 0.05;
        %         nTimepoints = 10000;
        %         trajLength = 500;
        nTraj = ceil(nTimepoints/trajLength);

        % remember random state for reproducibility
        randomState = rand('state');

        % generate arma model and fit
        arma(1:2) = struct('traj',[],'fitResults',[],'individual',[]);
        for i=1:2
            arma(i).traj = simSetarma([],[],...
                trueArma{i}(1),trueArma{i}(2),wnSigma,nTimepoints);
            % perturb model
            if addNoise
                arma(i).traj = perturbTrajectory(arma(i).traj,0.01,wnSigma);
            end
            % trajOut are the inputs for armaxFitKalman
            tmp(1:nTraj) = struct('observations',[]);
            arma(i).individual = tmp;
            for j = 1:nTraj
                arma(i).individual(j).observations = ...
                    arma(i).traj((j-1)*trajLength+1:min(nTimepoints,j*trajLength),:);
            end

            % get the initial guess
            initialGuess.arParamP0 = inverseLevDurbExpoAR(trueArma{i}(1));
            initialGuess.maParamP0 = inverseLevDurbExpoMA(trueArma{i}(2));

            % fit - overall fit
            arma(i).fitResults = ...
                armaxFitKalman(arma(i).individual,[],initialGuess,'tl');
            % individual fits
            for j=1:nTraj
                arma(i).individual(j).fitResults = ...
                    armaxFitKalman(arma(i).individual(j),[],initialGuess,'tl');
                % complete structure
                arma(i).individual(j).fitResults.numObserve = ...
                    size(arma(i).individual(j).observations,1);
                % name has to contain "_movie"
                arma(i).individual(j).fitResults.name = sprintf('arma_%i_traj_%i_movie',i,j);
                arma(i).individual(j).fitResults.orderLen = [length(trueArma{i}(1)),length(trueArma{i}(2))];
                arma(i).individual(j).fitResults.type = 'Len';
            end
        end

        % plot all the results
%         plotFigure = figure;
%         style = 'xo';
%         hold on
%         for i=1:2,
%             % plot center
%             plot(arma(i).fitResults.arParamK(1),...
%                 arma(i).fitResults.maParamK(1),[style(i) 'r']);
%             % plot error ellipse associated with center
%             h=errorEllipse(arma(i).fitResults.varCovMatF,...
%                 [arma(i).fitResults.arParamK(1),arma(i).fitResults.maParamK(1)],...
%                 'conf',-1);
%             set(h,'Color','r','Tag','')
% 
% 
%             for j=1:nTraj
%                 plot(arma(i).individual(j).fitResults.arParamK(1),...
%                     arma(i).individual(j).fitResults.maParamK(1),[style(i) 'k']);
%                 errorEllipse(arma(i).individual(j).fitResults.varCovMatF,...
%                     [arma(i).individual(j).fitResults.arParamK(1),...
%                     arma(i).individual(j).fitResults.maParamK(1)],...
%                     'conf',-1);
% 
%             end
%         end


        for i=2:-1:1
            for j=nTraj:-1:1
                clusterInput((i-1)*nTraj+j) = arma(i).individual(j).fitResults;
            end
        end
        for i=2:-1:1
            for j=nTraj:-1:1
                clusterInput((i-1)*nTraj+j).lenSeries.observations = arma(i).individual(j).observations;
            end
        end
        armaxModelComp(arma(1).fitResults,arma(2).fitResults)
        fuzzyGroupArmaTmp(clusterInput');
        
        
%         options.armaThreshold = 0.05;
%         options.verbose = 1;
%         
%         [attached, detached] = fuzzyGroupArma(clusterInput,options);
% 
%         figure(plotFigure)
%         for i=1:length(attached)
%             plot(attached(i).arParamK(1),...
%                 attached(i).maParamK(1),'xg');
%             % plot error ellipse associated with center
%             h=errorEllipse(attached(i).varCovMatF,...
%                 [attached(i).arParamK(1),attached(i).maParamK(1)],...
%                 'conf',-1);
%             set(h,'Color','g','Tag','')
%         end
% 
%         arma(1).randomState = randomState;
% 
%         varargout{1} = arma;
%         varargout{2} = attached;
%         varargout{3} = clusterInput;

    case 3
        % find distribution of ARMA coefficients from 100tp movies

        def_trueArma = {[0.9,-0.3] [0.8,-0.4]};
        def_addNoise = true;
        if nargin == 1 || isempty(varargin{1})
            trueArma = def_trueArma;
        else
            trueArma = varargin{1};
        end
        if nargin < 3 || isempty(varargin{2})
            addNoise = def_addNoise;
        else
            addNoise = varargin{2};
        end

        % generate model - it's fast
        %trueArma = {[0.9,-0.2] [0.3,-0.7],[0.9,0.6],[0.6,0.9]};
        nArma = length(trueArma);
        wnSigma = 0.05;
        nTimepoints = 500000;
        trajLength = 100;
        nTraj = ceil(nTimepoints/trajLength);

        % remember random state for reproducibility
        randomState = rand('state');

        % generate arma model and fit
        arma(1:nArma) = struct('traj',[],'fitResults',[],'individual',[]);
        for i=1:nArma
            arma(i).traj = simSetarma([],[],...
                trueArma{i}(1),trueArma{i}(2),wnSigma,nTimepoints);
            % perturb model
            if addNoise
                arma(i).traj = perturbTrajectory(arma(i).traj,0.01,wnSigma);
            end
            % trajOut are the inputs for armaxFitKalman
            tmp(1:nTraj) = struct('observations',[]);
            arma(i).individual = tmp;
            for j = 1:nTraj
                arma(i).individual(j).observations = ...
                    arma(i).traj((j-1)*trajLength+1:min(nTimepoints,j*trajLength),:);
            end

            % get the initial guess
            initialGuess.arParamP0 = inverseLevDurbExpoAR(trueArma{i}(1));
            initialGuess.maParamP0 = inverseLevDurbExpoMA(trueArma{i}(2));

            % fit - overall fit
            arma(i).fitResults = ...
                armaxFitKalman(arma(i).individual,[],initialGuess,'tl');
            % individual fits
            for j=1:nTraj
                arma(i).individual(j).fitResults = ...
                    armaxFitKalman(arma(i).individual(j),[],initialGuess,'tl');
                % complete structure
                arma(i).individual(j).fitResults.numObserve = ...
                    size(arma(i).individual(j).observations,1);
                % name has to contain "_movie"
                arma(i).individual(j).fitResults.name = sprintf('arma_%i_traj_%i_movie',i,j);
                arma(i).individual(j).fitResults.orderLen = [length(trueArma{i}(1)),length(trueArma{i}(2))];
            end
        end

        % plot all the results
        plotFigure = figure;
        style = {'.r','.b','.g','.m'};
        hold on
        for i=1:nArma
            % plot center
            plot(arma(i).fitResults.arParamK(1),...
                arma(i).fitResults.maParamK(1),['o',style{i}(2)]);
            % plot error ellipse associated with center
            %             h=errorEllipse(arma(i).fitResults.varCovMatF,...
            %                 [arma(i).fitResults.arParamK(1),arma(i).fitResults.maParamK(1)],...
            %                 'conf',-1);
            %             set(h,'Color','r','Tag','')


            for j=1:nTraj
                plot(arma(i).individual(j).fitResults.arParamK(1),...
                    arma(i).individual(j).fitResults.maParamK(1),style{i});
                %                 errorEllipse(arma(i).individual(j).fitResults.varCovMatF,...
                %                     [arma(i).individual(j).fitResults.arParamK(1),...
                %                     arma(i).individual(j).fitResults.maParamK(1)],...
                %                     'conf',-1);

            end
        end

        % get distance - from (1) trajectory of (3) center to (2) all
        % centers
        distances = zeros(nTraj,nArma,nArma);
        for myCenter=1:nArma
            for allCenters =1:nArma
                for traj = 1:nTraj
                    distances(traj,allCenters,myCenter) = ...
                        armaxModelComp(arma(myCenter).individual(traj).fitResults,...
                        arma(allCenters).fitResults);
                end
            end
        end

        figure,plot(reshape(distances,[nTraj*nArma,nArma]));

        varargout{1} = arma;
        varargout{2} = distances;


    case 4
        % one single cluster. Find the correct center, test the objective
        % function

        def_trueArma = {[0.8,-0.4]};
        def_addNoise = false;
        if nargin == 1 || isempty(varargin{1})
            trueArma = def_trueArma;
        else
            trueArma = varargin{1};
        end
        if nargin < 3 || isempty(varargin{2})
            addNoise = def_addNoise;
        else
            addNoise = varargin{2};
        end

        % generate model - it's fast
        %trueArma = {[0.9,-0.2] [0.3,-0.7],[0.9,0.6],[0.6,0.9]};
        nArma = length(trueArma);
        wnSigma = 0.05;
        nTimepoints = 50000;
        trajLength = 100;
        nTraj = ceil(nTimepoints/trajLength);

        % remember random state for reproducibility
        randomState = rand('state');

        % generate arma model and fit
        arma(1:nArma) = struct('traj',[],'fitResults',[],'individual',[]);
        for i=1:nArma
            arma(i).traj = simSetarma([],[],...
                trueArma{i}(1),trueArma{i}(2),wnSigma,nTimepoints);
            % perturb model
            if addNoise
                arma(i).traj = perturbTrajectory(arma(i).traj,0.01,wnSigma);
            end
            % trajOut are the inputs for armaxFitKalman
            tmp(1:nTraj) = struct('observations',[]);
            arma(i).individual = tmp;
            for j = 1:nTraj
                arma(i).individual(j).observations = ...
                    arma(i).traj((j-1)*trajLength+1:min(nTimepoints,j*trajLength),:);
            end

            % get the initial guess
            initialGuess.arParamP0 = inverseLevDurbExpoAR(trueArma{i}(1));
            initialGuess.maParamP0 = inverseLevDurbExpoMA(trueArma{i}(2));

            % fit - overall fit
            arma(i).fitResults = ...
                armaxFitKalman(arma(i).individual,[],initialGuess,'tl');
            % individual fits
            for j=1:nTraj
                arma(i).individual(j).fitResults = ...
                    armaxFitKalman(arma(i).individual(j),[],initialGuess,'tl');
                % complete structure
                arma(i).individual(j).fitResults.numObserve = ...
                    size(arma(i).individual(j).observations,1);
                % name has to contain "_movie"
                arma(i).individual(j).fitResults.name = sprintf('arma_%i_traj_%i_movie',i,j);
                arma(i).individual(j).fitResults.orderLen = [length(trueArma{i}(1)),length(trueArma{i}(2))];
                arma(i).individual(j).fitResults.type = 'Len';
            end
        end

        for i=nArma:-1:1
            for j=nTraj:-1:1
                clusterInput((i-1)*nTraj+j) = arma(i).individual(j).fitResults;
            end
        end
        for i=nArma:-1:1
            for j=nTraj:-1:1
                clusterInput((i-1)*nTraj+j).lenSeries.observations = arma(i).individual(j).observations;
            end
        end
        options.armaThreshold = 0.05;
        options.verbose = 1;
        % save current data
        save(sprintf('testFuzzyGroupArma_t4_%s',nowString))
        
        % add fitted ARMA model to distance matrix
        data = clusterInput;
        dataFull = arma.fitResults;
        dataFull.name = 'all';
        dataFull.orderLen = [1 1];
        dataFull.type = 'Len';
        dataFull.lengthSeries = [];
        data(end+1) = dataFull;

        % calculate distance matrix
        dm=zeros(nData);
        for iData = 2:nData,
            for jData=1:iData-1
                dm(iData,jData) = armaxModelComp(data(iData),...
                    data(jData));
                dm(jData,iData) = dm(iData,jData);
            end,
        end
        
        % plot distance matrix
        fuzzyGroupArma_visualize(data,[],1,[], [],5e-2,[],2);
        
        figure,histogram(dm(end,1:end-1))



        
        


        
    otherwise
        disp(sprintf('test %i not implemented yet',test));
end


function traj = perturbTrajectory(traj,minNoise,noiseRange)
% perturbTrajectory perturbs a trajectory

% perturb trajectory
trajSize =  size(traj(:,1));
traj  = [traj minNoise+noiseRange*rand(trajSize)];
traj(:,1) = traj(:,1) + randn(trajSize).*traj(:,2);


%-----
%%
% for testing
% for i=1:10,
%     [in{i},out{i},ci{i}]=testFuzzyGroupArma(2,[],[],2000,100);
%     saveas(1,regexprep(sprintf('armaCluster0213_%2.0f_all.fig',i),'\s','0')),
%     saveas(3,regexprep(sprintf('armaCluster0213_%2.0f_arma.fig',i),'\s','0')),
%     saveas(6,regexprep(sprintf('armaCluster0213_%2.0f_diss.fig',i),'\s','0')),
% end