function [couplingCoef,lineFit] = makiMakeScatterPlots(jobType,analysisStruct,...
    analysisStructA,distingPhase)
%MAKIMAKESCATTERPLOTS generates scatter plots of various kinetochore motion parameters
%
%SYNOPSIS [couplingCoef,lineFit] = makiMakeScatterPlots(jobType,analysisStruct,...
%    analysisStructA,distingPhase)
%
%INPUT  jobType: string which can take the values:
%               'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%               'MCAINSH'
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected
%                          movies. First column: file name, second column: file path.
%               .sepDispSpaceTime: Field resulting from running
%                          makiSepDispSpaceTime.
%                       Optional. If not input, GUI to load movies is launched.
%       analysisStructA: Same as analysisStruct but for movies that enter
%                        anaphase. Optional. Default: [].
%       distingPhase: 1 to distinguish between pre-anaphase phases, 0
%                     otherwise. Optional. default: 1.
%
%OUTPUT couplingCoef: Coupling coefficients for the different variables
%
%REMARKS Code assumes 7.5-sec sampling
%
%Khuloud Jaqaman, August 2008

%% input
if nargin < 3 || isempty(analysisStructA)
    analysisStructA = [];
end

if nargin < 4 || isempty(distingPhase)
    distingPhase = 1;
end

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
if distingPhase
    for iMovie = numMovies : -1 : 1
        dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
    end
end

%% divide pre-anaphase movies into the categories:
%(1) metaphase if no unaligned kinetochores
%(2) late prometaphase 1 if unaligned kinetochores are "close" to the plate
%(3) late prometaphase 2 if unaligned kinetochores are "far" from the plate

movieType = ones(numMovies,1);
if distingPhase
    for iMovie = 1 : numMovies

        %get frame phases
        moviePhases = vertcat(dataStruct(iMovie).updatedClass.phase);

        %if any of the phases are 'p', i.e. if there are unaligned kinetochores
        if any(moviePhases == 'p')

            %get the plate standard deviation for this movie
            plateStd = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosition(iMovie,2);

            %find which tracks are unaligned
            tracks = dataStruct(iMovie).tracks(dataStruct(iMovie).updatedClass(1).tracksUnaligned);

            %keep only tracks longer than 5 frames
            criteria.lifeTime.min = 5;
            tracks = tracks(chooseTracks(tracks,criteria));

            if ~isempty(tracks)

                %convert tracks from structure to matrix format
                [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracks);

                %get number of unaligned tracks and number of frames
                [numTracks,numFrames] = size(trackedFeatureIndx);

                %go over all unaligned tracks
                kinPos = NaN(numTracks,1);
                for iKin = 1 : numTracks

                    %find feature indices making up track
                    kinIndx = trackedFeatureIndx(iKin,:)';
                    kinIndx(kinIndx==0) = NaN;

                    %get aligned kinetochore coordinates
                    coordsKin = NaN(numFrames,6);
                    for iFrame = 1 : numFrames
                        if ~isnan(kinIndx(iFrame))
                            coordsKin(iFrame,:) = dataStruct(iMovie).frameAlignment(iFrame).alignedCoord(kinIndx(iFrame),:);
                        end
                    end

                    %calculate average position of kinetochore along the normal
                    %to the metaphase plate
                    kinPos(iKin) = nanmean(coordsKin(:,1)); %um

                end

                %get the furthest unaligned kinetochore position
                kinPosMax = max(abs(kinPos));

                %divide further kinetochore position by plate standard deviation
                maxPosStdRatio = kinPosMax / plateStd;

                %categorise movie based on this ratio
                if maxPosStdRatio < 10
                    movieType(iMovie) = 2;
                else
                    movieType(iMovie) = 3;
                end

            end %(if ~isempty(tracks))

        end %(if any(strcmp(moviePhases,'p')))

    end %(for iMovie = 1 : numMovies)
end %(if distingPhase)

mIndx = find(movieType==1);
lpm1Indx = find(movieType==2);
lpm2Indx = find(movieType==3);

%% collect data for scatter plots

%metaphase and late prometaphase data:

%center position std
plateStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosition(:,2,:);
plateStdM = plateStdM(:);

%center displacement std
dispStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosChange(:,2,:);
dispStdM = dispStdM(:);

%sister separation mean
sepMeanM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSeparation(:,1,:);
sepMeanM = sepMeanM(:);

%sister separation change std
sepChangeStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSepChange(:,2,:);
sepChangeStdM = sepChangeStdM(:);

%center switching time mean
centerPDispIntMean = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosPChangeInterval(:,1,:);
centerNDispIntMean = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosNChangeInterval(:,1,:);
centerDispIntMeanM = (centerPDispIntMean + centerNDispIntMean) / 2;
centerDispIntMeanM = centerDispIntMeanM(:);
% centerDispIntMeanM = rand(size(sepChangeStdM));

%product of center displacement and center switching time
prodDispTimeM = dispStdM .* centerDispIntMeanM / 7.5;

%angle with normal mean
angleNormMeanM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,1,:);
angleNormMeanM = angleNormMeanM(:);
% angleNormMeanM = rand(size(sepChangeStdM));

%angle with normal std
angleNormStdM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,2,:);
angleNormStdM = angleNormStdM(:);
% angleNormStdM = rand(size(sepChangeStdM));

%angular displacement mean
angVelMeanM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angularVel(:,1,:)*7.5;
angVelMeanM = angVelMeanM(:);
% angVelMeanM = rand(size(sepChangeStdM));

%anaphase data:

if isempty(analysisStructA)

    plateStdA = [];
    dispStdA = [];
    sepMeanA = [];
    sepChangeStdA = [];
    centerDispIntMeanA = [];
    prodDispTimeA = [];
    angleNormMeanA = [];
    angleNormStdA = [];
    angVelMeanA = [];
    
else
    
    %center position std
    plateStdA = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosition(:,2);

    %center displacement std
    dispStdA = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosChange(:,2);

    %sister separation mean
    sepMeanA = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSeparation(:,1);

    %sister separation change std
    sepChangeStdA = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSepChange(:,2);

    %center switching time mean
    centerPDispIntMean = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosPChangeInterval(:,1);
    centerNDispIntMean = analysisStructA.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosNChangeInterval(:,1);
    centerDispIntMeanA = (centerPDispIntMean + centerNDispIntMean) / 2;

    %product of center displacement and center switching time
    prodDispTimeA = dispStdA .* centerDispIntMeanA / 7.5;

    %angle with normal mean
    angleNormMeanA = analysisStructA.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,1);

    %angle with normal std
    angleNormStdA = analysisStructA.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,2);
    
    %angular displacement mean
    angVelMeanA = analysisStructA.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angularVel(:,1)*7.5;

end

%% line fit and coupling coefficients

%put all data together
plateStdAll = [plateStdM;plateStdA];
dispStdAll = [dispStdM;dispStdA];
sepMeanAll = [sepMeanM;sepMeanA];
sepChangeStdAll = [sepChangeStdM;sepChangeStdA];
centerDispIntMeanAll = [centerDispIntMeanM;centerDispIntMeanA];
prodDispTimeAll = [prodDispTimeM;prodDispTimeA];
angleNormMeanAll = [angleNormMeanM;angleNormMeanA];
angleNormStdAll = [angleNormStdM;angleNormStdA];
angVelMeanAll = [angVelMeanM;angVelMeanA];

%fit each scatter plot with 1 and 2 lines, and decide which fits better
indxKeep = find(~isnan(plateStdAll) & ~isnan(dispStdAll) & plateStdAll~=0 & dispStdAll~=0);
lineFit.spreadDisp = fit1Line2Lines(plateStdAll(indxKeep),dispStdAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(sepMeanAll) & plateStdAll~=0 & sepMeanAll~=0);
lineFit.spreadSep = fit1Line2Lines(plateStdAll(indxKeep),sepMeanAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(sepChangeStdAll) & plateStdAll~=0 & sepChangeStdAll~=0);
lineFit.spreadSepChange = fit1Line2Lines(plateStdAll(indxKeep),sepChangeStdAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(centerDispIntMeanAll) & plateStdAll~=0 & centerDispIntMeanAll~=0);
lineFit.spreadCenterDispInt = fit1Line2Lines(plateStdAll(indxKeep),centerDispIntMeanAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(prodDispTimeAll) & plateStdAll~=0 & prodDispTimeAll~=0);
lineFit.spreadProdDT = fit1Line2Lines(plateStdAll(indxKeep),prodDispTimeAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(angleNormMeanAll) & plateStdAll~=0 & angleNormMeanAll~=0);
lineFit.spreadAngle = fit1Line2Lines(plateStdAll(indxKeep),angleNormMeanAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(angleNormStdAll) & plateStdAll~=0 & angleNormStdAll~=0);
lineFit.spreadAngleStd = fit1Line2Lines(plateStdAll(indxKeep),angleNormStdAll(indxKeep),0.3);
indxKeep = find(~isnan(plateStdAll) & ~isnan(angVelMeanAll) & plateStdAll~=0 & angVelMeanAll~=0);
lineFit.spreadAngVel = fit1Line2Lines(plateStdAll(indxKeep),angVelMeanAll(indxKeep),0.3);

%calculate overall coupling coefficient
couplingCoef.all.spreadDisp = crossCorr(plateStdAll,dispStdAll,0);
couplingCoef.all.spreadSep = crossCorr(plateStdAll,sepMeanAll,0);
couplingCoef.all.spreadSepChange = crossCorr(plateStdAll,sepChangeStdAll,0);
couplingCoef.all.spreadCenterDispInt = crossCorr(plateStdAll,centerDispIntMeanAll,0);
couplingCoef.all.spreadProdDT = crossCorr(plateStdAll,prodDispTimeAll,0);
couplingCoef.all.spreadAngle = crossCorr(plateStdAll,angleNormMeanAll,0);
couplingCoef.all.spreadAngleStd = crossCorr(plateStdAll,angleNormStdAll,0);
couplingCoef.all.spreadAngVel = crossCorr(plateStdAll,angVelMeanAll,0);

%calculate coupling coefficient in the regime where metaphase and anaphase movies are
%max plate width = 0.5
indxMA = find(plateStdAll<0.5);
couplingCoef.maRegime.spreadDisp = crossCorr(plateStdAll(indxMA),dispStdAll(indxMA),0);
couplingCoef.maRegime.spreadSep = crossCorr(plateStdAll(indxMA),sepMeanAll(indxMA),0);
couplingCoef.maRegime.spreadSepChange = crossCorr(plateStdAll(indxMA),sepChangeStdAll(indxMA),0);
couplingCoef.maRegime.spreadCenterDispInt = crossCorr(plateStdAll(indxMA),centerDispIntMeanAll(indxMA),0);
couplingCoef.maRegime.spreadProdDT = crossCorr(plateStdAll(indxMA),prodDispTimeAll(indxMA),0);
couplingCoef.maRegime.spreadAngle = crossCorr(plateStdAll(indxMA),angleNormMeanAll(indxMA),0);
couplingCoef.maRegime.spreadAngleStd = crossCorr(plateStdAll(indxMA),angleNormStdAll(indxMA),0);
couplingCoef.maRegime.spreadAngVel = crossCorr(plateStdAll(indxMA),angVelMeanAll(indxMA),0);

%calculate the coupling coefficient in the two regimes as determined by the
%line fit

plateStdDiv = lineFit.spreadDisp.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadDisp = ...
    [crossCorr(plateStdAll(indxBelow),dispStdAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),dispStdAll(indxAbove),0)];

plateStdDiv = lineFit.spreadSep.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadSep = ...
    [crossCorr(plateStdAll(indxBelow),sepMeanAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),sepMeanAll(indxAbove),0)];

plateStdDiv = lineFit.spreadSepChange.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadSepChange = ...
    [crossCorr(plateStdAll(indxBelow),sepChangeStdAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),sepChangeStdAll(indxAbove),0)];

plateStdDiv = lineFit.spreadCenterDispInt.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadCenterDispInt = ...
    [crossCorr(plateStdAll(indxBelow),centerDispIntMeanAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),centerDispIntMeanAll(indxAbove),0)];

plateStdDiv = lineFit.spreadProdDT.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadProdDT = ...
    [crossCorr(plateStdAll(indxBelow),prodDispTimeAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),prodDispTimeAll(indxAbove),0)];

plateStdDiv = lineFit.spreadAngle.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadAngle = ...
    [crossCorr(plateStdAll(indxBelow),angleNormMeanAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),angleNormMeanAll(indxAbove),0)];

plateStdDiv = lineFit.spreadAngleStd.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadAngle = ...
    [crossCorr(plateStdAll(indxBelow),angleNormStdAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),angleNormStdAll(indxAbove),0)];

plateStdDiv = lineFit.spreadAngVel.line2.param(5);
indxBelow = find(plateStdAll<=plateStdDiv);
indxAbove = find(plateStdAll>plateStdDiv);
couplingCoef.lineRegime.spreadAngVel = ...
    [crossCorr(plateStdAll(indxBelow),angVelMeanAll(indxBelow),0) ...
    crossCorr(plateStdAll(indxAbove),angVelMeanAll(indxAbove),0)];

% % % %determine regimes by maximizing the difference between the coupling
% % % %coefficients
% % % couplingCoef.selfRegime.spreadDisp = couplCoef2Regimes(plateStdAll,dispStdAll);
% % % couplingCoef.selfRegime.spreadSep = couplCoef2Regimes(plateStdAll,sepMeanAll);
% % % couplingCoef.selfRegime.spreadSepChange = couplCoef2Regimes(plateStdAll,sepChangeStdAll);
% % % couplingCoef.selfRegime.spreadCenterDispInt = couplCoef2Regimes(plateStdAll,centerDispIntMeanAll);
% % % couplingCoef.selfRegime.spreadProdDT = couplCoef2Regimes(plateStdAll,prodDispTimeAll);
% % % couplingCoef.selfRegime.spreadAngle = couplCoef2Regimes(plateStdAll,angleNormMeanAll);
% % % couplingCoef.selfRegime.spreadAngleStd = couplCoef2Regimes(plateStdAll,angleNormStdAll);
% % % couplingCoef.selfRegime.spreadAngVel = couplCoef2Regimes(plateStdAll,angVelMeanAll);

%% scatter plots

%figure 1: center displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),dispStdM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),dispStdM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),dispStdM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,dispStdA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement std (um)');
if lineFit.spreadDisp.lines1or2 == 1
    lineParam = lineFit.spreadDisp.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadDisp.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 2: sister separation vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),sepMeanM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),sepMeanM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),sepMeanM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,sepMeanA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Sister separation mean (um)');
if lineFit.spreadSep.lines1or2 == 1
    lineParam = lineFit.spreadSep.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadSep.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 3: sister separation change vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),sepChangeStdM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),sepChangeStdM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),sepChangeStdM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,sepChangeStdA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Sister separation change std (um)');
if lineFit.spreadSepChange.lines1or2 == 1
    lineParam = lineFit.spreadSepChange.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadSepChange.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 4: center displacement switching time vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),centerDispIntMeanM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),centerDispIntMeanM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),centerDispIntMeanM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,centerDispIntMeanA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement switching time (s)');
if lineFit.spreadCenterDispInt.lines1or2 == 1
    lineParam = lineFit.spreadCenterDispInt.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadCenterDispInt.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 5: displacement switching time * displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),prodDispTimeM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),prodDispTimeM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),prodDispTimeM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,prodDispTimeA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement x switching time (um)');
if lineFit.spreadProdDT.lines1or2 == 1
    lineParam = lineFit.spreadProdDT.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadProdDT.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 6: angle vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angleNormMeanM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angleNormMeanM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),angleNormMeanM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angleNormMeanA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Mean angle with normal (degrees)');
if lineFit.spreadAngle.lines1or2 == 1
    lineParam = lineFit.spreadAngle.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadAngle.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 7: angle std vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angleNormStdM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angleNormStdM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),angleNormStdM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angleNormStdA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Standard deviation of angle with normal (degrees)');
if lineFit.spreadAngleStd.lines1or2 == 1
    lineParam = lineFit.spreadAngleStd.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadAngleStd.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off

%figure 8: angular displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angVelMeanM(lpm1Indx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angVelMeanM(lpm2Indx),'ko','MarkerSize',10,'LineWidth',1)
plot(plateStdM(mIndx),angVelMeanM(mIndx),'MarkerSize',10,'Marker','o',...
    'LineWidth',1,'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angVelMeanA,'rx','LineWidth',1,'MarkerSize',15);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Angular displacement (degrees)');
if lineFit.spreadAngVel.lines1or2 == 1
    lineParam = lineFit.spreadAngVel.line1.param;
    xvalues = [min([plateStdM;plateStdA]) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
else
    lineParam = lineFit.spreadAngVel.line2.param;
    xvalues = [min([plateStdM;plateStdA]) lineParam(5) max([plateStdM;plateStdA])];
    yvalues = lineParam(1) * xvalues + lineParam(2);
    yvalues(3) = lineParam(3) * xvalues(3) + lineParam(4);
end
plot(xvalues,yvalues,'k','LineWidth',1)
hold off


%% subfunctions

function results = couplCoef2Regimes(xData,yData)

%initial guess
initialGuess = 0.3;
options = optimset('Display','on');

%get divider
divider = fminunc(@couplCoefDiff,initialGuess,options,xData,yData);

%calculate coupling coefficients in the two regimes
couplCoefBef = crossCorr(xData(xData<=divider),yData(xData<=divider),0);
couplCoefAft = crossCorr(xData(xData>divider),yData(xData>divider),0);

%assign results
results = [divider couplCoefBef couplCoefAft];

function coefDiff = couplCoefDiff(divider,xData,yData)

%calculate coupling coefficients in the two regimes
couplCoefBef = crossCorr(xData(xData<=divider),yData(xData<=divider),0);
couplCoefAft = crossCorr(xData(xData>divider),yData(xData>divider),0);

%calculate the absolute difference and give it a minus sign
coefDiff = abs(couplCoefAft(1)/couplCoefBef(1));

%% ~~~ the end ~~~
