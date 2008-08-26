function couplingCoef = makiMakeScatterPlots(jobType,analysisStruct,analysisStructA)
%MAKIMAKESCATTERPLOTS generates scatter plots of various kinetochore motion parameters
%
%SYNOPSIS makiMakeScatterPlots(jobType,analysisStruct)
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

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%% divide movies into the categories:
%(1) metaphase if no unaligned kinetochores
%(2) late prometaphase 1 if unaligned kinetochores are "close" to the plate
%(3) late prometaphase 2 if unaligned kinetochores are "far" from the plate

movieType = ones(numMovies,1);
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

mIndx = find(movieType==1);
lpm1Indx = find(movieType==2);
lpm2Indx = find(movieType==3);

%% collect data for scatter plots

%metaphase and late prometaphase data:

%center position std
plateStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosition(:,2);

%center displacement std
dispStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosChange(:,2);

%sister separation mean
sepMeanM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSeparation(:,1);

%sister separation change std
sepChangeStdM = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.sisterSepChange(:,2);

%center switching time mean
centerPDispIntMean = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosPChangeInterval(:,1);
centerNDispIntMean = analysisStruct.sepDispSpaceTime.Inlier.meanStdMin25P50P75PMax.indcell.centerPosNChangeInterval(:,1);
centerDispIntMeanM = (centerPDispIntMean + centerNDispIntMean) / 2;

%product of center displacement and center switching time
prodDispTimeM = dispStdM .* centerDispIntMeanM / 7.5;

%angle with normal mean
angleNormMeanM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,1);

%angle with normal std
angleNormStdM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angleWithNormal(:,2);

%angular displacement mean
angVelMeanM = analysisStruct.sisterConnection.Inlier.meanStdMin25P50P75PMax.indcell.angularVel(:,1)*7.5;

%remove the 4 outlier movies that have plate thinkness std > 1.1 um
%and the 1 outlier movie with angular displacement > 22 degrees
indxRemove = [find(plateStdM>1.1); find(angVelMeanM>22)];
plateStdM(indxRemove) = NaN;
dispStdM(indxRemove) = NaN;
sepMeanM(indxRemove) = NaN;
sepChangeStdM(indxRemove) = NaN;
centerDispIntMeanM(indxRemove) = NaN;
prodDispTimeM(indxRemove) = NaN;
angleNormMeanM(indxRemove) = NaN;
angleNormStdM(indxRemove) = NaN;
angVelMeanM(indxRemove) = NaN;

%anaphase data:

if ~isempty(analysisStructA)

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

%remove the 1 outlier movie that has a plate thinkness std > 0.6 um
indxRemove = find(plateStdA>0.6);
plateStdA(indxRemove) = NaN;
dispStdA(indxRemove) = NaN;
sepMeanA(indxRemove) = NaN;
sepChangeStdA(indxRemove) = NaN;
centerDispIntMeanA(indxRemove) = NaN;
prodDispTimeA(indxRemove) = NaN;
angleNormMeanA(indxRemove) = NaN;
angleNormStdA(indxRemove) = NaN;
angVelMeanA(indxRemove) = NaN;

%% scatter plots

%figure 1: center displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),dispStdM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),dispStdM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),dispStdM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,dispStdA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement std (um)');
hold off

%figure 2: sister separation vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),sepMeanM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),sepMeanM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),sepMeanM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,sepMeanA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Sister separation mean (um)');
hold off

%figure 3: sister separation change vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),sepChangeStdM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),sepChangeStdM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),sepChangeStdM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,sepChangeStdA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Sister separation change std (um)');
hold off

%figure 4: center displacement switching time vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),centerDispIntMeanM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),centerDispIntMeanM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),centerDispIntMeanM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,centerDispIntMeanA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement switching time (s)');
hold off

%figure 5: displacement switching time * displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),prodDispTimeM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),prodDispTimeM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),prodDispTimeM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,prodDispTimeA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Center displacement x switching time (um)');
hold off

%figure 6: angle vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angleNormMeanM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angleNormMeanM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),angleNormMeanM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angleNormMeanA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Mean angle with normal (degrees)');
hold off

%figure 7: angle std vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angleNormStdM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angleNormStdM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),angleNormStdM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angleNormStdA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Standard deviation of angle with normal (degrees)');
hold off

%figure 8: angular displacement vs. spread
figure, hold on
plot(plateStdM(lpm1Indx),angVelMeanM(lpm1Indx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.5 0.5 0.5])
plot(plateStdM(lpm2Indx),angVelMeanM(lpm2Indx),'k.','MarkerSize',15)
plot(plateStdM(mIndx),angVelMeanM(mIndx),'MarkerSize',15,'Marker','.',...
    'LineStyle','none','Color',[0.4784 0.06275 0.8941]);
if ~isempty(analysisStructA)
    plot(plateStdA,angVelMeanA,'rx','LineWidth',2,'MarkerSize',8);
    legend('LPM near','LPM far','M','A')
else
    legend('LPM near','LPM far','M')
end
xlabel('Center position std (um)');
ylabel('Angular displacement (degrees)');
hold off

%% coupling coefficients

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

%all
couplingCoef.all.spreadDisp = crossCorr(plateStdAll,dispStdAll,0);
couplingCoef.all.spreadSep = crossCorr(plateStdAll,sepMeanAll,0);
couplingCoef.all.spreadSepChange = crossCorr(plateStdAll,sepChangeStdAll,0);
couplingCoef.all.spreadCenterDispInt = crossCorr(plateStdAll,centerDispIntMeanAll,0);
couplingCoef.all.spreadProdDT = crossCorr(plateStdAll,prodDispTimeAll,0);
couplingCoef.all.spreadAngle = crossCorr(plateStdAll,angleNormMeanAll,0);
couplingCoef.all.spreadAngleStd = crossCorr(plateStdAll,angleNormStdAll,0);
couplingCoef.all.spreadAngVel = crossCorr(plateStdAll,angVelMeanAll,0);

%look only in the regime where there are only metaphase and anaphase movies
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

%% ~~~ the end ~~~
