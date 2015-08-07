function []=analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath)
%analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath)
%extracts an average behavior of traction forces from tracks after
%supervised classification.
% e.g.:
% pathForTheMovieDataFile = '/project/cellbiology/gdanuser/adhesion/Sangyoon/Alexia/2015-04-29/cell2/Int4sec';
% outputPath = 'analysis2';
% analyzeForceAfterClassifyingAdhesionTracks(pathForTheMovieDataFile,outputPath)

%% Run colocalizationAdhesionsWithTFM if it's not run
band = 4;
% showAllTracks = false;
% plotEachTrack = false;
% tmaxEach = [];
tmax=2500;
outputFilePath = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
dataPath = [outputFilePath filesep 'data'];
outputFile=strcat(dataPath,filesep,'tracksNA.mat');
% See if you can use existing tracks
if exist(outputFile,'file')
    disp('tracksNA file is found. Using it ... If you do not want to reuse this, please backup the file and rename it to other name than tracksNA.mat.')
    tracksNAFile = load(outputFile,'tracksNA');
    tracksNA = tracksNAFile.tracksNA;
else
    % run analyzeAdhesionMaturation for obtaining tracks from paxillin channel
    tracksNA=colocalizationAdhesionsWithTFM(pathForTheMovieDataFile,band,tmax,false,false,[],...
        'onlyEdge',true,'outputPath',outputPath,'saveAnalysis',true,'matchWithFA',true);
end

%% filter out matured adhesion
pathForColocalization = [pathForTheMovieDataFile filesep 'Colocalization' filesep outputPath];
% load([pathForColocalization filesep 'data'  filesep 'tracksNA.mat'])
[tracksNA, ~]=separateMatureAdhesionTracks(tracksNA);
%% Classify with ampTotal and distToEdge
[idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9]= classifyNascentAdhesionTracks(pathForColocalization,'tracksNA',tracksNA);
 
%% turn all zeros to NaNs in ampTotal
for ii=1:numel(tracksNA)
    tracksNA(ii).ampTotal(tracksNA(ii).ampTotal==0)=NaN;
end
%% Filtering for group1
% lateAmpSlopeG1 = arrayfun(@(x) x.lateAmpSlope,tracksNA(idGroup1));
% idxLateAmpSlopeG1 = lateAmpSlopeG1<0;
% We filter out tracks whose ampTotal is too high, bigger than mean value
% of ampTotal maxima
meanAmpMaximum = mean(arrayfun(@(x) nanmax(x.ampTotal),tracksNA(idGroup1)));
lateAmpTotalG1 = arrayfun(@(x) x.ampTotal(x.endingFrameExtra),tracksNA(idGroup1));
idxLateAmpLow = lateAmpTotalG1<meanAmpMaximum;
idGroup1nl = find(idGroup1); %converting to non-logical index
idGroup1f = idGroup1nl(idxLateAmpLow);

%% drawing group1
epsPath=[pathForColocalization filesep 'eps'];
fileStore = [epsPath filesep 'ampForcePlotG1.eps'];
plotIntensityForce(tracksNA(idGroup1f),fileStore)
%% Filtering for group2
ampSlopeG2 = arrayfun(@(x) x.ampSlope,tracksNA(idGroup2));
idxIncreasingAmpG2 = ampSlopeG2>0;
idGroup2nl = find(idGroup2); %converting to non-logical index
idGroup2f = idGroup2nl(idxIncreasingAmpG2);
%% group 2
fileStoreG2 = [epsPath filesep 'ampForcePlotG2.eps'];
plotIntensityForce(tracksNA(idGroup2f),fileStoreG2)

%% group 3
% lifeTimeG3 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup3));
% [~,longestID]=max(lifeTimeG3);
% % longestID = idGroup3(longestID);
% idxGroup3=find(idGroup3);
% ii=idxGroup3(longestID);
% figure
% plot(1:(tracksNA(ii).lifeTime),tracksNA(ii).amp(tracksNA(ii).presence))
%% group 3 plotting

fileStoreG3 = [epsPath filesep 'ampForcePlotG3.eps'];
plotIntensityForce(tracksNA(idGroup3),fileStoreG3)

%% group4 plotting
fileStoreG4 = [epsPath filesep 'ampForcePlotG4.eps'];
plotIntensityForce(tracksNA(idGroup4),fileStoreG4)
%% group5 plotting
fileStoreG5 = [epsPath filesep 'ampForcePlotG5.eps'];
plotIntensityForce(tracksNA(idGroup5),fileStoreG5)
%% group6 plotting
fileStoreG6 = [epsPath filesep 'ampForcePlotG6.eps'];
plotIntensityForce(tracksNA(idGroup6),fileStoreG6)
%% group7 plotting
fileStoreG7 = [epsPath filesep 'ampForcePlotG7.eps'];
plotIntensityForce(tracksNA(idGroup7),fileStoreG7)
%% group8 plotting
fileStoreG8 = [epsPath filesep 'ampForcePlotG8.eps'];
plotIntensityForce(tracksNA(idGroup8),fileStoreG8)
%% group9 plotting
fileStoreG9 = [epsPath filesep 'ampForcePlotG9.eps'];
plotIntensityForce(tracksNA(idGroup9),fileStoreG9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For group 1,3,7, look at the cross-correlation score, and collect those with high scores and see if they represent adhesion/edge protrusion  behavior
% group 1
ccG1 = arrayfun(@(x) x.CCscore,tracksNA(idGroup1f));
% if we pick up adhesions wth high correlation, there should be an
% increasing phase in the force too...
idxHighCCG1=ccG1>0.5;
% delete very short life time
LTG1 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup1f));
idxLongLTG1=LTG1>20;
fileStoreG1_hiCC = [epsPath filesep 'ampForcePlotG1_hiCC.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCC)
%% those with intermediate cc
idxIntmedCCG1=ccG1>0.1 & ccG1<=0.5;
% delete very short life time
fileStoreG1_middleCC = [epsPath filesep 'ampForcePlotG1_middleCC.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_middleCC)
%% those with low and negative cc
idxLowCCG1=ccG1<0.1;
% delete very short life time
fileStoreG1_lowCC = [epsPath filesep 'ampForcePlotG1_lowCC.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_lowCC)
%% Which group has tracks with high cc-score?
ccAll = arrayfun(@(x) x.CCscore,tracksNA);
hiCCNAs = ccAll>0.5;
idGroups = {idGroup1,idGroup2,idGroup3,idGroup4,idGroup5,idGroup6,idGroup7,idGroup8,idGroup9};
numhiCC = zeros(9,1);
ratiohiCC = zeros(9,1);
numTracksPerGroups = zeros(9,1);
for ii=1:9
    numTracksPerGroups(ii) = sum(idGroups{ii});
    numhiCC(ii) = sum(idGroups{ii} & hiCCNAs);
    ratiohiCC(ii) = numhiCC(ii) /  numTracksPerGroups(ii);
end
[~,indMaxCC]=max(ratiohiCC);
disp(['A group containing the most ratio of high CC score: group ' num2str(indMaxCC) ' with ratio of ' num2str(ratiohiCC(indMaxCC)) '.'])
iGroups = (1:9)';

ratioT = table(iGroups,numhiCC,numTracksPerGroups,ratiohiCC);
disp(ratioT)
writetable(ratioT,[pathForColocalization filesep 'data' filesep 'ratioTracksHiCCscore.csv'])
%% Event alignment for Group1
% since there is variance in peaks of vinculin, the rise-and-fall behavior
% is not picked up in the average. So I'll align these event with the
% maximum of the intensity timeseries
fileStoreG1_hiCCshifted = [epsPath filesep 'ampForcePlotG1_hiCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCCshifted,1)

fileStoreG1_middleCCshifted = [epsPath filesep 'ampForcePlotG1_middleCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_middleCCshifted,1)

fileStoreG1_lowCCshifted = [epsPath filesep 'ampForcePlotG1_lowCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_lowCCshifted,1)

fileStoreG1_Allshifted = [epsPath filesep 'ampForcePlotG1_All_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxLongLTG1)),fileStoreG1_Allshifted,1)
%% Event alignment for Group3
% since there is variance in peaks of vinculin, the rise-and-fall behavior
% is not picked up in the average. So I'll align these event with the
% maximum of the intensity timeseries
ccG3 = arrayfun(@(x) x.CCscore,tracksNA(idGroup3));
LTG3 = arrayfun(@(x) x.lifeTime,tracksNA(idGroup3));
meanForceG3 = arrayfun(@(x) mean(x.forceMag),tracksNA(idGroup3));
idxLongLTG3=LTG3>20;
idxHighCCG3=ccG3>0.5;
idxIntmedCCG3=ccG3>0.1 & ccG3<=0.5;
idxLowCCG3=ccG3<0.1;
lowForceG3 = meanForceG3<20;
idxGroup3 = find(idGroup3);
fileStoreG3_hiCCshifted = [epsPath filesep 'ampForcePlotG3_hiCC_shifted.eps'];
plotIntensityForce(tracksNA(idxGroup3(idxHighCCG3 & idxLongLTG3)),fileStoreG3_hiCCshifted,1)

fileStoreG3_hiCCshiftedLowForce = [epsPath filesep 'ampForcePlotG3_hiCC_shifted_lowForce.eps'];
plotIntensityForce(tracksNA(idxGroup3(idxHighCCG3 & idxLongLTG3 & lowForceG3)),fileStoreG3_hiCCshiftedLowForce,1)

fileStoreG3_middleCCshifted = [epsPath filesep 'ampForcePlotG3_middleCC_shifted.eps'];
plotIntensityForce(tracksNA(idxGroup3(idxIntmedCCG3 & idxLongLTG3)),fileStoreG3_middleCCshifted,1)

fileStoreG3_lowCCshifted = [epsPath filesep 'ampForcePlotG3_lowCC_shifted.eps'];
plotIntensityForce(tracksNA(idxGroup3(idxLowCCG3 & idxLongLTG3)),fileStoreG3_lowCCshifted,1)

fileStoreG3_Allshifted = [epsPath filesep 'ampForcePlotG3_All_shifted.eps'];
plotIntensityForce(tracksNA(idxGroup3(idxLongLTG3)),fileStoreG3_Allshifted,1)


%% Where are they?
% Now we have to make a function that shows in the last vinculin image, the
% tracks of interest with labels 
imgMap = load([pathForColocalization filesep 'pax'  filesep 'paxImgStack.mat'],'paxImgStack');
imgMap = imgMap.paxImgStack;
colors = distinguishable_colors(3,'k');

figure, imshow(imgMap(:,:,end),[])
hold on
htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)));
htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)));
htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),'UniformOutput',false);
arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)));
legend([htrackG1{1} htrackG2{1} htrackG3{1}],{'G1 with high CC','G1 with intermediate CC','G1 with low CC'},'TextColor','w','Location','best')
legend('boxoff')
print('-depsc2', '-r300', [pathForColocalization filesep 'eps' filesep 'FluorescenceChannelWithTracksInGroup1.eps']);
%% See if there is any time lag
avgTimeLag = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLongLTG1))));
avgTimeLagHiCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1))));
stdTimeLagHiCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1))));
avgTimeLagMidCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1))));
stdTimeLagMidCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1))));
avgTimeLagLoCC = mean(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1))));
stdTimeLagLoCC = std(arrayfun(@(x) x.CCmaxLag,tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1))));
%% Plot edge protrusion distance
figure, hold all
fileStoreG1_hiCC_edge = [epsPath filesep 'edgeDistPlotG1_hiCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxHighCCG1 & idxLongLTG1)),fileStoreG1_hiCC_edge,true,'Source', 'edgeAdvanceDist')
fileStoreG1_midCC_edge = [epsPath filesep 'edgeDistPlotG1_midCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxIntmedCCG1 & idxLongLTG1)),fileStoreG1_midCC_edge,true,'Source', 'edgeAdvanceDist')
fileStoreG1_loCC_edge = [epsPath filesep 'edgeDistPlotG1_loCC_shifted.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_loCC_edge,true,'Source', 'edgeAdvanceDist')
fileStoreG1_loCC_edgeNoshift = [epsPath filesep 'edgeDistPlotG1_loCC.eps'];
plotIntensityForce(tracksNA(idGroup1f(idxLowCCG1 & idxLongLTG1)),fileStoreG1_loCC_edgeNoshift,false,'Source', 'edgeAdvanceDist')
% ii=idGroup1f(4);
% fileStoreG1_1 = [epsPath filesep 'ampForcePlotG1_' num2str(ii) '.eps'];
% plotIntensityForce(tracksNA(ii),fileStoreG1_1)
save([pathForColocalization filesep 'data' filesep 'allDataAfterClassification.mat'])

end
