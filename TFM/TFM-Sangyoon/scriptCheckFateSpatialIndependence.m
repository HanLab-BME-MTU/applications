% scriptCheckFateSpatialIndependce
% load('movieData.mat')
load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/Colocalization/Region1_gc5/data/tracksNA.mat')
MD=MovieData.load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/movieData.mat');
%% separateMatureAdhesionTracks
outputPath = 'trackSeparation';
tracksNA = separateMatureAdhesionTracks(tracksNA, MD, outputPath);
%% checkAdhesionFateSpatialIndependence
r=50;
[popStat,failNeiStat,matureNeiStat] = checkAdhesionFateSpatialIndependence(tracksNA, r, outputPath);
%% plotting
figure, errorbar([mean(popStat.mRatio),mean(failNeiStat.mRatio),mean(matureNeiStat.mRatio)],...
                            [std(popStat.mRatio),std(failNeiStat.mRatio)/sqrt(length(failNeiStat.mRatio)),...
                            std(matureNeiStat.mRatio)/sqrt(length(matureNeiStat.mRatio))])
%% Showing maturing vs. failing adhesion tracks on cell image in different colors
% showing cell channel
matureID = arrayfun(@(x) (x.maturing==1),tracksNA);
failID = arrayfun(@(x) (x.maturing==0),tracksNA);
FAID = arrayfun(@(x) (x.maturing==3),tracksNA);
longNAID = arrayfun(@(x) (x.maturing==2),tracksNA);

matureIDidx = find(matureID);
failIDidx = find(failID);
FAIDidx = find(FAID);
longNAIDidx = find(longNAID);

paxChannel = MD.getChannel(2);
% For focal adhesions (maturing= 3) or nonmaturing long NA(maturing= 2) 
%% Showing itself
load('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM/130429_cell11_100f/ROI/Colocalization/test/data/cropInfo.mat')
f50 = 88;
figure, imshow(paxChannel.loadImage(f50),[])
hold on
for ii=matureIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(tracksNA(ii).startingFrame:f50),tracksNA(ii).yCoord(tracksNA(ii).startingFrame:f50), 'g')
    end
end
for ii=failIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'r')
    end
end
for ii=FAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        try
            plot(tracksNA(ii).adhBoundary{f50}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50}(:,1)+cropInfo(2), 'k')
            plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
        catch
            try
                plot(tracksNA(ii).adhBoundary{f50-1}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50-1}(:,1)+cropInfo(2), 'k')
                plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
            catch
                try
                    plot(tracksNA(ii).adhBoundary{f50-2}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50-2}(:,1)+cropInfo(2), 'k')
                    plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
                catch
                    try
                        plot(tracksNA(ii).adhBoundary{f50-3}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50-3}(:,1)+cropInfo(2), 'k')
                        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
                    catch
                        try
                            plot(tracksNA(ii).adhBoundary{f50-4}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50-4}(:,1)+cropInfo(2), 'k')
                            plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
                        catch
                            plot(tracksNA(ii).adhBoundary{f50-5}(:,2)+cropInfo(1), tracksNA(ii).adhBoundary{f50-5}(:,1)+cropInfo(2), 'k')
                            plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'k')
                        end
                    end
                end
            end
        end
    end
end
for ii=longNAIDidx'
    if tracksNA(ii).startingFrame <=f50 && tracksNA(ii).endingFrame >=f50
        plot(tracksNA(ii).xCoord(1:f50),tracksNA(ii).yCoord(1:f50), 'y')
    end
end