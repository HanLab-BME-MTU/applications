function [ output_args ] = untitled2( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Elevation vs Azimuth
EB3MarkerSize=10;
KinMarkerSize=50;
cmapKin=jet(600);
cmapEB3=summer(150); %cmapEB3=cmapEB3(1:40,:);
temporalWindow=1; %Number of Frames used for integration
%CellIdx={[1:],[1:5],[1:3]};
for MLIdx=1:length(aMovieListArray)
    ML=aMovieListArray(k); 
    for cIdx=CellIdx{MLIdx}
        MD=MovieData.loadMatFile(ML.movieDataFile_{cIdx});
        handles=setupFigure(1,4,'Name',[ conditionName{MLIdx} ' Cell ' num2str(cIdx) ],'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
        plotSpindleSphericalProjection(handles,cumulAzi{MLIdx,cIdx},cumulElev{MLIdx,cIdx},cumulPoleId{MLIdx,cIdx},cumulTimePt{MLIdx,cIdx}, ...
            cumulAziKin{MLIdx,cIdx},cumulElevKin{MLIdx,cIdx},cumulTimePtKin{MLIdx,cIdx},cumulTrackIdKin{MLIdx,cIdx},[0 MD.nFrames_*MD.timeInterval_],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin,[])
        outpurDir=[MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius)];
        mkdir(outpurDir);
        print([outpurDir filesep 'cumulative.eps'],'-depsc');
        mkdir([outpurDir filesep 'png']);
        mkdir([outpurDir filesep 'eps']);
        for t=20:temporalWindow:60%(MD.nFrames_-temporalWindow)
            [handles,~,fhandle]=setupFigure(1,4,'Name',[ conditionName{MLIdx} ' Cell ' num2str(cIdx) ],'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
            set(fhandle,'Visible','off');
            plotSpindleSphericalProjection(handles,cumulAzi{MLIdx,cIdx},cumulElev{MLIdx,cIdx},cumulPoleId{MLIdx,cIdx},cumulTimePt{MLIdx,cIdx}, ... 
                cumulAziKin{MLIdx,cIdx},cumulElevKin{MLIdx,cIdx},cumulTimePtKin{MLIdx,cIdx},cumulTrackIdKin{MLIdx,cIdx}, ... 
                MD.timeInterval_*[t t+temporalWindow],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin, ...
                EB3CatchingId{MLIdx,cIdx});
            
            print([outpurDir filesep 'png' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.png'],'-dpng');
            print([outpurDir filesep 'eps' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.eps'],'-depsc');
        end 
    end
end



end

