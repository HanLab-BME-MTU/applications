function testTrackabilityMitosis()
%% Check tracability between early and late mitosis on a few Cells
outputFolder='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/trackability/interpolarSlice';
%% loading timing based movie table
allMovieToAnalyse=readtable('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/phaseProgression/analysis/movieTables/allMovieToAnalyse.xlsx');
blurrPoleCheckedMoviesIdx=(~(allMovieToAnalyse.blurred|allMovieToAnalyse.doubleCell));
blurrPoleCheckedMovies=allMovieToAnalyse(blurrPoleCheckedMoviesIdx,:);
goodAndOKSNRIdx=ismember(allMovieToAnalyse.EB3SNR,'OK')|ismember(allMovieToAnalyse.EB3SNR,'Good');
blurrPoleCheckedMoviesHighSNR=allMovieToAnalyse(goodAndOKSNRIdx&blurrPoleCheckedMoviesIdx,:);

% Process Collection
minutesBinning=[-2 0.5 2.5 10 16];
processOverlayCell=cell(1,length(minutesBinning)-1);
processProj=cell(1,length(minutesBinning)-1);
maxSpeedCell=cell(1,length(minutesBinning)-1);
[counts,binIdx]=histc(blurrPoleCheckedMoviesHighSNR.starts_min_,minutesBinning);
condName=arrayfun(@(m) [num2str(m) ' min'],minutesBinning,'unif',0);


TestType='twoCells';

switch TestType
case 'twoCells'
    conditionsIdx=[1 length(minutesBinning)-1];
    cellIdx=zeros(1,length(minutesBinning)-1);
    cellIdx(1)=2;
    cellIdx(length(minutesBinning)-1)=3;
case 'AllProcessed'
    conditionsIdx=[1:length(minutesBinning)-1];
    cellIdx=arrayfun(@(bIdx) find(binIdx==bIdx)',1:length(minutesBinning)-1);
end

% Selecting number 3 and 2 on early and late
for bIdx=conditionsIdx % note only begining and end.
    condMovieIndex=find(binIdx==bIdx)';
    for mIdx=condMovieIndex(cellIdx(bIdx))
        if(binIdx(mIdx)>0)
            MD=MovieData.loadMatFile(blurrPoleCheckedMoviesHighSNR.analPath{mIdx});
            if(~isempty(MD.getPackage(160)))
                MD.getPackage(160).displayProcessInfo
                processProj=MD.getPackage(160).getProcess(8);
                processDetection=MD.getPackage(160).getProcess(2);
                processBuildRef=MD.getPackage(160).getProcess(4);
                [maxSpeed,processOverlay]=trackabilityAnalysis(MD,processDetection,processProj,processBuildRef);
                processOverlayCell{bIdx}=[processOverlayCell{bIdx} processOverlay];
                maxSpeedCell{bIdx}=[maxSpeedCell{bIdx} maxSpeed];
            end
        end
    end
end
showFrame=5;
figure();
imshow(imread(sprintfPath(processOverlayCell{1}.outFilePaths_{2},showFrame)),'InitialMagnification','fit');
truesize
figure();
imshow(imread(sprintfPath(processOverlayCell{length(minutesBinning)-1}.outFilePaths_{2},showFrame)),'InitialMagnification','fit');
truesize
figure();

subplot(1,2,1);
histogram(vertcat(maxSpeedCell{1}{showFrame}),0:0.01:1.1);
xlabel('maxSpeed');
subplot(1,2,2);
histogram(vertcat(maxSpeedCell{length(minutesBinning)-1}{showFrame}),0:0.01:1.1);
xlabel('maxSpeed');
%% display
printProcMIPArray([num2cell(processOverlayCell{1}) num2cell(processOverlayCell{length(minutesBinning)-1})],...
    [outputFolder filesep 'earlyLateSingle'],'MIPIndex',2,'MIPSize',400,'maxWidth',820,'maxHeight',1200);