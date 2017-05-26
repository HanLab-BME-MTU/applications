tmp=load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\EB3\track\plustiptrackerio\track\trackResults.mat');

tracksTestSpeed=TracksHandle(tmp.tracksFinal);
tracksTestSpeedStruct=TracksStruct(tmp.tracksFinal);
numTrack=10;
numTrack=length(tracksTestSpeed);
tic;
for ebIdx=1:numTrack
    tracksTestSpeed(ebIdx).x=tracksTestSpeed(ebIdx).z;
end
toc;

tic;
for ebIdx=1:numTrack
    tracksTestSpeed(ebIdx).x=1;
end
toc;
%%
clear testScratch;
testScratch(numTrack)=TracksHandle;
tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).x=tracksTestSpeed(ebIdx).z;
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).y;
    testScratch(ebIdx).startFrame=tracksTestSpeed(ebIdx).startFrame;
    testScratch(ebIdx).endFrame=tracksTestSpeed(ebIdx).endFrame;   
end
toc;
%%
clear testScratch;
testScratch=tracksTestSpeed.copy();
tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).x=tracksTestSpeed(ebIdx).z;
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).y;   
    testScratch(ebIdx).startFrame=tracksTestSpeed(ebIdx).startFrame;
    testScratch(ebIdx).endFrame=tracksTestSpeed(ebIdx).endFrame;   
end
toc;

%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeed(ebIdx).z(end);
    b=tracksTestSpeed(ebIdx).y(end);
    c=tracksTestSpeed(ebIdx).x(end);
end
toc;


%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeed(ebIdx).z(end);
    b=tracksTestSpeed(ebIdx).y(end);
    c=tracksTestSpeed(ebIdx).x(end);
end
toc;


%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeedStruct(ebIdx).z(end);
    b=tracksTestSpeedStruct(ebIdx).y(end);
    c=tracksTestSpeedStruct(ebIdx).x(end);
end
toc;


