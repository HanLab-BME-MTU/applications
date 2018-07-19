tmp=load('C:\Users\Philippe\project-local\externBetzig\analysis\adavid\smallSample\prometaphase\earlyCell1_12\EB3\track\plustiptrackerio\track\trackResults.mat');

tracksTestSpeed=TracksHandle(tmp.tracksFinal);
tracksTestSpeedStruct=TracksStruct(tmp.tracksFinal);

numTrack=1000;

%% Testing Swaping
tic;
for ebIdx=1:numTrack
    tracksTestSpeed(ebIdx).x=tracksTestSpeed(ebIdx).z;
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
end
toc;

%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).x=tracksTestSpeed(ebIdx).z;
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).y;   
end
toc;


%% Testing reading from object. 
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeed(ebIdx).z(1);
    b=tracksTestSpeed(ebIdx).y(1);
    c=tracksTestSpeed(ebIdx).x(1);
end
toc;


%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeedStruct(ebIdx).getX(1);
    b=tracksTestSpeedStruct(ebIdx).getX(1);
    c=tracksTestSpeedStruct(ebIdx).getX(1);
end
toc;


%%
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    a=tracksTestSpeedStruct(ebIdx).z(1);
    b=tracksTestSpeedStruct(ebIdx).y(1);
    c=tracksTestSpeedStruct(ebIdx).x(1);
end
toc;

%% Testing reading from struct instead of object. 
clear testScratch;
testScratch(numTrack)=struct('x',[],'y',[],'z',[],'f',[]);
tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).x=tracksTestSpeed(ebIdx).y;
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).y;
    testScratch(ebIdx).z=tracksTestSpeed(ebIdx).z;
    testScratch(ebIdx).f=tracksTestSpeed(ebIdx).f;
    
end
toc;

tic;
for ebIdx=1:numTrack
    a=tracksTestSpeed(ebIdx).z(end);
    b=tracksTestSpeed(ebIdx).y(end);
    c=tracksTestSpeed(ebIdx).x(end);
end
toc;

tic;
for ebIdx=1:numTrack
    a=testScratch(ebIdx).z(end);
    b=testScratch(ebIdx).y(end);
    c=testScratch(ebIdx).x(end);
end
toc;


