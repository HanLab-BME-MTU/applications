tracksTestSpeed=EB3TracksISOREF.copy();
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
testScratch(numTrack)=struct('x',[],'y',[],'f',[]);
tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).x=tracksTestSpeed(ebIdx).z;
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).y;
    testScratch(ebIdx).f=tracksTestSpeed(ebIdx).f;
end
toc;

tic;
for ebIdx=1:numTrack
    testScratch(ebIdx).y=tracksTestSpeed(ebIdx).x;
    testScratch(ebIdx).f=tracksTestSpeed(ebIdx).f;
end
toc;
