clc;
clear;

imaqInfo = imaqhwinfo;

hwInfo = imaqhwinfo('hamamatsu');

hwInfo.DeviceInfo(1);

vidobj = videoinput('hamamatsu',1);

triggerconfig(vidobj,'manual');

%%

start(vidobj);
for i = 1:1:1;
    [snapshot, metadata] = getsnapshot(vidobj);
    
    figure; imshow(snapshot,[]); colormap jet
end

stop(vidobj)

%%
delete(vidobj)