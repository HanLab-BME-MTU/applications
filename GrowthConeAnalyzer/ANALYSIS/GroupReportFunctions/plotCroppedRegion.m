function [ output_args ] = plotCroppedRegion(projPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% load the crop overlay
filename  = dir([projPath filesep 'OriginalStack']);
filename = filename(3:end);
if length(filename)>1
    display('Maria you put some other shit in there check it');
end
img = double(imread([projPath filesep 'OriginalStack' filesep filename(1).name]));

imshow(-img,[]);
if exist([projPath filesep 'cropRegion.mat'],'file')~=0;
    
    
    load([projPath filesep 'cropRegion.mat']);
    % load position
    
    x = pos(1);
    x2 = pos(1)+pos(3);
    y = pos(2);
    y2= pos(2)+pos(4);
    line([x,x],[y +y2],'color','r','Linewidth',2);
    line([x,x2],[y,y],'color','r','Linewidth',2);
    line([x,x2],[y2,y2],'color','r','Linewidth',2);
    line([x2,x2],[y,y2],'color','r','Linewidth',2);
    
end
end

