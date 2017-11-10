function h=donutFilter(outerR,innerR)
% Quick and dirty donut filter (edges are not perfect):
fltOuterDisk = fspecial('disk',outerR);
fltInnerDisk = fspecial('disk',innerR);

% find the right scaling:
outerMax=max(fltOuterDisk(:));
innerMax=max(fltInnerDisk(:));

% scale the inner filter:
fltInnerDiskNrm=fltInnerDisk*outerMax/innerMax;

% bring the inner filter to the same dimension as the outer filter:
fltInnerDiskNrmExt=padarray(fltInnerDiskNrm,[outerR-innerR outerR-innerR]);

% subtract the two, to obtain the donut:
fltDonut=fltOuterDisk-fltInnerDiskNrmExt;

% Normalize the filter:
h=fltDonut/sum(fltDonut(:));