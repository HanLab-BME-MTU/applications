function []=checkDetection(path_image,candsOLD,candsNEW)

%read in image:
if nargin <1 || isempty(path_image)
    [ref_filename, ref_pathname] = uigetfile({'*.TIF';'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select the image to be overlayed with cands');
     path_image=[ ref_pathname filesep ref_filename];
end

I = double(imread(path_image));

imagesc(I);
colormap('gray')
hold on
for i=1:length(candsOLD)
    plot(candsOLD(i).Lmax(2),candsOLD(i).Lmax(1),'ob','MarkerSize',2)
end
for i=1:length(candsNEW)
    plot(candsNEW(i).Lmax(2),candsNEW(i).Lmax(1),'.r','MarkerSize',2)
end
hold off