function [Masks]=fourEdgeMask(pixsiz, sigma)
%fourEdgeMask generates a set of masks for finding features that look like the 4edge
%DNA orgami pattern (100 nm by 70 nm). 
%pixsiz -> effective pixel size of the image
%sigma  -> size of guassian psf
%
% Jeffrey Werbin
% Harvard Medical School
%
% last update: 12/13/12

 Masks = cell([18,1]);
 pad = 5; % size of 'log' mask
 s = fix(pad/2);%shifting to center of log
 
 % list of dye positions in nm delta X = 62.4 delta Y = 41.6
 pos = [ 15,15; 15, 56.6; 77.4, 15; 77.4,56.6];
 pos = round(pos/pixsiz+pad); %converts to coordinates of mask
 
 %creates a centered lapalcian of gaussian
 h = fspecial('log',pad,sigma/pixsiz);
 
 % create a field approximately the size od the origami
 x= ceil(100.0/pixsiz);
 y= ceil(70.0/pixsiz);
 m = zeros([x+2*pad,y+2*pad]);
 
 for i = 1:numel(pos(:,1))
     m(pos(i,1)-s:pos(i,1)-s+pad-1,pos(i,2)-s:pos(i,2)-s+pad-1) = m(pos(i,1)-s:pos(i,1)-s+pad-1,pos(i,2)-s:pos(i,2)-s+pad-1)+h;
 end
 
 %crops the mask out of the padding
 %note this crops out by an extra pixel
 % this should only change where the maxium response pixel is relative to
 % the actual feature
 
 mask = m(pad:pad+1+x,pad:pad+1+y);
 
 %rotate masks by increments of 5 degrees from 0 to 85 
 for i=1:18
     Masks{i}= imrotate(mask,(i-1)*5,'bicubic','crop');
 end



end