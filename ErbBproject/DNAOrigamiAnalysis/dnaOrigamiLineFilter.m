function Masks=dnaOrigamiLineFilter(len,wid, pixsiz, amp)
% dnaOrigamiLineFilter generates a set of mask for detecting DNA origami
% lines. Uses a simple positive block for the line and negative outside of
% the line. 
%
%inputs : 
%       len -> lenght of the line
%       wid -> width of line
%       pixsiz -> pix size of image
%       amp   -> amplitude of the line
%       
%outputs:
%       Masks -> a cell array of masks
%
%
% Jeffrey Werbin
% Harvard Medical School
%
% last update 12/21/12


 Masks = cell([36,1]);
 
 %place a 50% padding around the line
 pad = fix(len/pixsiz); % total padding
 s = fix(pad/2);%shifting to center of log
 
 %approximates line witdh
 wid = floor(wid/pixsiz)+1;
 
 mask = zeros([2*pad,2*wid+1]);
 mask(s+1:s+1+pad,s+1:s+wid) = amp;
 
 ind = mask == amp;
 
 total = sum(sum(mask));
 
 mask(~ind) = -1*total/(sum(sum(~ind)));
 
 %rotate masks by increments of 5 degrees from 0 to 175 
 for i=1:36
     Masks{i}= imrotate(mask,(i-1)*5,'bicubic','crop');
 end
 
 
end


 
 
 
 
 