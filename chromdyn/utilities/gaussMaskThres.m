function gmask=gaussMaskThres(amplitude,sigma,fSze,cent,thres);
% guassMaskThres	create a gaussian 3D mask where only values above thres are listed
%
%    gmask=gaussMaskThres(amplitude,sigma,fSze,cent,thres);
%
%    INPUT: amplitude  of gaussian
%           sigma      of gauss mask [sigmaX sigmaY sigmaZ]
%           fSze       size of the gauss mask [sizeX sizeY sizeZ]
%           cent       center postion of gaussian in matrix coordinates
%                       (1,1 is at the center of the top left pixel in an
%                       array!)
%           thres      threshold value
%
%    OUTPUT: gmask   two column matrix where col_1 contains the indices and col_2 the values
%                    note: the size fSze is not stored

% c: 6/08/02 dT

% create matrix containing gauss values. Maximum of Gauss (but not
% necessarily of pixels) is amplitude
gaussImg = amplitude * GaussMask3D(sigma,fSze,cent,0,1);

% find index and values over threshold. Can't combine with find, because
% the find command evaluates the expression before returning values.
idx = find(gaussImg(:)>thres);
val = gaussImg(idx);

% assign output
gmask = [idx, val];

% Dom's old code
% part=thres/amplitude;
% reqSze=ceil(2*sigma.*sqrt(-2*log(part)));
% %grow to odd size
% reqSze=reqSze+~rem(reqSze,2);
% cenReqSze=ceil(reqSze/2);
% 
% shift=round(cent);
% gauss=GaussMask3D(sigma,reqSze,cent-shift,0);
% i=find(gauss>part);
% [ix iy iz]=ind2sub(reqSze,i);
% % put into big picture & cut if necessary
% ix=ix+shift(1)-cenReqSze(1);
% iy=iy+shift(2)-cenReqSze(2);
% iz=iz+shift(3)-cenReqSze(3);
% ix=max(ix,1);
% ix=min(ix,fSze(1));
% iy=max(iy,1);
% iy=min(iy,fSze(2));
% iz=max(iz,1);
% iz=min(iz,fSze(3));
% 
% iBig=sub2ind(fSze,ix,iy,iz);
% gmask=[iBig amplitude*gauss(i)];
