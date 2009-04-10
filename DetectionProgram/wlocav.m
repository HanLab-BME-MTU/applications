function [av,sg]=wlocav(x,edge)

% REPLACE THIS... (Sylvain)

% [liy,lix]=size(x);
% 
% av = zeros(liy, lix);
% sg = zeros(liy, lix);
% 
% %------------------------------
% %     START THE ALGORITHM
% %------------------------------
% 
% 
% 
% lmax=lix+2*edge;
% lmay=liy+2*edge;
% % Copy the vector to transform.
% t =zeros(lmay,lmax);
% 
% for jy=1:lmay
%     for jx=1:lmax        % apply continuous boundary conditions
%         if ((jx < edge+1) && (jy < edge+1))
%             t(jy,jx)=x(1,1);
%         elseif (((jx < lix+edge+1) && (jy < edge+1)) && (jx>edge))
%             t(jy,jx)=x(1,jx-edge);
%         elseif ((jx > lix+edge) && (jy < edge+1))
%             t(jy,jx)=x(1,lix);
%         elseif (((jy < liy+edge+1) && (jx < edge+1)) && (jy>edge))
%             t(jy,jx)=x(jy-edge,1);
%         elseif (((jx > lix+edge) && (jy < edge+liy+1)) && (jy>edge))
%             t(jy,jx)=x(jy-edge,lix);
%         elseif (((jy > liy+edge) && (jx < edge+lix+1)) && (jx>edge))
%             t(jy,jx)=x(liy,jx-edge);
% 
%         elseif ((jy > liy+edge) && (jx < edge+1))
%             t(jy,jx)=x(liy,1);
%         elseif ((jy > liy+edge) && (jx > lix+edge))
%             t(jy,jx)=x(liy,lix);
%         else
%             t(jy,jx)=x(jy-edge,jx-edge);
%         end
%     end
% end
% 
% for iy=1+edge:liy+edge
%     for ix=1+edge:lix+edge
%         tm=t(iy-edge:iy+edge,ix-edge:ix+edge);
%         av(iy-edge,ix-edge)=mean(tm(:));
%         sg(iy-edge,ix-edge)=std(tm(:));
%     end
% end

% ...BY THIS:
h = fspecial('average', [2*edge+1, 2*edge+1]);
av = imfilter(x, h, 'replicate');
av2 = imfilter(x .* x, h, 'replicate');
sg = av2 - av.^2;
sg(sg < 0) = 0;
n = numel(h);
sg = sqrt((n / (n - 1)) * sg);
% END MODIF (Sylvain)
