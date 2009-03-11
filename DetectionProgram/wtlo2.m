function y=wtlo2(x,k)
% WTLO2 is the low pass transformation for either row or column used for
% wavelet decomposition.  
%  Y = WTLO2(X,K), where X is the image to be processed.
%    K is the level of decomposition for wavelet transformation.
%    Y is the processed result (with low pass transformation).
%
%    Reference:  "Olivo-Marin J.C. 2002. Extraction of spots in biological
%    images using multiscale products. Pattern Recognit. 35: 1989–1996."
%
% Last updated:  Shann-Ching Sam Chen, April 8, 2008
%                add description to the function and rewrote the inner loop
%                in a vectorized form to speed up computation    

[liy,lix]=size(x);

%------------------------------
%     START THE ALGORITHM
%------------------------------
edge=2^(k+1);
lmax=lix+2*edge;
% Copy the vector to transform.
tzeros =zeros(1,lmax);
y = zeros(size(x));
ix=[1+edge:lix+edge];

for it=1:liy
    % For every row of the input matrix...
    % (this makes one wavelet transform
    % for each of the rows of the input matrix)

    % These codes were moved outside the for loop to speed up computation, Sam, 04/08/2008
    %     edge=2^(k+1);
    %     lmax=lix+2*edge;
    %     % Copy the vector to transform.
    %     t =zeros(1,lmax);

    t = tzeros;
    for jx=1:lmax   % apply continuous boundary conditions
        if jx < edge+1
            %t(jx)=x(it,lix-edge+jx);
            t(jx)=x(it,1);
        elseif jx > lix+edge
            %t(jx)=x(it,jx-lix-edge);
            t(jx)=x(it,lix);
        else
            t(jx)=x(it,jx-edge);
        end
    end

    %     These codes were moved outside the for loop to speed up computation, Sam, 04/08/2008
    %     sy(it,:)=yl;		       	% Wavelet vector (1 row vector)
    
    %     for ix=1+edge:lix+edge
    %         % Do lowpass filtering ...
    %         % first line below gives the better B3-spline a trous filter
    %         yl(ix-edge)=0.375*t(ix)+0.25*t(ix+2^(k-1))+0.25*t(ix-2^(k-1))+0.0625*t(ix+2^k)+0.0625*t(ix-2^k);
    %         % next line gives simplest (triangle) a trous filter
    %         %yl(ix-edge)=0.5*t(ix)+0.25*t(ix+2^(k-1))+0.25*t(ix-2^(k-1));
    %     end
    
    % Wavelet vector (1 row vector)
    y(it,:)=0.375*t(ix)+0.25*t(ix+2^(k-1))+0.25*t(ix-2^(k-1))+0.0625*t(ix+2^k)+0.0625*t(ix-2^k);
end				% End of the "rows" loop.

%------------------------------
%    END OF THE ALGORITHM
%------------------------------