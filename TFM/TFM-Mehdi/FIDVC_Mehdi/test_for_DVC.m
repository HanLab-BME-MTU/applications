% displacements between two volumetric images through digital volume
% correlation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I: cell containing the undeformed, I{1}, and deformed, I{2} 3-D images
%   sSize: interrogation window (subset) size
%   sSpacing: interrogation window (subset) spacing.  Determines window
%             overlap factor
%   ccThreshold: threshold value that defines a bad cross-correlation
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field (u{1:3} = {u_x, u_y, u_z})
%   cc: peak values of the cross-correlation for each interrogation
% 
% NOTES
% -------------------------------------------------------------------------
% To run you need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
% 
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

clear; close all;

I={{1},{2}}; % put your own undeformed(I{1}), and deformed(I{2}) 3D-images
sSize = [128 128 64]; % subset size
sSpacing = sSize/2; 
DVCPadSize = sSpacing/2;
ccThreshold = 1.0000e-04;


I0{1} = (filename{1});
I0{2} = (filename{2});
sizeI0 = size(I0{1});
sizeI = ceil(sizeI0./sSpacing).*sSpacing;
prePad = ceil((sizeI - sizeI0)/2);
postPad = floor((sizeI - sizeI0)/2);
I{1} = padarray(I0{1},prePad,0,'pre');
I{1} = padarray(I{1},postPad,0,'post');

I{2} = padarray(I0{2},prePad,0,'pre');
I{2} = padarray(I{2},postPad,0,'post');

% Estimate displacements via IDVC
[u, cc] = DVC_Mehdi(I,sSize,sSpacing,DVCPadSize,ccThreshold);

save('resultsDVC.mat','u','cc');



