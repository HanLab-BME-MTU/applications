function focalAdhesionDetector(I, BW, sigmaPSF)

% Unknown parameters:
% type:    nescent/mature(elongated) focal adhesion
% n:       number of focal adhesions
% theta:   orientation of focal adhesion
% l:       length of focal adhesion
% w:       width of focal adhesion
% A:       average amplitude along focal adhesion
% xc,yc:   center of focal adhesion

% Known or assumed:
% width:   the width of every focal adhesion is assumed to be the same:
%          Each focal adhesion are diffraction limited structure where the
%          width = 2 * sigmaPSF
%
% theta:   we assumed theta constant within a single pack of focal
%          adhesions.

I = double(I);

% filtering
[~, T, NMS] = steerableFiltering(I, 'UnserM2', sigmaPSF);

NMS(NMS < 0) = 0;

idx = find(NMS & BW);
[y x] = ind2sub(size(NMS), idx);
theta = T(idx);

imshow(I, []); hold on;
quiver(x, y, cos(theta), sin(theta));
