function [sxx, syy, sxy, svm, meta] = constructInternalStressFromTraction(FF, outDir, varargin)
% constructInternalStressFromTraction
% Monolayer Stress Microscopy (MSM) using least-norm Fourier solution.
%
% INPUT
%   FF(frame).pos : [2xN] or [Nx2] grid coordinates
%   FF(frame).vec : [2xN] or [Nx2] traction vectors (Tx,Ty) in Pa
%
% OUTPUT (units)
%   sxx, syy, sxy : stress resultants (Pa*m = N/m)
%
% OPTIONS (name/value)
%   'frameIdx'        : which frame to plot (default 1)
%   'pixelSizeNm'     : default 325 nm/pixel
%   'posUnits'        : 'pixel' | 'um' | 'm'  (default 'pixel')
%   'lam'             : regularization parameter (default 1e-2)
%   'savePrefix'      : filename prefix (default 'msm')
%
% NOTE
%   Governing equations:
%     d/dx sxx + d/dy sxy = -Tx
%     d/dx sxy + d/dy syy = -Ty
%   solved in Fourier domain with minimum-norm regularization.
% % EXAMPLE
% load('forceField.mat');  % gives variable forceField
% outDir = pwd;
% 
% [sxx, syy, sxy] = constructInternalStressFromTraction(forceField, outDir, ...
%     'frameIdx', 1, ...
%     'posUnits', 'pixel', ...
%     'pixelSizeNm', 325, ...
%     'lam', 1e-2, ...
%     'savePrefix', 'Fig10_MSM');

if nargin < 2 || isempty(outDir)
    outDir = pwd;
end
if ~exist(outDir,'dir'); mkdir(outDir); end

% Defaults
frameIdx    = 1;
pixelSizeNm = 325;
posUnits    = 'pixel';   % IMPORTANT: change to 'um' or 'm' if your pos is not pixels
lam         = 1e-2;
savePrefix  = 'msm';

% Parse varargin
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'frameidx',    frameIdx    = varargin{k+1};
        case 'pixelsizenm', pixelSizeNm = varargin{k+1};
        case 'posunits',    posUnits    = lower(varargin{k+1});
        case 'lam',         lam         = varargin{k+1};
        case 'saveprefix',  savePrefix  = varargin{k+1};
        otherwise
            error('Unknown option: %s', varargin{k});
    end
end

% --- Extract and standardize ---
pos = FF(frameIdx).pos;
vec = FF(frameIdx).vec;

[posXY, vecXY] = standardizePosVec(pos, vec); % Nx2

x = posXY(:,1);
y = posXY(:,2);
tx = vecXY(:,1);  % Pa
ty = vecXY(:,2);  % Pa

ux = unique(x); uy = unique(y);
nx = numel(ux); ny = numel(uy);
assert(nx*ny == numel(x), 'Grid size mismatch: nx*ny != N. Check pos/vec formatting.');

dx = median(diff(ux));
dy = median(diff(uy));

% Convert spacing to meters
switch posUnits
    case 'pixel'
        dx_m = dx * (pixelSizeNm*1e-9);
        dy_m = dy * (pixelSizeNm*1e-9);
    case 'um'
        dx_m = dx * 1e-6;
        dy_m = dy * 1e-6;
    case 'm'
        dx_m = dx;
        dy_m = dy;
    otherwise
        error('posUnits must be ''pixel'', ''um'', or ''m''.');
end

% Reshape into [ny x nx] assuming x varies fastest
Tx = reshape(tx, [nx, ny])';  % Pa
Ty = reshape(ty, [nx, ny])';  % Pa

% Solve MSM -> resultants (N/m)
[sxx, syy, sxy] = msmLeastNormFFT_resultant(Tx, Ty, dx_m, dy_m, lam);

% ---- von Mises (plane stress) ----
svm = sqrt( sxx.^2 - sxx.*syy + syy.^2 + 3*(sxy.^2) );   % units: N/m

% --- Save plots (use µm axes for readability) ---
switch posUnits
    case 'pixel'
        Xplot = (ux-ux(1)) * (pixelSizeNm/1000); % µm relative
        Yplot = (uy-uy(1)) * (pixelSizeNm/1000); % µm relative
        xlab = 'x (\mum)'; ylab = 'y (\mum)';
    case 'um'
        Xplot = ux-ux(1); Yplot = uy-uy(1);
        xlab = 'x (\mum)'; ylab = 'y (\mum)';
    case 'm'
        Xplot = (ux-ux(1))*1e6; Yplot = (uy-uy(1))*1e6;
        xlab = 'x (\mum)'; ylab = 'y (\mum)';
end

saveMap(Xplot, Yplot, sxx, '\sigma_{xx} (N/m)', fullfile(outDir, savePrefix + "_sigma_xx.png"), xlab, ylab);
saveMap(Xplot, Yplot, syy, '\sigma_{yy} (N/m)', fullfile(outDir, savePrefix + "_sigma_yy.png"), xlab, ylab);
saveMap(Xplot, Yplot, sxy, '\sigma_{xy} (N/m)', fullfile(outDir, savePrefix + "_sigma_xy.png"), xlab, ylab);
saveMap(Xplot, Yplot, svm, '\sigma_{vM} (N/m)', ...
    fullfile(outDir, savePrefix + "_sigma_vM.png"), xlab, ylab);

meta = struct();
meta.frameIdx = frameIdx;
meta.dx_m = dx_m; meta.dy_m = dy_m;
meta.posUnits = posUnits;
meta.pixelSizeNm = pixelSizeNm;
meta.lam = lam;
end

function [posXY, vecXY] = standardizePosVec(pos, vec)
% Returns Nx2 for both pos and vec

pos = double(pos);
vec = double(vec);

if size(pos,1)==2 && size(pos,2)~=2
    posXY = pos.';     % 2xN -> Nx2
elseif size(pos,2)==2
    posXY = pos;       % Nx2
else
    error('pos must be 2xN or Nx2.');
end

if size(vec,1)==2 && size(vec,2)~=2
    vecXY = vec.';     % 2xN -> Nx2
elseif size(vec,2)==2
    vecXY = vec;       % Nx2
else
    error('vec must be 2xN or Nx2.');
end
end

function saveMap(xv, yv, Z, ttl, outPng, xlab, ylab)
figure('Color','w');
imagesc(xv, yv, Z); axis image; set(gca,'YDir','normal');
colorbar; title(ttl);
xlabel(xlab); ylabel(ylab);
set(gca,'FontSize',12);
exportgraphics(gcf, outPng, 'Resolution', 300);
close(gcf);
end

function [sxx, syy, sxy] = msmLeastNormFFT_resultant(Tx, Ty, dx_m, dy_m, lam)
% Least-norm MSM in Fourier space returning stress resultants (N/m)
% Governing:
%   i*kx*Sxx + i*ky*Sxy = -Tx
%   i*kx*Sxy + i*ky*Syy = -Ty
%
% With Tx,Ty in Pa and k in rad/m -> S has units Pa*m = N/m.

[ny, nx] = size(Tx);
Txk = fft2(Tx);
Tyk = fft2(Ty);

kx = 2*pi*ifftshift( (-floor(nx/2):ceil(nx/2)-1) / (nx*dx_m) ); % rad/m
ky = 2*pi*ifftshift( (-floor(ny/2):ceil(ny/2)-1) / (ny*dy_m) ); % rad/m
[KX, KY] = meshgrid(kx, ky);

ikx = 1i*KX;
iky = 1i*KY;

% A*A' components
a11 = (KX.^2 + KY.^2);
a22 = (KX.^2 + KY.^2);
a12 = (KX.*KY);
a21 = a12;

% regularize
detM = (a11+lam).*(a22+lam) - a12.*a21;
inv11 = (a22+lam)./detM;
inv22 = (a11+lam)./detM;
inv12 = (-a12)./detM;
inv21 = (-a21)./detM;

b1 = -Txk;
b2 = -Tyk;

y1 = inv11.*b1 + inv12.*b2;
y2 = inv21.*b1 + inv22.*b2;

cikx = conj(ikx);
ciky = conj(iky);

Sxxk = cikx.*y1;
Sxyk = ciky.*y1 + cikx.*y2;
Syyk = ciky.*y2;

% remove k=0 (mean stress is undetermined)
Sxxk(1,1)=0; Sxyk(1,1)=0; Syyk(1,1)=0;

sxx = real(ifft2(Sxxk));
sxy = real(ifft2(Sxyk));
syy = real(ifft2(Syyk));
end