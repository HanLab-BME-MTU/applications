function [frame, xv, yv, sv, Av] = simGaussianBeads(nx, ny, sigma, varargin)
% SIMGAUSSIANBEADS generates a given number of beads using 2D Gaussians in an image.
% The generated Gaussian signals do not overlap with the image boundaries.
% The difference from SIMGAUSSIANSPOTS are
% 1) intensities do not superimpose when overlapped
% 2) beads are created with variable amplitudes
% 3) there are shadows around the new bead (intensity diminishes by 50% if
% there is any positive pixel values
%
%   Usage: [frame,xv,yv,sv,Av]=simGaussianSpots(nx,ny,sigma,varargin)
%   Input:
%           nx:    image size in x direction
%           ny:    image size in y direction
%           sigma: 2D Gaussian standard deviation (in pixels)
%   Options:
%           'x': x coordinates of centers of 2D Gaussians (can be subpixel)
%           'y': y coordinates of centers of 2D Gaussians (can be subpixel)
%           'A': amplitudes of 2D Gaussians
%           'npoints'   : number of 2D Gaussian to be generated
%           'background': value of background
%           'Border' : border conditions: 'padded' (default), 'periodic', or 'truncated'
%           'Normalization: {'on' | 'off' (default)} divides Gaussians by 2*pi*sigma^2 when 'on'
%           'Verbose': {'on' | 'off' (default)}
%   Output:
%           frame: image with 2D Gaussian signals
%           xv:    x coordinates of centers of Gaussians (can be subpixel)
%           vy:    y coordinates of centers of Gaussians (can be subpixel)
%           sv:    vector of standard deviations
%           Av:    vector of amplitudes
%
% Example:
% img = simGaussianSpots(200, 100, 2, 'npoints', 50, 'Border', 'periodic');

% Francois Aguet, last modified July 30, 2012

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('nx',@isnumeric);
ip.addRequired('ny',@isnumeric);
ip.addRequired('sigma',@isnumeric);
ip.addParamValue('x', []);
ip.addParamValue('y', []);
ip.addParamValue('A', []);
ip.addParamValue('npoints', 1);
ip.addParamValue('background', 0);
ip.addParamValue('verbose', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('Border', 'padded', @(x) any(strcmpi(x, {'padded', 'periodic', 'truncated'})));
ip.addParamValue('Normalization', 'off', @(x) any(strcmpi(x, {'analytical', 'sum', 'off'})));
ip.parse(nx, ny, sigma, varargin{:});

np = ip.Results.npoints;
c = ip.Results.background;
xv = ip.Results.x(:);
yv = ip.Results.y(:);
Av = ip.Results.A(:);
sv = ip.Results.sigma(:);

if numel(xv) ~= numel(yv)
    error('''x'' and ''y'' must have the same size.');
end

if ~isempty(xv)
    np = length(xv);
end

% feature vectors
if numel(sv)==1
    sv = sv*ones(np,1);
end
wv = ceil(4*sv);

w_max = 2*max(wv)+1;
if (w_max>nx || w_max>ny)
    error(['For a sigma of ' num2str(max(sv), '%.1f') ', nx and ny must be greater than ' num2str(2*w_max+1)]);
end

if numel(Av)==1
    Av = Av*ones(np,1);
end
if isempty(Av)
    Av = ones(np,1);
end

% Generate point coordinates, if input is empty
if strcmpi(ip.Results.Border, 'padded');
    % discard signals close to image border
    if ~isempty(xv)
        idx = xv <= wv+1 | yv <= wv+1 | xv > nx-wv-1 | yv > ny-wv-1;
        xv(idx) = [];
        yv(idx) = [];
        Av(idx) = [];
        sv(idx) = [];
        if strcmpi(ip.Results.verbose, 'on')
            fprintf('Number of discarded points: %d\n', numel(idx));
        end
        np = length(xv);
    else
        xv = (nx-2*wv-1).*rand(np,1) + wv+1;
        yv = (ny-2*wv-1).*rand(np,1) + wv+1;
        [~, idx] = sort(xv+yv*nx); % sort spots according to row position
        xv = xv(idx);
        yv = yv(idx);
    end
else
%     if ~isempty(xv) && (any(min([xv; yv])<0.5) || any(max(xv)>=nx+0.5) || any(max(yv)>=ny+0.5))
%         error('All points must lie within x:[0.5 nx+0.5), y:[0.5 ny+0.5).');
%     end
    if isempty(xv)
        xv = nx*rand(np,1)+0.5;
        yv = ny*rand(np,1)+0.5;
    end
end

% background image
frame = c*ones(ny, nx);

xi = round(xv);
yi = round(yv);
dx = xv-xi;
dy = yv-yi;

if strcmpi(ip.Results.Normalization, 'analytical')
    Av = Av ./ (2*pi*sv.^2);
end

switch ip.Results.Border
    case 'padded'
        for k = 1:np
            wi = wv(k);
            xa = xi(k)-wi:xi(k)+wi;
            ya = yi(k)-wi:yi(k)+wi;
            [xg,yg] = meshgrid(-wi:wi);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g =  Av(k)*g;
            end
            frame(ya,xa) = frame(ya,xa) + g;
        end
    case 'periodic'
        lbx = xi-wv;
        ubx = xi+wv;
        lby = yi-wv;
        uby = yi+wv;
        
        for k = 1:np
            shifts = [0 0];
            if lbx(k)<1
                shifts(2) = 1-lbx(k);
            elseif ubx(k)>nx
                shifts(2) = nx-ubx(k);
            end
            if lby(k)<1
                shifts(1) = 1-lby(k);
            elseif uby(k)>ny
                shifts(1) = ny-uby(k);
            end
            wi = -wv(k):wv(k);
            [xg,yg] = meshgrid(wi,wi);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g =  Av(k)*g;
            end
            xa = (xi(k)-wv(k):xi(k)+wv(k)) + shifts(2);
            ya = (yi(k)-wv(k):yi(k)+wv(k)) + shifts(1);
            if all(shifts==0)
                frame(ya,xa) = frame(ya,xa) + g;
            else
                frame = circshift(frame, shifts);
                frame(ya,xa) = frame(ya,xa) + g;
                frame = circshift(frame, -shifts);
            end
        end
    case 'truncated'
        lbx = max(xi-wv,1);
        ubx = min(xi+wv,nx);
        lby = max(yi-wv,1);
        uby = min(yi+wv,ny);
        
        for k = 1:np
            wx = (lbx(k):ubx(k)) - xi(k);
            wy = (lby(k):uby(k)) - yi(k);
            [xg,yg] = meshgrid(wx,wy);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g =  Av(k)*g;
            end
            xa = lbx(k):ubx(k);
            ya = lby(k):uby(k);
            frame(ya,xa) = frame(ya,xa) + g;
        end
end

% old code - this is wrong
% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('nx',@isnumeric);
% ip.addRequired('ny',@isnumeric);
% ip.addRequired('sigma',@isnumeric);
% ip.addParamValue('x', []);
% ip.addParamValue('y', []);
% ip.addParamValue('A', []);
% ip.addParamValue('npoints', 1);
% ip.addParamValue('background', 0);
% ip.addParamValue('verbose', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
% ip.addParamValue('Border', 'padded', @(x) any(strcmpi(x, {'padded', 'periodic', 'truncated'})));
% ip.addParamValue('Normalization', 'off', @(x) any(strcmpi(x, {'analytical', 'sum', 'off'})));
% ip.parse(nx, ny, sigma, varargin{:});
% 
% np = ip.Results.npoints;
% c = ip.Results.background;
% xv = ip.Results.x(:);
% yv = ip.Results.y(:);
% Av = ip.Results.A(:);
% sv = ip.Results.sigma(:);
% 
% if numel(xv) ~= numel(yv)
%     error('''x'' and ''y'' must have the same size.');
% end
% 
% if ~isempty(xv)
%     np = length(xv);
% end
% 
% % feature vectors
% if numel(sv)==1
%     sv = sv*ones(np,1);
% end
% wv = ceil(4*sv);
% 
% w_max = 2*max(wv)+1;
% if (w_max>nx || w_max>ny)
%     error(['For a sigma of ' num2str(max(sv), '%.1f') ', nx and ny must be greater than ' num2str(2*w_max+1)]);
% end
% 
% if numel(Av)==1
%     Av = Av*ones(np,1);
% end
% if isempty(Av)
% %     Av = ones(np,1);
%     Av = 0.4+0.6*rand(np,1);
% end
% 
% % Generate point coordinates, if input is empty
% if isempty(xv)
%     xv = nx*rand(np,1)+0.5;
%     yv = ny*rand(np,1)+0.5;
% end
% 
% % Here beads are created in magnified frame then will be shrunk to
% % the real scale
% mf = 10;
% % background image
% frame = c*ones(ny, nx);
% framemag = c*ones(ny*mf, nx*mf);
% 
% xi = round(xv*mf);
% yi = round(yv*mf);
% dx = xv*mf-xi;
% dy = yv*mf-yi;
% wvm = ceil(4*sv*mf);
% 
% switch ip.Results.Border
%     case 'truncated'
%         
%         lbx = max(xi-wvm,1);
%         ubx = min(xi+wvm,nx*mf);
%         lby = max(yi-wvm,1);
%         uby = min(yi+wvm,ny*mf);
%         
%         for k = 1:np
%             wx = (lbx(k):ubx(k)) - xi(k);
%             wy = (lby(k):uby(k)) - yi(k);
%             [xg,yg] = meshgrid(wx,wy);
%             g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*(mf*sv(k))^2));
%             g =  Av(k)*g;
%             xa = lbx(k):ubx(k);
%             ya = lby(k):uby(k);
%             % Here, we don't just superimpose intensities but overwrite
%             % them on the location of the bead
%             ge = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*(mf*sv(k))^2));
%             beadMask = (ge>0.19);
%             beadNeigh = ~beadMask;
%             beadBoundary = bwperim(beadMask);
%             beadBoundary2 = (bwdist(beadNeigh)>0 &  bwdist(beadNeigh)<=4).*((5-bwdist(beadBoundary))/5);
%             framemag(ya,xa) = framemag(ya,xa).*beadNeigh + framemag(ya,xa).*beadBoundary2 + g;% - 0.5*frame(ya,xa).*frame(ya,xa).*(g<0.0001).*(g>0);
% %             framemag(ya,xa) = framemag(ya,xa).*beadNeigh + g;
%         end
%         frame = imresize(framemag,1/mf);
% end
