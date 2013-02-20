function crbV = gaussianFitCRB(x, y, x0, y0, A, c, sigma, sigma_n, noiseDistr)

if nargin<9
    noiseDistr = 'Gaussian';
end

g = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));
dA = g;
dc = ones(size(g));
dx0 = (x-x0)/sigma^2*A.*g; 
dy0 = (y-y0)/sigma^2*A.*g;
g = A*g + c;

crbV = zeros(4,numel(sigma_n));

if nargin<8 || strcmpi(noiseDistr, 'Gaussian')
    crbV(1,:) = sigma_n*sqrt(2/pi)/A;
    crbV(2,:) = sigma_n*sqrt(2/pi)/A;
    D = [sum(dA(:).^2) sum(dA(:).*dc(:));
        sum(dA(:).*dc(:)) sum(dc(:).^2)];
    D = sqrt(diag(inv(D)));
    crbV(3,:) = sigma_n*D(1);
    crbV(4,:) = sigma_n*D(2);
    % F = [sum(dx0(:).^2) 0 0 0;
    %      0 sum(dy0(:).^2) 0 0;
    %      0 0 sum(dA(:).^2) sum(dA(:).*dc(:));
    %      0 0 sum(dA(:).*dc(:)) sum(dc(:).^2)];
    % crbV = sqrt(diag(inv(F)));
else
    F = [sum(dx0(:).^2./g(:)) 0 0 0;
         0 sum(dy0(:).^2./g(:)) 0 0;
         0 0 sum(dA(:).^2./g(:)) sum(dA(:).*dc(:)./g(:));
         0 0 sum(dA(:).*dc(:)./g(:)) sum(dc(:).^2./g(:))];
    F = sqrt(diag(inv(F)));
    crbV(1,:) = F(1)./sqrt(sigma_n);
    crbV(2,:) = F(2)./sqrt(sigma_n);
    crbV(3,:) = F(3)./sqrt(sigma_n);
    crbV(4,:) = F(4)./sqrt(sigma_n);
end
