%[crbV] = gaussianFitCRB returns the Cramer-Rao bounds for the parameters of a 2D Gaussian fit

function crbV = gaussianFitCRB(X, x0, A, c, sigma, sigma_n, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('NoiseDistr', 'Gaussian', @(x) any(strcmpi(x, {'Gaussian', 'Poisson'})));
% ip.addParamValue('Dims', 2);
ip.parse(varargin{:});

switch size(X,2)
    case 2
        g = exp(-((X(:,1)-x0(1)).^2+(X(:,2)-x0(1)).^2)/(2*sigma^2));
        dA = g;
        dc = ones(size(g));
        dx0 = (X(:,1)-x0(1))/sigma^2*A.*g;
        dy0 = (X(:,2)-x0(2))/sigma^2*A.*g;
        g = A*g + c;
        N = numel(g);
    
        % 4 parameters: x, y, A, c
        crbV = zeros(4,numel(sigma_n));
        if strcmpi(ip.Results.NoiseDistr, 'Gaussian')
            crbV(1,:) = sigma_n*sqrt(2/pi)/A;
            crbV(2,:) = sigma_n*sqrt(2/pi)/A;
            %D = [sum(dA(:).^2) sum(dA(:).*dc(:));
            %    sum(dA(:).*dc(:)) sum(dc(:).^2)];
            D = [pi*sigma^2 2*pi*sigma^2;
                 2*pi*sigma^2 N];
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
    case 3
        if numel(sigma)==1
            sigma = [sigma sigma];
        end
        g = exp(-((X(:,1)-x0(1)).^2+(X(:,2)-x0(1)).^2)/(2*sigma(1)^2)) ...
            .* exp(-((X(:,3)-x0(3)).^2)/(2*sigma(2)^2));
        dA = g;
        dc = ones(size(g));
        dx0 = (X(:,1)-x0(1))/sigma(1)^2*A.*g;
        dy0 = (X(:,2)-x0(2))/sigma(1)^2*A.*g;
        dz0 = (X(:,3)-x0(3))/sigma(2)^2*A.*g;
        g = A*g + c;
        N = numel(g);
        
        % 5 parameters: x, y, z, A, c
        crbV = zeros(5,numel(sigma_n));
        if strcmpi(ip.Results.NoiseDistr, 'Gaussian')
            
            crbV(1,:) = sigma_n*sqrt(2/(A^2*pi^1.5*sigma(2)));
            crbV(2,:) = sigma_n*sqrt(2/(A^2*pi^1.5*sigma(2)));
            crbV(3,:) = sigma_n*sqrt(2*sigma(2)/(A^2*pi^1.5*sigma(1)^2));
            D = [pi^1.5*sigma(1)^2*sigma(2) 2*sqrt(2)*pi^1.5*sigma(1)^2*sigma(2);
                2*sqrt(2)*pi^1.5*sigma(1)^2*sigma(2) N];
            
            %D = [sum(dA(:).^2) sum(dA(:).*dc(:));
            %    sum(dA(:).*dc(:)) sum(dc(:).^2)];
            D = sqrt(diag(inv(D)));
            crbV(4,:) = sigma_n*D(1);
            crbV(5,:) = sigma_n*D(2);
            % F = [sum(dx0(:).^2) 0 0 0;
            %      0 sum(dy0(:).^2) 0 0;
            %      0 0 sum(dA(:).^2) sum(dA(:).*dc(:));
            %      0 0 sum(dA(:).*dc(:)) sum(dc(:).^2)];
            % crbV = sqrt(diag(inv(F)));
        else
            
        end
end
        


