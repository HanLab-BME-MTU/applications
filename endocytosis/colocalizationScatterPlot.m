% Francois Aguet, June 2010

function R = colocalizationScatterPlot(x, y, av, name1, name2)

if nargin<3
    av = [];
end
if nargin<4
    name1 = 'channel 1';
end
if nargin<5
    name2 = 'channel 2';
end

mux = mean(x);
muy = mean(y);
N = length(x);

sxx = sum(x.^2) - N*mux^2; % N*var(x,1)
syy = sum(y.^2) - N*muy^2;
sxy = sum(x.*y) - N*mux*muy; % N*cov(x,y)

R = sxy^2/(sxx*syy);


m = sxy/sxx;
b = 1/N*sum(y-m*x);

u = [1; m]/sqrt(1+m^2);
v = [x; y-b]; % position vectors
normv = sqrt(sum(v.^2,1));

% directional angle:
alpha = atan2(u'*[-v(2,:); v(1,:)], u'*v); % tan -> no normalization needed

theta = acos(u(1));
dx = x*sin(theta) - (y-b)*cos(theta);

% % variance of the distance
% sin(theta)^2*sxx/N + cos(theta)^2*syy/N - sin(2*theta)*sxy/N
% (sxx^2*syy - sxy^2*sxx)/(sxx^2+sxy^2) / N
% (1-R) / (1/syy+R/sxx) / N
% (1-R) / (1/var(y,1) + R/var(x,1))
var_d = (1-R) / (1/var(y,1) + R/var(x,1));


figure;
plot(x, y, 'ko');
axis equal;
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
xlabel(name1, 'FontName', 'Helvetica', 'FontSize', 16);
ylabel(name2, 'FontName', 'Helvetica', 'FontSize', 16);
axis(av);
%title(['R^2 = ' num2str(R) ', \sigma_d = ' num2str(sqrt(var_d))]);
title(['R^2 = ' num2str(R)]);