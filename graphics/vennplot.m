% Francois Aguet, 04/28/2013

function vennplot(a, b, ab, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('a');
ip.addRequired('b');
ip.addRequired('ab');
ip.addOptional('labels', {'A', 'B'});
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Hues', [0.33 0.55]);
ip.addParamValue('Handle', gca, @ishandle);
ip.addParamValue('DisplayMode', '');
ip.addParamValue('Font', []);
ip.parse(a, b, ab, varargin{:});
hues = ip.Results.Hues;

[r1, r2, d] = getVennParameters(a, b, ab);
tfont = ip.Results.Font;
if isempty(tfont)
    tfont = {'FontName', get(0, 'DefaultAxesFontName'),...
        'FontSize', get(0, 'DefaultAxesFontSize')};
end

aw = sqrt(1/(2*pi))*2;
ah = sqrt(1/pi);
axis equal;
axis([d/2-aw d/2+aw -ah ah]);
text(d/2-aw/2, ah, ip.Results.labels{1}, tfont{:}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(d/2+aw/2, ah, ip.Results.labels{2}, tfont{:}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

dt = pi/20;
theta = (0:dt:2*pi)';

xi = (r1^2-r2^2+d^2)/(2*d);
yi = real(sqrt((4*d^2*r1^2 - (d^2-r2^2+r1^2)^2)/(4*d^2)));


tl1 = atan2(yi, xi);
th1 = atan2(-yi, xi)+2*pi;
tl2 = atan2(yi, -d+xi);
th2 = atan2(-yi, -d+xi)+2*pi;

theta = sort([theta; th1; th2; tl1; tl2]);

% coordinates of the two circles
x1 = r1*cos(theta);
y1 = r1*sin(theta);
x2 = r2*cos(theta)+d;
y2 = r2*sin(theta);

t1 = theta(theta>=tl1 & theta<=th1);
t1i = [theta(theta>=th1); theta(theta<=tl1)];
t2i = theta(theta>=tl2 & theta<=th2);
t2 = [theta(theta>=th2); theta(theta<=tl2)];

plot(x1, y1);
plot(x2, y2);


c1 = hsv2rgb([hues(1) 0.8 1]);
c2 = hsv2rgb([hues(2) 0.8 1]);
c12 = (c1+c2)/2;

fill([r1*cos(t1); r2*cos(t2i(end:-1:1))+d], [r1*sin(t1); r2*sin(t2i(end:-1:1))], c1, 'EdgeColor', 'none');
fill([r1*cos(t1i); r2*cos(t2i)+d], [r1*sin(t1i); r2*sin(t2i)], c12, 'EdgeColor', 'none');
fill([r2*cos(t2)+d; r1*cos(t1i)], [r2*sin(t2); r1*sin(t1i)], c2, 'EdgeColor', 'none');



function [r1, r2, d] = getVennParameters(a, b, ab)

r1 = sqrt((a+ab)/pi);
r2 = sqrt((b+ab)/pi);


opts = optimset('Jacobian', 'off', ...
    'MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-6, ...
    'Tolfun', 1e-6);


d = fminbnd(@(x) areaCost(x, r1, r2, ab), 0, r1+r2, opts);

function a = iArea(r1, r2, d)
a = r2^2*acos((d.^2+r2^2-r1^2)./(2*d*r2)) + r1^2*acos((d.^2+r1^2-r2^2)./(2*d*r1)) -...
    0.5*sqrt((-d+r2+r1).*(d+r2-r1).*(d-r2+r1).*(d+r2+r1));

function v = areaCost(d, r1, r2, a)
v = abs(iArea(r1, r2, d) - a)/a;

