function R = focalAdhesionDetector(I, BW, sigma)

hside = ceil(3 * sigma);
x = -hside:hside;
g0 = (1/(sqrt(2*pi) * sigma)) * exp(-x.^2 / (2 * sigma^2));
g1 = - (x/(sqrt(2*pi) * sigma^3)) .* exp(-x.^2 / (2 * sigma^2));
g2 = exp(-x.^2 / (2 * sigma^2)) .* (x.^2 - sigma^2) / (sqrt(2*pi) * sigma^5);

Ixy = conv2(g1,g1,I, 'same');
Ixx = conv2(g0,g2,I, 'same');
Iyy = conv2(g2,g0,I, 'same');

L = bwlabel(BW);
labels = 1:max(L(:));

Sxy = arrayfun(@(l) sum(Ixy(L == l)), labels');
Sxx = arrayfun(@(l) sum(Ixx(L == l)), labels');
Syy = arrayfun(@(l) sum(Iyy(L == l)), labels');

theta = -pi/2:pi/64:pi/2-pi/64;
ct = cos(theta);
st = sin(theta);

R = Syy * ct.^2 - Sxy * (2 * ct .* st) + Sxx * st.^2;
[R T] = max(R,[],2);

imshow(I,[]); hold on;
for l = 1:numel(labels)
    ind = find(L == l);
    [y x] = ind2sub(size(L), ind);
    ct = ones(numel(ind), 1) * cos(T(l));
    st = ones(numel(ind), 1) * sin(T(l));
    quiver(x,y,ct,st);
end
hold off;
