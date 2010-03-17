function testFM

I = imread('/home/sb234/Desktop/testFM.tif');

R = steerableFiltering(double(I), 'UnserM4', 1.5);

X1 = [34,123;30,98;20,93;12,87];
X2 = [56,55;44,34;24,11;50,43];

F = R - min(R(:)) + 1;

D = zeros(size(X1,1), size(X2,1));
P = cell(size(X1,1),size(X2,1));

for i = 1:size(X1,1)
    x = X1(i,:);
    
    [U E S R] = fastMarching(x, F);
    
    for j = 1:size(X2, 1)
        p = backPropagate(X2(j,:), X1(i,:), U);

        l = sum(sqrt((p(1:end-1,1) - p(2:end,1)).^2 + ...
            (p(1:end-1,2) - p(2:end,2)).^2));
        
        P{i,j} = p;
        
        D(i,j) = l;%U(X2(j,1),X2(j,2));% / l; % mean speed        
    end
end

[idx1,idx2] = lap(D);

figure, imshow(I,[]); hold on;

for i = 1:numel(idx1)
    p = P{idx1(i),idx2(i)};
    plot(p(:, 1), p(:, 2), 'Color', rand(3, 1), 'LineWidth', 2);
end
hold off;
