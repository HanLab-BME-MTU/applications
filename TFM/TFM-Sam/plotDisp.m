
% Code plots the displacement vectors calculated by TPT.mat
% exampleRunFile.m creates the 'x' and 'track' vectors used in this script

%\\ Calculate Displacement Vectors
relP=zeros(length(x{1}{1}(:,1)),3);
for i = 1:length(x{1}{1})
    for k = 1:3
        if track{1}{1}(i) ~= 0
          relP(i,k) = x{2}{1}(track{1}{1}(i),k)-x{1}{1}(i,k);
        elseif track{1}{1}(i) == 0
          relP(i,k) = 0;
        else
        
        end
    end
end

%Extract coordinates of reference image points
x1=x{1}{1}(:,1);
y1=x{1}{1}(:,2);
z1=x{1}{1}(:,3);

%Calculate magnitude of each vector
mags=sqrt(relP(:,1).^2+relP(:,2).^2+relP(:,3).^3);
mags=real(mags);

%Enables colormapping section, 0 = no color mapping
mapcolor=1;
%Displays the full axes, starting from 0. 0 = cut axes
fullaxes=1;

%\\ Plotting
hold on
%Set plot of every nth point
n = 1;
figure(1)
q=quiver3(x1(n:n:end),y1(n:n:end),z1(n:n:end),relP(n:n:end,1),relP(n:n:end,2),relP(n:n:end,3));
if fullaxes == 1
    xlim([0 256])
    ylim([0 256])
    zlim([0 20])
end

%\\ Colormapping
if mapcolor == 1
%Get the current colormap
currentColormap = colormap(gca);
%Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
get(q.Head)
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
get(q.Head)

%We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
get(q.Tail)
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
get(q.Tail)

set(gca,'Color','k')

title('T-PT Tracking Displacement Field - 0.5kPa')
xlabel('Y COORD')
ylabel('X COORD')
zlabel('Z COORD')
hold off
end