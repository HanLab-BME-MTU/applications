function  make_adhesion_data(bild_datei, adh_mask_datei, TFM_results,frame, grid_mat, fnorm, xshift,yshift)

if nargin < 8
    xshift = 0;
    yshift = 0;
end

data = load(adh_mask_datei);
C = contourc(double(data.adh_mask),1);
bild = imread(bild_datei);

n = 1;
i = 1;
sites = [];
while i < size(C,2)
    sites(n).contour(1:C(2,i),1:2) = C(1:2,i+(1:C(2,i)))';
    sites(n).contour(1:C(2,i),1) = sites(n).contour(1:C(2,i),1)+ xshift;
    sites(n).contour(1:C(2,i),2) = sites(n).contour(1:C(2,i),2)+ yshift;
    i = i+C(2,i)+1;
    n = n +1;
end

for i=1:size(sites,2)
    inpos = inpolygon(TFM_results(frame).pos(:,1)+ TFM_results(frame).vec(:,1),TFM_results(frame).pos(:,2) + TFM_results(frame).vec(:,2),sites(i).contour(:,1),sites(i).contour(:,2));
    if nnz(inpos) > 0
        traction_data(i).max = max((TFM_results(frame).force(inpos,1).^2 + TFM_results(frame).force(inpos,2).^2).^0.5);
        traction_data(i).vec(1:2) = [mean(TFM_results(frame).force(inpos,1)),mean(TFM_results(frame).force(inpos,2))];
        traction_data(i).abs = norm(traction_data(i).vec(1:2));
        traction_data(i).pos(1:2) =  [mean(sites(i).contour(:,1)),mean(sites(i).contour(:,2))];
    else
        traction_data(i).pos = nan;
        traction_data(i).max = nan;
        traction_data(i).vec(1:2) = nan;
        traction_data(i).max = nan;
    end
    
    in_fa = inpolygon(data.vert_pos,data.hor_pos,sites(i).contour(:,1),sites(i).contour(:,2));
    sites(i).mean_brightness = mean(mean(bild(in_fa)));
    sites(i).area = polyarea(sites(i).contour(:,1),sites(i).contour(:,2));
end

figure; imagesc(bild); hold; colormap gray;
contour(data.vert_pos, data.hor_pos,double(data.adh_mask),1,'y');
for n=1:size(sites,2)
        text(mean(sites(n).contour(:,1)),mean(sites(n).contour(:,2)), num2str(n),'Color','g');
end
hold off;

figure;
plot([sites(~isnan([traction_data.max])).mean_brightness],[traction_data(~isnan([traction_data.max])).max],'r*');
figure;
plot([sites(~isnan([traction_data.max])).area],[traction_data(~isnan([traction_data.max])).max],'*');

figure, hold on; colormap jet;
surf(grid_mat(:,:,1), grid_mat(:,:,2),fnorm),view(0,90),shading interp, axis equal;
set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse'),
maxfnorm = max(max(fnorm));
for n=1:size(sites,2)
        plot3(maxfnorm*ones(size(sites(i).contour(:,1))),sites(i).contour(:,1),sites(i).contour(:,2),'-k')
end
hold off;

save([bild_datei(1:end-4),'_adhesion_data.mat'], 'sites','traction_data','-mat');