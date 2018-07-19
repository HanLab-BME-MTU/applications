function [focals,forces] = point_forces(bild_datei, pos, vec,cutoff)

bild = imread(bild_datei);
figure;
imagesc(bild), axis equal, colormap gray; hold;

focals = [];
while 1
    [x,y] = ginput(1);
     if isempty(x)
            break;
     else
        scatter(x,y,'*y');
        focals(end+1,1:2) = [x,y];
     end
end  
hold off;

for i = 1:size(pos,1)
   for j = 1:size(focals,1)
      G(2*i-1:2*i,2*j-1:2*j) = greens(pos(i,1:2)-focals(j,1:2),cutoff);
   end
end
disp('SVD');
[U,s,V] = csvd(G);

u(1:2:size(vec,1)*2,1) = vec(:,1);
u(2:2:size(vec,1)*2,1) = vec(:,2);

close(gcf);
figure;
subplot(2,1,1);
    disp(['Vorschlag f�r L-kurven Knick: Lambda=', num2str(l_curve(U,s,u,'Tikh'))]);
subplot(2,1,2);
    disp(['Vorschlag f�r das Gcv Minimum: Lambda=', num2str(gcv(U,s,u,'Tikh'))]);


while 1
    lambda = input('Regularization parameter ? (0 to exit) ');
    
    if lambda == 0
        break
    end
    
    [F,rho,eta]   = tikhonov(U,s,V,u,lambda);
    
    % forces as matrix
    forces(1:length(F)/2,1) = F(1:2:end); 
    forces(1:length(F)/2,2) = F(2:2:end); 
    
    figure;
    imagesc(bild), axis equal, hold;
    quiver(pos(:,1), pos(:,2), vec(:,1), vec(:,2),'g');
    quiver(focals(:,1), focals(:,2), forces(:,1), forces(:,2),'r');
end

function G = greens(r,cutoff)

    x = r(1);
    y = r(2);
    mag = sqrt(x^2+y^2);

    if (mag <= cutoff)
       G(1:2,1:2) = 0;
    else

        G(1,1) = 3/(4*pi*mag)*(1+x*x/mag/mag);
        G(1,2) = 3*x*y/(4*pi*mag*mag*mag);
        G(2,1) = G(1,2);
        G(2,2) = 3/(4*pi*mag)*(1+y*y/mag/mag);
    end
end
end

