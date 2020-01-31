% Code plots the ground truth displacement field from dispField.mat file

%\\ Extracting Displacement Values of Ground Truth Field
data=load('displField.mat');
ux=data.ux(:,:,32);
uy=data.uy(:,:,32);
uz=data.uz(:,:,32);

%Meshgrid of Required Size
[X,Y,Z]=meshgrid(0:1:197,0:1:197,1);

%\\ Plotting Ground Truth
n=5; %plot every nth point
figure(2)
q=quiver3(X(n:n:end),Y(n:n:end),Z(n:n:end),ux(n:n:end),uy(n:n:end),uz(n:n:end));

%Calculate magnitude of each vector
%Future Use: Coloring vectors based on magnitude
mags=sqrt(ux.^2+uy.^2+uz.^3);
mags=real(mags);