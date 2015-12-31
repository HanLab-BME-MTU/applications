%% Simulate an EDF Line-Illumination


clc
clear

MATRIX_SIZE=600;            %Size of Matrix
IMAGING_MODE=1;             %1 for 1-photon, 2 for 2-photon.


%%%%%%%%%%%%%%%%%%%%%%%% Illumination descriptor.

if IMAGING_MODE == 2
    EXCITATION_WAVELENGTH=0.900; %In microns.
else
    EXCITATION_WAVELENGTH=0.488;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Objective descriptors

PUPIL_PLANE_MAX_RADIUS=4.7945; %Back aperture radius in millimeters

FOCAL_LENGTH_OBJECTIVE=7; %Special Optics, focal length in millimeters.

REFRACTIVE_INDEX_WATER=1.3333;

THETA_MAX=asind(PUPIL_PLANE_MAX_RADIUS./(FOCAL_LENGTH_OBJECTIVE.*REFRACTIVE_INDEX_WATER));


%%%%%%%%%%%%%%%%%%%%%%%%%%% Spherical shell descriptors

SPHERICAL_SHELL_RADIUS=MATRIX_SIZE/6;     %Spherical shell radius, which for the excitation OTF would be 1/940, but for the emission OTF 1/540?

DELTA_K=REFRACTIVE_INDEX_WATER/(EXCITATION_WAVELENGTH*SPHERICAL_SHELL_RADIUS);     %1.33 is the refractive index of water.

DELTA_X=1/(MATRIX_SIZE*DELTA_K); %In nanometers, 176.69 nm

IMAGE_FIELD_OF_VIEW=ones(MATRIX_SIZE,1);

IMAGE_FIELD_OF_VIEW=[DELTA_X:DELTA_X:DELTA_X.*MATRIX_SIZE];

BesselFourierSphericalShell=zeros(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);
%Preallocate the matrices to increase speed.


for k=1:MATRIX_SIZE
    
    for l=1:MATRIX_SIZE
        
        for j=1:MATRIX_SIZE
            
            x=k-MATRIX_SIZE/2;
            
            y=l-MATRIX_SIZE/2;
            
            z=j-MATRIX_SIZE/2;
            
            q=sqrt(x^2+y^2+z^2);
            
            if q<SPHERICAL_SHELL_RADIUS && q > SPHERICAL_SHELL_RADIUS-3
                
                r2=sqrt(x^2+y^2);
                
                alpha=atan2d(r2,z);
                
                if alpha < THETA_MAX ;
                    
                    BesselFourierSphericalShell(k,l,j)=1;
                    
                end
            end
        end
    end
end

clear k l j r2 x y z

%%
lineMatrix=zeros(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);
lineWidth=2;
start=MATRIX_SIZE/2-lineWidth/2;
stop=MATRIX_SIZE/2+lineWidth/2;

lineMatrix(start:stop,:,:)=1;

BesselFourierSphericalShell2=lineMatrix.*BesselFourierSphericalShell;

clear start stop lineWidth
%Take the spherical shell that you created previously, perform a
%multidimensional FFT on it, and then shift the Fourier transform to the DC
%component.  Specify your imaging type here; 1=1-photon, 2=2-photon.

%%

BesselEField=fftshift(fftn(BesselFourierSphericalShell2));

BesselIllumination=abs(BesselEField).^(2*IMAGING_MODE);

%%

BesselCrossSection1=BesselIllumination(:,round(MATRIX_SIZE/2),:);

BesselCrossSection1=squeeze(BesselCrossSection1);

figure; 

imshow(BesselCrossSection1,[]); colormap jet;


BesselCrossSection2=BesselIllumination(round(MATRIX_SIZE/2),:,:);

BesselCrossSection2=squeeze(BesselCrossSection2);

figure;

imshow(BesselCrossSection2,[]); colormap jet;


BesselCrossSection3=BesselIllumination(:,:,round(MATRIX_SIZE/2));

BesselCrossSection3=squeeze(BesselCrossSection3);

figure;

imshow(BesselCrossSection3,[]); colormap jet;

%%
MAXIMUM_VALUE=max(max(max(BesselIllumination)));

BesselIlluminationNormalized=BesselIllumination./MAXIMUM_VALUE;

[faces,verts,colors]=isosurface(BesselIlluminationNormalized,0.8,BesselIlluminationNormalized);

p1 = patch('Vertices',verts,'Faces', faces, ...
    'FaceVertexCData', colors, ...
    'FaceColor', 'interp', ...
    'edgecolor', 'interp');

isonormals(BesselIlluminationNormalized, p1);

view(30,-15);
axis vis3d;
lighting gouraud
camlight left

%%
[x,y,z,v] = flow;

xmin = min(x(:)); 
ymin = min(y(:)); 
zmin = min(z(:));

xmax = max(x(:)); 
ymax = max(y(:)); 
zmax = max(z(:));

hslice=surf(linspace(0,1,100),linspace(0,1,100),zeros(100));

rotate(hslice,[-1,0,0],-60);
xd = get(hslice,'XData');
yd = get(hslice,'YData');
zd = get(hslice,'ZData');

delete hslice

figure
colormap(jet)
h = slice(x,y,z,v,xd,yd,zd);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

hold on
hx = slice(x,y,z,v,xmax,[],[]);
hx.FaceColor = 'interp';
hx.EdgeColor = 'none';

hy = slice(x,y,z,v,[],ymax,[]);
hy.FaceColor = 'interp';
hy.EdgeColor = 'none';

hz = slice(x,y,z,v,[],[],zmin);
hz.FaceColor = 'interp';
hz.EdgeColor = 'none';

daspect([1,1,1])
axis tight
view(-38.5,16)
camzoom(1.4)
camproj perspective
lightangle(-45,45)
colormap (jet(24))

%%

