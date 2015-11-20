%% Program calculates the excitation protfile of a Bessel beam.
%Specify the imaging parameters.
clc

clear

MATRIX_SIZE=600;            %Size of Matrix

EXCITATION_LAMBDA=488;      %Excitation wavelength [nm]

IMAGING_MODE=1;             %1 for 1-photon, 2 for 2-photon.

PUPIL_PLANE_MAX_RADIUS=4.7945; %Back aperture radius in millimeters

PUPIL_PLANE_MIN_RADIUS=PUPIL_PLANE_MAX_RADIUS-.05;

FOCAL_LENGTH_OBJECTIVE=7; %Special Optics, focal length in millimeters.

REFRACTIVE_INDEX_WATER=1.3333;

THETA_MAX=asind(PUPIL_PLANE_MAX_RADIUS./(FOCAL_LENGTH_OBJECTIVE.*REFRACTIVE_INDEX_WATER));

THETA_MIN=asind(PUPIL_PLANE_MIN_RADIUS./(FOCAL_LENGTH_OBJECTIVE.*REFRACTIVE_INDEX_WATER));



%%

SPHERICAL_SHELL_RADIUS=0.3.*MATRIX_SIZE;     %Spherical shell radius, which for the excitation OTF would be 1/940, but for the emission OTF 1/540?

DELTA_K=REFRACTIVE_INDEX_WATER/(EXCITATION_LAMBDA*SPHERICAL_SHELL_RADIUS);     %1.33 is the refractive index of water.

DELTA_X=1/(MATRIX_SIZE*DELTA_K); %In nanometers, 176.69 nm

IMAGE_FIELD_OF_VIEW=ones(MATRIX_SIZE,1);

IMAGE_FIELD_OF_VIEW=[DELTA_X:DELTA_X:DELTA_X.*MATRIX_SIZE];


%% Paint the Fourier Shell of the high-NA Gaussian beam.

gaussianFourierSphericalShellMaxNA=zeros(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);  %Preallocate the matrices to increase speed.

for k=1:MATRIX_SIZE
    
    for l=1:MATRIX_SIZE
        
        for j=1:MATRIX_SIZE
            
            x=k-MATRIX_SIZE/2;
            
            y=l-MATRIX_SIZE/2;
            
            z=j-MATRIX_SIZE/2;
            
            q=sqrt(x^2+y^2+z^2);
            
            if q<SPHERICAL_SHELL_RADIUS && q>SPHERICAL_SHELL_RADIUS-3
                
                r2=sqrt(x^2+y^2);
                
                alpha=atan2(r2,z);
                
                if alpha<THETA_MAX*pi/180;
                   
                    gaussianFourierSphericalShellMaxNA(k,l,j)=1;
               
                end
                                
            end
            
        end
        
    end
    
end

clear z y x r2 r q j k l

%% Paint the Fourier Shell of the Low-NA Gaussian beam.

gaussianFourierSphericalShellMinNA=zeros(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);  %Preallocate the matrices to increase speed.

for k=1:MATRIX_SIZE
    
    for l=1:MATRIX_SIZE
        
        for j=1:MATRIX_SIZE
            
            x=k-MATRIX_SIZE/2;
            
            y=l-MATRIX_SIZE/2;
            
            z=j-MATRIX_SIZE/2;
            
            q=sqrt(x^2+y^2+z^2);
            
            if q<SPHERICAL_SHELL_RADIUS && q>SPHERICAL_SHELL_RADIUS-3
                
                r2=sqrt(x^2+y^2);
                
                alpha=atan2(r2,z);
                
                if alpha<THETA_MIN*pi/180;
                   
                    gaussianFourierSphericalShellMinNA(k,l,j)=1;
               
                end
                
            end
            
        end
        
    end
    
end



%Essentially, we are creating a sphere of radius 1/lambda with a finite thickness, and then
%painting a section of that sphere, determined by the numerical aperture of
%the system in place.  Shell thickness of 3, which is related to the coherence length of the laser.  Smaller as the light becomes more monochromatic.
%Anything smaller than 20, relative to the z-axis, is included.  
%NA=2*Sin(Alpha)

clear z y x r2 r q j k l

%%
%Take the spherical shell that you created previously, perform a
%multidimensional FFT on it, and then shift the Fourier transform to the DC
%component.  Specify your imaging type here; 1=1-photon, 2=2-photon.

GaussianEFieldMaxNA=fftshift(fftn(gaussianFourierSphericalShellMaxNA)); 

GaussianEFieldMinNA=fftshift(fftn(gaussianFourierSphericalShellMinNA));

GaussianDifField=GaussianEFieldMaxNA-GaussianEFieldMinNA;



%% Bessel-Gauss, Single Point Illumination.

%Returns array the same size as the input (B), but with all singleton dimensions removed...
%A singleton only has one element.
%Figure; imshow(abs(B(:,:,XXX)),[]) --> %Here, XXX is the z-height of the
%Gaussian beam as it goes propagates in space.


plotBessel3D(GaussianDifField,MATRIX_SIZE,IMAGING_MODE,IMAGE_FIELD_OF_VIEW)

%% Comb illumination.

%COMB_WIDTH=50E-6; %Width of comb, in meters.

%COMB=zeros(size(GaussianDifField));




