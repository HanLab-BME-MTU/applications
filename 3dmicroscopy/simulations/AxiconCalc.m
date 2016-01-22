;clc
clear

%% Determine the refractive index for UV fused silica.

lambda=.561; %In microns.

p1=(0.6962.*lambda^2)./((lambda.^2)-0.0047);
p2=(0.4079.*lambda^2)./((lambda.^2)-.0135);
p3=(0.8975.*lambda^2)./((lambda.^2)-97.9340);

refractiveIndex=sqrt(p1+p2+p3+1)

%% Determines the depth of field of the axicon immediately following the axicon.

inputBeamDiameter=1.8;

ringWidth=inputBeamDiameter./2;

alpha=[1 2 5 10]

DepthOfField=(0.5*inputBeamDiameter)./((refractiveIndex-1).*alpha)




%% Determines the focal length of the lens following the axicon necessary for creating a 9.8 mm O.D. ring.


AnnulusOuterDiameter=9.8E-3;

focalLength=AnnulusOuterDiameter./(2.*tand((refractiveIndex-1).*alpha))



%% Determines the beam waist at the annulus, assuming standard trig relations.

beamWaist=(4.*(lambda./1000)./pi).*(focalLength./ringWidth)
 




%%

fLength=7E-3;
fTheta=32.2;

diameter=fLength.*4.*(lambda./1000)./(2.*pi.*beamWaist);

yPropagationLength=diameter/(sind(fTheta).*1000)
