function objImage = generateSyntheticImageGeometry(varargin)


%% ---------------- Input ------------------- %%

ip = inputParser;
ip.FunctionName = mfilename;

%Image parameters
ip.addParamValue('ImageSize',[256 256 256],@(x)(numel(x)==3 && all(x>1)));
ip.addParamValue('PixelSizeXY',353.1e-9,@(x)(numel(x)==1 && x > 0 && x < 1e-5));%Make sure they input it in meters!
ip.addParamValue('PixelSizeZ',500e-9,@(x)(numel(x)==1 && x > 0 && x < 1e-5));%Make sure they input it in meters!
ip.addParamValue('OverSampling',3,@(x)(numel(x)==1 && isposint(x)));%Factor by which to over-sample each voxel to avoid aliasing

%Genearal object parameters
ip.addParamValue('ObjectType','Sphere',@(x)(ischar(x)));%Name of object type to create
ip.addParamValue('ObjectRadius',20e-6,@(x)(numel(x)==1 && x > 0 && x < 1e-3));
ip.addParamValue('ObjectLength',20e-6,@(x)(numel(x)==1 && x > 0 && x < 1e-3));
ip.addParamValue('ObjectOrientation','z',@(x)(ischar(x) && any(strcmpi(x,{'x','z'}))));
ip.addParamValue('MembraneThickness',[],@(x)(numel(x)==1 && x > 0 && x < 1e-5));%Thickness of "membrane" (object wall) in meters

ip.parse(varargin{:});
p = ip.Results;

%If not input, set the membrane thickness based on the final voxel size.
if isempty(p.MembraneThickness)
    p.MembraneThickness = max(p.PixelSizeXY,p.PixelSizeZ);
end

pixAspect = p.PixelSizeZ / p.PixelSizeXY;%pixel aspect ratio
radInVox = p.ObjectRadius / p.PixelSizeXY;%Object radius in XY voxels
lenInVox = p.ObjectLength / p.PixelSizeXY;%Object length in XY voxels
thickInVox = p.MembraneThickness / p.PixelSizeXY;%membrane thickness in voxels

%% ------------- Parameters ---------------- %%

% %Scale factor so we can use integers and still have non-cuboidal voxels
% scFact = 1e2;

%% --------------- Geometry Creation --------------- %%
%Create the higher-resolution geometry which will be down-sampled to the
%final image size

%Create coordinate matrices, accounting for the origin and using single to save memory.

% [X,Y,Z] = meshgrid(single(-(p.ImageSize(2)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(2)/2),...
%                    single(-(p.ImageSize(1)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(1)/2),...
%                    single(-(p.ImageSize(3)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(3)/2));
                              
[X,Y,Z] = meshgrid(single(-(p.ImageSize(2)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(2)/2) * p.PixelSizeXY,...
                   single(-(p.ImageSize(1)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(1)/2) * p.PixelSizeXY,...
                   single(-(p.ImageSize(3)/2 - 1/p.OverSampling):1/p.OverSampling:p.ImageSize(3)/2) * p.PixelSizeZ);


switch p.ObjectType


   case 'Sphere'
              
        objImage = X .^ 2 + Y .^2 + Z .^2  > p.ObjectRadius ^2 & ...
                   X .^ 2 + Y .^2 + Z .^2  <= (p.ObjectRadius+p.MembraneThickness) ^2;
                              
        
    case 'Cylinder'
        
        switch p.ObjectOrientation
            
            case 'x'
                
                objImage = (Z .^ 2 + Y .^2 > p.ObjectRadius ^2 & ... %The walls of the cylinder
                           Z .^ 2 + Y .^2 <= (p.ObjectRadius+p.MembraneThickness) ^2 & ...
                           abs(X) < p.ObjectLength ) | ...
                           (Z .^ 2 + Y .^2 <= (p.ObjectRadius+p.MembraneThickness) ^2 & ... % and the caps
                           abs(X) > p.ObjectLength & abs(X) < (p.ObjectLength + p.MembraneThickness));
        
                
            case 'z'
               
                objImage = (X .^ 2 + Y .^2 > p.ObjectRadius ^2 & ... %The walls of the cylinder
                           X .^ 2 + Y .^2 <= (p.ObjectRadius+p.MembraneThickness) ^2 & ...
                           abs(Z) < p.ObjectLength ) | ...
                           (X .^ 2 + Y .^2 <= (p.ObjectRadius+p.MembraneThickness) ^2 & ... % and the caps
                           abs(Z) > p.ObjectLength & abs(Z) < (p.ObjectLength + p.MembraneThickness));
                       
                        
                
                
            otherwise 
                error('Object orientation must be "x" or "z"')
               
        end
                        
   otherwise
       error([ p.ObjectType ' is not a supported object type!'])
end

%Free up the memory of the coord matrices
clear X Y Z

%% --------------- Down-Sampling ---------------- %%


%This is not the ideal way to do this but it'll work: Local-mean filter
%before down-sampling.
if p.OverSampling > 1
    objImage = imfilter(single(objImage),ones(ones(1,3) * p.OverSampling) / (p.OverSampling^3));
    objImage = objImage(1:p.OverSampling:end,:,:);
    objImage = objImage(:,1:p.OverSampling:end,:);
    objImage = objImage(:,:,1:p.OverSampling:end);
end







               