function dn = fsmSimAssemblyODEFun(t,n)
%fsmSimAssemblyODEFun auxiliary function for fsmSimAssembly

global MESH CHEM PSF INRATE OUTRATE;

dT = t - MESH.time(end);
nMonoIn = round(sum(INRATE(:) * dT));
nMonoOut = round(sum(OUTRATE(:) * dT));
   
% add the monomers that polymerize
tempLbl = lblField(INRATE,nMonoIn,MESH.LR);
MESH.mesh = MESH.mesh + tempLbl;
% remove the monomers that drop out
tempLbl = lblField(OUTRATE,nMonoOut,MESH.LR);
MESH.mesh = MESH.mesh - tempLbl;

% clip negative label values (this means that locally all labels have been
% fallen out of the meshwork -> no signal produced at this location
MESH.mesh(find(MESH.mesh < 0)) = 0;

% add a second field with free background monomer concentration 
nFreeMono = round((CHEM.cM*1e-6)*(6.02e23/1e15)...
   *prod(size(MESH.img))*MESH.pix^2*MESH.hMesh);

fmLbs = lblField(ones(size(MESH.mesh)),nFreeMono,MESH.LR);
tempMesh = MESH.mesh + fmLbs;

% generate image
% convolve with PSF
hResLbDistr = conv2(tempMesh,PSF,'same');

% pixelation
pixLbDistr = pixelate(hResLbDistr,MESH.res);

% add a noise model here 
pixNel = pixLbDistr*MESH.nElpNlabel;
noisyPixNel = poissrnd(pixNel);   % add poisson distributed noise

if MESH.writeFields
   pNoise = noisyPixNel - pixNel;
end;

% generate an additive dark noise field with zero mean, 
% MESH.fWC * (MESH.sDn / 2^MESH.nBit)
dNoise = randn(size(noisyPixNel))*MESH.fWC * (MESH.sDn / 2^MESH.nBit);
noisyPixNel = noisyPixNel + dNoise;

% clip negative values (can happen for low electron numbers)
noisyPixNel(find(noisyPixNel < 0)) = 0;

% scale the electron field to an image
MESH.img = noisyPixNel / MESH.fWC;

% clip values larger than 1 (saturated images)
MESH.img(find(MESH.img > 1)) = 1;

% Display
% surf(MESH.img); axis([0 30 0 30 0 1]);
% refresh;

% add speckle detector here
parameters=[2, MESH.sDn/2^MESH.nBit, 2e-4]; % TEMP
[nb,speckMap] = speckleDetector(Gauss2D(MESH.img,1),parameters);

% display speckle map 
% imshow(MESH.img,[]);
% imOverlayPointMap(speckMap,'r.');

% dn
dn = nb - MESH.nSpeck(end);

if MESH.writeFields
   timePoint = t;
   lblMap = MESH.mesh;
   img = MESH.img;
   indString = sprintf('%.6d',length(MESH.time));
   save(strcat('fsmSimMaps',indString),...
      'timePoint','lblMap','img','pixNel','pNoise','dNoise',...
      'nb','speckMap');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% local service functions

function lbs = lblField(rateField,nMono,lr)

% folding two MC steps into each other to distribute labels
   
% 1.) generate a distribution of new monomers to be added
%     The rate field acts as a weight for the spatial distribution   
monoField = rand(size(rateField)).*rateField;
monoField = monoField / sum(monoField(:));

% 2.) generate a string of labelled and unlabelled new monomers that will be incorporated
newMono = (rand(nMono,1) < lr);

% 3.) based on the monomer distribution, generate an address space into 
%     string of new monomers:
%     cumsum(mesh(:))*nMono defines the number of monomers that have to 
%     be included up to the coumn-wise ordered mesh field; the initial 0 added
%     accounts for the starting point of monomer density integration, i.e. 
%     diff(aSpace) would immediately tell how many monomners are to be placed. 
%     Add +1 to make  sure that the addresses start with 1 not violating the address conventions
%     of MATLAB
aSpace = [0;round(cumsum(monoField(:))*nMono)]+1;

% 4.) add an inital no-label to the label string and take the cumsum of it; 
%     analogous to step 3, this means that difference of cBin tells what kind
%     of label is placed
cBin = [0; cumsum(newMono)];

% 5.) grab the labels based on the positions in aSpace and take the difference;
%     reshape the column vector back to a the dimensions of the mesh and
%     add the new labels to existing mesh
lbs = reshape(diff(cBin(aSpace)), size(rateField));