function dn = fsmSimAssemblyODEFun(t,n)
%fsmSimAssemblyODEFun auxiliary function for fsmSimAssembly

global MESH PSF ARATE CFM LASTTIME LASTN NPASSES;

dT = t - LASTTIME;
nMono = round(ARATE * dT * length(MESH.img)*MESH.pix^2);


if nMono > 0
   % folding two MC steps into each other to distribute labels
   
   % 1.) generate a distribution of new monomers to be added
   monoField = rand(size(MESH.mesh));            
   monoField = monoField / sum(monoField(:));
   
   % 2.) generate a string of labelled and unlabelled new monomers that will be incorporated
   newMono = (rand(nMono,1) < MESH.LR);
   
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
   lbs = reshape(diff(cBin(aSpace)), size(MESH.mesh));
   MESH.mesh = MESH.mesh + lbs;
   
   % add a second field with free monomer concentration 
   fmLbs = freeMonomerDensity(MESH,CFM);
   tempMesh = MESH.mesh + fmLbs;
   
   % generate image
   % convolve with PSF
   hResLbDistr = conv2(tempMesh,PSF,'same');
   
   % pixelation
   pixLbDistr = pixelate(hResLbDistr,MESH.res);
      
   % add a noise model here 
   pixNel = pixLbDistr*MESH.nElpNlabel;
   noisyPixNel = poissrnd(pixNel);   % add poisson distributed noise
   
   % generate an additive dark noise field with zero mean, 
   % MESH.fWC * (MESH.sDn / 2^MESH.nBit)
   noisyPixNel = noisyPixNel + ...
      randn(size(noisyPixNel))*MESH.fWC * (MESH.sDn / 2^MESH.nBit);
   
   % clip negative values (can happen for low electron numbers)
   noisyPixNel(find(noisyPixNel < 0)) = 0;
   
   % scale the electron field to an image
   MESH.img = noisyPixNel / MESH.fWC;

   
   % Display
   surf(MESH.img); axis([0 30 0 30 0 0.05]);
   pause(1);
   refresh;
   
   % add speckle detector here
   %I=locmax2D(MESH.img,[5 5]);            % To be replaced!
   %nI=I>0; nb=length(find(nI));
   
   % dn
   %dn=nb-LASTN;
   %LASTN= nb;
   
   dn = 0;
else   
   dn = 0;
end;

LASTTIME = t;
NPASSES = NPASSES + 1;

if MESH.writeFields
   
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% local service functions

function lbs = freeMonomerDensity(MESH,cM)

% Constant definition: Avogadro's number Na
Na=6.02e23;

% cMono = number of actin monomers per unit volume (MESH.pix^2*h)
cMono=round((cM*1e-6)*(Na/1e15)*prod(size(MESH.img))*MESH.pix^2*MESH.hMesh);

% 1) Generate a distribution of new monomers (in space)
monoField = rand(size(MESH.mesh));
monoField = monoField/sum(monoField(:));

% 2) Generate a string of labelled and unlabelled monomers
newMono = (rand(cMono,1) < MESH.LR);

% 3) Generate an address space
aSpace = [0;round(cumsum(monoField(:))*cMono)]+1;

% 4) Labels
cBin = [0; cumsum(newMono)];

% 5) Reshape labels
lbs=reshape(diff(cBin(aSpace)), size(MESH.mesh));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%