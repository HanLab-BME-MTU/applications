function [outM] = fsmSimAssembly(iniM)
%FSMASSEMBLY simulates the assembly of a speckling actin mesh
%
% SYNOPSIS
%
% INPUT  iniM  : data structure of meshwork with fields (optional)
%                .img  : pixelated image of meshwork
%                .mesh : meshwork (can have higher resolution)
%                .res  : resolution increase between .mesh and .img (must be >=1)
%                .pix  : pixel size in object space (um)
%                .hMesh: mean height of the meshwork (e.g. the height of a lamellum)
%                .NA   : numerical aperture of the microscope
%                .nElpNlabel : how many electrons per label are generated?
%                .fWC  : full well capacity of the camera
%                .nBit : number of bits in the data
%                .sDn  : sigma of dark noise (in grayvalues)
%                .LR   : labeling ratio
%                .nSpeck : vector with speckle number (typically from a previous run)
%                .time   : time points at which .nSpeck has been calculated
%                .writeFields : 1 if all the fields should be stored in files 
%                               0 if not (default)

global MESH PSF ARATE CFM LASTTIME LASTN NPASSES;

% default control parameter section 
dim = 30;           % sidelength of simulated area in pixels
pix = 0.05;         % pixel size in um
res = 1;            % resolution increase between meshwork and pixel image
hMesh = 0.170;      % mean height of the meshwork (e.g. the height of a lamellum)
NA  = 1.4;          % numerical aperture of the microscope
tSpan =  600;       % time span of period over which the development is observed
nElpNlabel = 255;   % see calculation protocol 29-5-01
fWC = 30000;        % specs for ORCA2 camera 
nBit = 14;
sDn  = 4.5;         % dark noise in greyvalues

writeFields = 0;    % fields are not stored

% chemical parameter section 

LR = 0.005;      % labelling ratio

cM = 25;         % concentration of free monomers [ uM ]
kPon = 11.6;     % association rate at free plus ends [ # / (uM s)] Previous: 11.3
kPoff = 1.4;     % dissociation rate at free plus ends [ # / s]     Previous: 1.6
kMon = 0;        % association rate at free minus ends [ # / (uM s)]
kMoff = cM*kPon - kPoff;     
% dissociation rate at free minus ends (for the moment we assume perfect 
% treadmilling, i.e. on average the material that is incorporated at the plus end
% will fall out at the minus end; assembly / disassembly in the meshwork 
% is only controlled by the ratio free plus ends / free minus ends

acFPE = 50;     % area concentration of free plus ends [ # / um^2 ]
acFME = 20;     % area concentration of free minus ends [ # / um^2 ]

% Initialize global variable LASTN, which stores the number of speckles of 
% the previous iteration of ode23
LASTN=0;

% check the incoming data structure, if there is any
if nargin == 1
   MESH = iniM;
   if ~isfield(MESH,'img')
      MESH.img = zeros(dim);
   end;
   if ~isfield(MESH,'res')
      MESH.res = res;
   else
      MESH.res = round(MESH.res);
      if MESH.res < 1
         MESH.res = res;
      end;
   end;
   if ~isfield(MESH,'mesh')
      MESH.mesh = zeros(size(MESH.img)*MESH.res);
   else
      if(any((size(MESH.img)*MESH.res)~=size(MESH.mesh)))
         MESH.img = zeros(size(MESH.img));
         MESH.mesh = zeros(size(MESH.img)*MESH.res);
      end;
   end;
      
   if ~isfield(MESH,'NA')
      MESH.NA = NA;
   end;
   
   if ~isfield(MESH,'nElpNlabel ')
      MESH.nElpNlabel = nElpNlabel ;
   end;
   
   if ~isfield(MESH,'fWC')
      MESH.fWC = fWC;
   end;
   
   if ~isfield(MESH,'nBit')
      MESH.nBit = nBit;
   end;
   
   if ~isfield(MESH,'sDn')
      MESH.sDn = sDn;
   end;
   
   if ~isfield(MESH,'pix')
      MESH.pix = pix;
   end;
   
   if ~isfield(MESH,'hMesh')
      MESH.hMesh = hMesh;
   end;

   if ~isfield(MESH,'LR')
      MESH.LR = LR;
   end;
   
   if ~isfield(MESH,'time')
      tCtrl = [0 tSpan];
      LASTTIME = 0;
      initialNSpeck = 0;
      NPASSES = 0;
   else 
      tCtrl = [0 tSpan] + iniM.time(end);   % set the start time as the old stop time
      LASTTIME = iniM.time(end);
      initialNSpeck = iniM.nSpeck(end);
      % chop off the last element; it defines the initial state of the new run and 
      % will be concatenated again at the end.
      prevTime = iniM.time(1:end-1);     
      prevNSpeck = iniM.nSpeck(1:end-1);
      NPASSES = length(iniM.time);
   end;
   
   if ~isfield(MESH,'writeFields')
      MESH.writeFields = writeFields;
   end;
      
else
   MESH.img = zeros(dim);
   MESH.mesh = zeros(size(MESH.img)*MESH.res);
   MESH.res = res;
   MESH.NA = NA;
   MESH.pix = pix;
   MESH.hMesh = hMesh;
   MESH.LR = LR;
   MESH.nElpNlabel = nElpNlabel;
   MESH.fWC = fWC;
   MESH.nBit = nBit;
   MESH.sDn = sDn;
   MESH.writeFields = writeFields;
   
   % prepare the parameters for the ODE solver
   tCtrl = [0 tSpan];
   LASTTIME = 0;
   initialNSpeck = 0;
   prevTime = [];
   prevNSpeck = [];
   NPASSES = 0;
end;

% calculate PSF
PSF = psf2D(MESH.pix,MESH.NA,0.5); 

% calculate the assembly rate 
ARATE = ( kPon * acFPE + kMon  * acFME ) * cM -  ( kPoff * acFPE + kMoff  * acFME );

% copy the concentration of free monomers and meshwork height into a global variables
CFM = cM;

% prepare parameters of the ODE solver
% the time step should not be larger than the time needed to add on 
% average 1000 monomers to the entire frame
tStep = 1000 / (ARATE * length(MESH.img)*MESH.pix^2);
options = odeset('RelTol',1,'AbsTol',1,'InitialStep',tStep);
tCtrl = [tCtrl(1):tStep:tCtrl(2)];

% run the odesolver
[MESH.time,MESH.nSpeck] = ode23('fsmSimAssemblyODEFun',tCtrl,initialNSpeck);

% copy the results to the output structure
outM = MESH;

% append the new time points to the incoming ones
outM.time = [prevTime ; MESH.time];
outM.nSpeck = [prevNSpeck ; MESH.nSpeck];

   
clear global MESH PSF ARATE CFM HMESH LASTTIME LASTN NPASSES;   