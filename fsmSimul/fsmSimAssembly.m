function [outM, outCtrls] = fsmSimAssembly(varargin)
%FSMASSEMBLY simulates the assembly of a speckling actin mesh
%
% SYNOPSIS [outM, outCtrls] = fsmSimAssembly(varargin)
%
% varargin = {iniM, inCtrls}
%
% INPUT  iniM  : data structure of meshwork with fields (optional)
%                .img  : pixelated image of meshwork (typically from a prec. run; 
%                        default: 30 X 30 pixel images will be generated)
%                .mesh : meshwork (can have higher resolution), (typically from a prec. run)
%                .fpe  : distribution of free plus ends (by setting this field []
%                        a new random distribution is generated with a mean
%                        specified by acFPE (see inCtrls)
%                .fme  : distribution of free minus ends (by setting this field []
%                        a new random distribution is generated with a mean
%                        specified by acFME (see inCtrls)
%                .res  : resolution increase between .mesh and .img (must be >=1) (default: 1)
%                .pix  : pixel size in object space (um) (default: 0.05 um)
%                .hMesh: mean height of the meshwork (e.g. the height of a lamellum) (default: 0.17 um)
%                .NA   : numerical aperture of the microscope (default: 1.4)
%                .nElpNlabel : how many electrons per label are generated? (default: 255; cf. calculation 
%                                                                           protocol 29-5-01)
%                .fWC  : full well capacity of the camera (default: 30000; specs of ORCA2)
%                .nBit : number of bits in the data (default: 14)
%                .sDn  : sigma of dark noise (in grayvalues) (default: 4.5; camera calibration mai 2001)
%                .LR   : labeling ratio (default: 0.005)
%                .nSpeck : vector with speckle number (typically from a prec. run)
%                .time   : time points at which .nSpeck has been calculated (typically from a prec. run)
%                .writeFields : 1 if all the fields should be stored in files 
%                               0 if not (default)
%
%        inCtrls : data structure with chemical controls (optional)
%                 .tSpan: period in [s] for which the process is excecuted 
%                 .tMinStep : minimum time step in [s] for calculation (if not set,
%                             the program determines the time step such that not more
%                             than 100 monomers are added or removed on average
%                 .cM  : concentration of free monomer
%                 .kPon: association rate at free plus ends [ # / (uM s)] (default: 11.6)
%                 .kPoff:dissociation rate at free plus ends [ # / s] (default: 1.4)
%                 .kMon: association rate at free minus ends (default: 0)
%                 .kMoff:dissociation rate at free minus ends (default: cM*kPon - kPoff, i.e. perfect 
%                                                                       treadmilling)
%                 .acFPE:area concentration of free plus ends [ # / um^2 ] (default: 50)
%                 .acFME:area concentration of free minus ends [ # / um^2 ] (default: 20)
%
%         Both input parameters are optional and can be filled with whatever fields are necessary



global MESH CHEM PSF INRATE OUTRATE;

% check for the availbility and correctness of input variables
fillOptions(varargin);

% calculate PSF
PSF = psf2D(MESH.pix,MESH.NA,0.5); 

% if there are no fields with plus and minuse free ends, generate a randon field 
if isempty(MESH.fpe)
   MESH.fpe = (rand(size(MESH.mesh))+0.5)*CHEM.acFPE*(MESH.pix/MESH.res)^2;
end;
if isempty(MESH.fme)
   MESH.fme = (rand(size(MESH.mesh))+0.5)*CHEM.acFME*(MESH.pix/MESH.res)^2;
end;

% calculate INRATE and OUTRATE fields
INRATE = ( CHEM.kPon * MESH.fpe + CHEM.kMon  * MESH.fme ) * CHEM.cM;
OUTRATE = ( CHEM.kPoff * MESH.fpe + CHEM.kMoff  * MESH.fme);

% prepare parameters of the ODE solver
% the time step should not be larger than the time needed to add or remove on 
% average 100 labelled monomers to the entire frame
tStep = 100 / ( abs( sum(INRATE(:)) - sum(OUTRATE(:)) ) * MESH.LR);
if(isfield(CHEM,'tMinStep'))
   if tStep < CHEM.tMinStep
      tStep = CHEM.tMinStep;
   end;
end;

% options = odeset('RelTol',1,'AbsTol',1,'InitialStep',tStep);
tCtrl = [0:tStep:CHEM.tSpan]+MESH.time(end);

% run the odesolver (usage not appropriate because the solver can also 
% make steps in negative direction, which are non-deterministic
% [MESH.time,MESH.nSpeck] = ode23('fsmSimAssemblyODEFun',tCtrl,initialNSpeck);
for i = 2:length(tCtrl)
   dN = feval('fsmSimAssemblyODEFun',tCtrl(i),0);
   MESH.nSpeck(end+1) = MESH.nSpeck(end)+dN;
   MESH.time(end+1) = tCtrl(i);
   disp(sprintf('runtime = %6.2f; N = %6d; dN = %6d; incub.time  = %6.2f',...
      tCtrl(i)-tCtrl(1),MESH.nSpeck(end),dN,tCtrl(i)));
end;

% copy the results to the output structure
outM = MESH;

% copy the chemical controls to the output structure
outCtrls = CHEM;

clear global MESH CHEM PSF INRATE OUTRATE;   
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local service function

function fillOptions(inputs)

global MESH CHEM

pMeshFields = {'img','mesh','fpe','fme','res','pix','hMesh','NA','nElpNlabel','fWC','nBit',...
      'sDn','LR','nSpeck','time','writeFields'};
pChemFields = {'tSpan','cM','kPon','kPoff','kMon','kMoff','acFPE','acFME'};

% default mesh parameter section 
dim = 30;           % sidelength of simulated area in pixels
pix = 0.05;         % pixel size in um
res = 1;            % resolution increase between meshwork and pixel image
hMesh = 0.170;      % mean height of the meshwork (e.g. the height of a lamellum)
NA  = 1.4;          % numerical aperture of the microscope
LR = 0.005;         % labelling ratio
nElpNlabel = 255;   % see calculation protocol 29-5-01
fWC = 30000;        % specs for ORCA2 camera 
nBit = 14;
sDn  = 4.5;         % dark noise in greyvalues

writeFields = 0;    % fields are not stored

% default chemical parameter section 
tSpan =  120;       % incubation time for which the simulation is run under the conditions 
                    % defined by the parameters set below (default: 120 s)
cM = 25;            % concentration of free monomers [ uM ]
kPon = 11.6;        % association rate at free plus ends [ # / (uM s)] Previous: 11.3
kPoff = 1.4;        % dissociation rate at free plus ends [ # / s]     Previous: 1.6
kMon = 0;           % association rate at free minus ends [ # / (uM s)]
kMoff = cM*kPon - kPoff;     
% dissociation rate at free minus ends (for the moment we assume perfect 
% treadmilling, i.e. on average the material that is incorporated at the plus end
% will fall out at the minus end; assembly / disassembly in the meshwork 
% is only controlled by the ratio free plus ends / free minus ends

acFPE = 50;     % area concentration of free plus ends [ # / um^2 ]
acFME = 20;     % area concentration of free minus ends [ # / um^2 ]

% is any of the inputs a mesh control structure ?
i = 1;
MESH = [];
while (i <= length(inputs)) & isempty(MESH)
   for j = 1:length(pMeshFields)
      fFields(j)=any(strcmp(pMeshFields(j),fieldnames(inputs{i})));
   end
   if any(fFields)
      MESH = inputs{i};
   end;
   i = i + 1;
end;
% is any of the inputs a chem control structure ?
i = 1;
CHEM = [];
fFields = [];
while (i <= length(inputs)) & isempty(CHEM)
   for j = 1:length(pChemFields)
      fFields(j)=any(strcmp(pChemFields(j),fieldnames(inputs{i})));
   end;
   if any(fFields)
      CHEM = inputs{i};
   end;
   i = i + 1;
end;

% check the validity of the CHEM structure
if ~isfield(CHEM,'tSpan')
   CHEM.tSpan = tSpan;
end;

if ~isfield(CHEM,'cM')
   CHEM.cM = cM;
end;

if ~isfield(CHEM,'kPon')
   CHEM.kPon = kPon;
end;

if ~isfield(CHEM,'kPoff')
   CHEM.kPoff = kPoff;
end;

if ~isfield(CHEM,'kMon')
   CHEM.kMon = kMon;
end;

if ~isfield(CHEM,'kMoff')
   CHEM.kMoff = CHEM.cM*CHEM.kPon - CHEM.kPoff;
end;

if ~isfield(CHEM,'acFPE')
   CHEM.acFPE = acFPE;
end;

if ~isfield(CHEM,'acFME')
   CHEM.acFME = acFME;
end;


% check the validity of the MESH structure
if ~isempty(MESH)
   % check resolution ratio between mesh and img first
   if ~isfield(MESH,'res')
      MESH.res = res;
   else
      MESH.res = round(MESH.res);
      if MESH.res < 1
         MESH.res = res;
      end;
   end;
   % the time field is the most important, .nSpeck has to be of the same length
   % and if the time field exists, .img and .mesh represent the current label distribution 
   % and image thereof.
   if ~isfield(MESH,'time')
      MESH.nSpeck = 0;
      MESH.time = 0;
      MESH.img = zeros(dim);
      MESH.mesh = zeros(size(MESH.img)*MESH.res);
   else       
      if ~isfield(MESH,'nSpeck') 
         error('missing .nSpeck field associated with input .time field')
      end;
      if (length(MESH.time) ~= length(MESH.nSpeck))
         error('.time and .nSpeck are incompatible');
      end;
      if ~isfield(MESH,'mesh')
         error('a time point > 0 is entered without .mesh field');
      end;
      if ~isfield(MESH,'img') 
         MESH.img = zeros(size(MESH.mesh)/MESH.res);
      else
         if (any((size(MESH.img)*MESH.res)~=size(MESH.mesh)))
            error('.img , .mesh , and .res fields are incompatible');
         end;
      end;
   end;
   
   % check the free end fields, and generate new ones if necessary
   if ~isfield(MESH,'fpe')
      MESH.fpe = [];
   else
      if ~isempty(MESH.fpe)
         if(any(size(MESH.fpe)~=size(MESH.mesh)))
            error('.fpe and .mesh fields are incompatible');
         end;
      end;
   end;
   
   if ~isfield(MESH,'fme')
      MESH.fme = [];
   else
      if ~isempty(MESH.fme)
         if(any(size(MESH.fme)~=size(MESH.mesh)))
            error('.fme and .mesh fields are incompatible');
         end;
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
   
   if ~isfield(MESH,'writeFields')
      MESH.writeFields = writeFields;
   end;
      
else
   % the structure has to be filled with default values
   MESH.nSpeck = 0;
   MESH.time = 0;
   MESH.res = res;
   MESH.pix = pix;
   MESH.img = zeros(dim);
   MESH.mesh = zeros(size(MESH.img)*MESH.res);
   MESH.fpe = [];
   MESH.fme = [];
   MESH.NA = NA;
   MESH.hMesh = hMesh;
   MESH.LR = LR;
   MESH.nElpNlabel = nElpNlabel;
   MESH.fWC = fWC;
   MESH.nBit = nBit;
   MESH.sDn = sDn;
   MESH.writeFields = writeFields;
end;