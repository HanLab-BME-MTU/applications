function fsmDemo()
%FSMDEMO simulation of fluorescent speckling
%
% SYNOPSIS fsmDemo()
%
% INPUT screen : screen settings; supported are
%                'normal' (default)
%                'vaio'   (small SONY Vaio LapTop)

% STARTED GD 27-Nov-1999

% global variable 
global WVL__;
WVL__ = 500;
global FIRSTTIME__
FIRSTTIME__ = 1;

% start control and graphics window
% if nargin == 0
%   screen = 'normal';
% end;

% switch lower(screen)
% case 'normal'
%    ctrlWinH = fsmDemoCtrlWin;
% case 'vaio'
%    ctrlWinH  = fsmDemoCtrlWinVaio;
% otherwise 
%    error('unknown screen setting');
% end;

ctrlWinH = fsmDemoCtrlWin;

fH = ctrlWinH;
graphWinH = figure('Units','points', ...
	'MenuBar','none', ...
	'Name','FSM Demo', ...
	'NumberTitle','off', ...
   'Resize','off');
ctrlPos = get(ctrlWinH,'Position');
graphPos = get(graphWinH,'Position');

% switch lower(screen)
% case 'normal'
% case 'vaio'
%    ctrlPos = [10, 10, 0, 0];
%    
%    
% end;

% dependent on the resolution of the projector, one has to shift windows up
screenSize = get(0,'Screensize');  % 0 is the address of the root; for resolutions less than the graphics
                                   % card, the first 2 numbers specifify the offset   
ctrlPos(1) = ctrlPos(1)-screenSize(1)+1; 
ctrlPos(2) = ctrlPos(2)-screenSize(2)+1; 
set(ctrlWinH,'Position',ctrlPos);
set(graphWinH,'Position',...
   [ctrlPos(1),ctrlPos(2)+ctrlPos(4)+20,graphPos(3:4)]);

firstTime = 1;

psf2D=[];
winPosPolym=[ctrlPos(1)+graphPos(3)+20,ctrlPos(2)+ctrlPos(4)+20,20,graphPos(4)];
nTimes = 1;
nEl=0;
nPix = 0;
NA = 0;
lab = 0;
meanSN = 0;
while(fH == ctrlWinH)
   
   % get parameters
   oldParams = [nEl,nPix,NA,lab];
   [nEl,nPix,mSze,pSze,NA,wvl,lab,camRng,nElPerLabel,fwc,gain,nse,speed] = getParams(ctrlWinH,meanSN);
   chg = any(([nEl,nPix,NA,lab] - oldParams)~=0);
   
   % disp(sprintf('%f %d',NA,nPix));
   if(chg)
      lastChg = nTimes;
   end;
   
   % add nPix monomers at the front end
   if(firstTime)
      polym = zeros(nEl,1);
   end;
   polym = [sum(rand(nPix,mSze(2))<lab,2);polym(1:end-nPix)];
   
   % convolve with PSF and sample in pixels
   
   % 1.Step: create PSF kernel
   %         defintion between -/+ the 2nd root of the Bessel function 
   y = -7.02*(1.22 * wvl/(2*NA))/3.83 + 0.01 : mSze(1) : ...
      7.02*(1.22 * wvl/(2*NA))/3.83;
   ys = y /(1.22 * wvl/(2*NA))*3.83;
   psfs = (besselj(1,ys)./ys);
   psf = psfs.*conj(psfs);
   psf = psf / sum(psf);  % make sure that the optics is an energy 
                          % preserving filter

   % 2. Step: convolution
   picPolym = conv(polym,psf);
   
   % remove border area first
   hspsf = floor(length(psf)/2);
   picPolym = picPolym(hspsf+1:end-hspsf);
   
   % Show polymer
   % $$$$$$$$$$$$$$$ Comment out here if 2D view is superfluous
   % psf2D=fsm2dPolym(firstTime,psf,polym,nPix,psf2D,winPosPolym);
   
   % generate pixelated picture
   pixBdIndx = 1:nPix:length(picPolym);
   cumPicPolym = cumsum(picPolym);
   clear pixPicPolym;
   for(i=2:length(pixBdIndx))
      pixPicPolym(i-1) = cumPicPolym(pixBdIndx(i))- ...
         cumPicPolym(pixBdIndx(i-1));
   end;
   
   % convert the label distribution into electrons
   pixElPolym = pixPicPolym * nElPerLabel;
   
   % add dark noise and shot noise
   dNoise = randn(size(pixElPolym))*nse(1);
   
   if(nse(2))
      % shot noise contaminated 
      sncPixElPolym = poissrnd(pixElPolym);
      indxNaN = find(isnan(sncPixElPolym));
      sncPixElPolym(indxNaN) = pixElPolym(indxNaN);
      sNoise = sncPixElPolym - pixElPolym;
      meanSN = mean(abs(sNoise))/fwc; 
   else
      sNoise = zeros(size(pixElPolym));
   end;
      
   pixElNsePolym = pixElPolym + dNoise + sNoise;
   
   % calculate signal-to-noise
   lMaxEl = locmax1d(pixElPolym);
   lMinEl = locmin1d(pixElPolym);
   % mean aplitude (improve this later with a real peak height analysis)
   if isempty(lMinEl) | isempty(lMaxEl)
      snr = 0;
   else
      mAmpl = mean(pixElPolym(lMaxEl)) - mean(pixElPolym(lMinEl));
      nseEstim = sqrt(mean(sNoise.^2)+mean(dNoise.^2));
      if nseEstim > 0;
         snr = mAmpl / nseEstim;
      else
         snr = inf;
      end;
   end;
   
   % generate an image within the camera range and with an offset
   % and correct for gain
   pixImgNsePolym = gain*(pixElNsePolym / fwc *(camRng(2) - camRng(1)))+ camRng(1); 
   pixImgPolym = gain*(pixElPolym / fwc *(camRng(2) - camRng(1)))+ camRng(1); 

   
   % clip values smaller than 0 and larger than camRng2
   pixImgPolym(find(pixImgPolym < 0)) = 0;
   pixImgPolym(find(pixImgPolym > camRng(2))) = camRng(2);
   pixImgNsePolym(find(pixImgNsePolym < 0)) = 0;
   pixImgNsePolym(find(pixImgNsePolym > camRng(2))) = camRng(2);

   
   % adjust the mean intensity
%   maxPicPolym = max(pixPicPolym);
%   if(maxPicPolym > 0)
%      pixPicPolym = maxI*pixPicPolym/maxPicPolym;
%   end;
      
   % compute a contrast measure
   %cm = getContrastMeas('meanModulation',pixImgPolym);
   cm = getContrastMeas('peakModulation',pixImgPolym);
   
   % plot functions
   figure(graphWinH);
   plot(1:length(pixImgPolym),pixImgPolym,'b-',1:length(pixImgNsePolym),pixImgNsePolym,'r-');
   
   % iMax = max(pixPicPolym);
   % if(iMax<0.1)
   %    set(gca,'YLim',[0,0.1]);
   % else
   %    set(gca,'YLim',[0,ceil(10*iMax)/10]);
   % end;
   set(gca,'YLim',[0,camRng(2)]);
   set(gca,'XLim',[1,length(pixImgPolym)]);
   str = sprintf('Modulation: %6.3f   SNR: %6.3f',cm,snr);
   if( (nTimes-lastChg) > length(pixImgPolym))
      title(str,'FontWeight','bold','FontSize',12);
   else
      title(str,'FontWeight','normal','FontSize',10);
   end;
   
   %refresh(graphWinH);
   %refresh(ctrlWinH);
   wait(speed);
   
   fH = findobj(0,'Type','figure','Tag','figFSMDemo');
   firstTime = 0;
   nTimes = nTimes + 1;
end;

clear global WVL__;
clear global FIRSTTIME__;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%

function [nEl,nPix,monomSze,polymSze,NA,wvl,lab,camRng,nElPerLabel,fwc,gain,noise,speed] = getParams(fH,sn)
% gets the control parameters from the control window

global WVL__;

% definition of key quantities (units nm)
polymSze = 5000;     % (size of the polymer)
monomSze = [7.9, 13]; % [size of the monomer (identical with size of a full helical ring), 
                      %  number of protofilaments in helix]
% pixSze = 12000;     % (pixel size in the image plane)
mag = 100;            % (microscope magnification)
% NA  = 1.4;          % (numerical aperture)
wvl = WVL__;          % (wavelength)
% lab = 5/100;        % (percentage of labelled free monomers)
camRng = [100,2^14];  % simulation of a 14 camera
nElPerLabel = 255;    % value taken from a calibration of the ORCA II
fwc = 30000;          % value taken from ORCA II specs
                      
speed = 0.1;          % (frame rate control for graphical display)


% get pixel size from slider
pixSze = get(findobj(fH,'Tag','sliderPixSze'),'Value')*1000;

% get NA from list
indx = get(findobj(fH,'Tag','listNA'),'Value');
str = get(findobj(fH,'Tag','listNA'),'String');
NA = str2num(str(indx,:));

% get mag from list
indx = get(findobj(fH,'Tag','listMag'),'Value');
str = get(findobj(fH,'Tag','listMag'),'String');
mag = str2num(str(indx,:));

% get lab from slider
lab = (101^get(findobj(fH,'Tag','sliderLab'),'Value')-1)/100;


nEl = round(polymSze / monomSze(1));   
nPix = round(pixSze / (mag * monomSze(1)));  % how many full helical rings do fall
                                             % inside a pixel

% get gain from slider
gain = 51^get(findobj(fH,'Tag','sliderGain'),'Value') - 1;


% get dark and shot noise values
if( get(findobj(fH,'Tag','checkDN'),'Value'))
   noise(1) = fwc*(51^get(findobj(fH,'Tag','sliderDN'),'Value')-1)/100;
else
   noise(1) = 0;
end;

noise(2) = get(findobj(fH,'Tag','checkSN'),'Value');
if noise(2)
   str = sprintf('shot noise: %5.2f %%',sn*100);
   set(findobj(fH,'Tag','textSN'),'String',str);
end;
  

function wait(speed)

t0 = clock;
j = 0;
while(etime(clock,t0)<speed)
   pause(0.2);
   j = j + 1;
end;

%------------------------------------------------------------------------
% Wrapper to various measures for contrast 
function ret = getContrastMeas(funName,im);
cmdStr = [funName,'(im)'];
ret = eval(cmdStr);

function ret=meanModulation(im);
% simple mean modulation

indxLmax = locmax1d(im);
if(isempty(indxLmax))
   mlmax = 0;
else
   mlmax = mean(im(indxLmax));
end;

indxLmin = locmin1d(im);
if(isempty(indxLmin))
   mlmin = 0;
else
   mlmin = mean(im(indxLmin));
end;

if( (mlmax+mlmin) > 0 ) 
   ret = (mlmax - mlmin) / (mlmax+mlmin);
else
   ret = 0;
end;


function ret=peakModulation(im);
% simple mean modulation

% just for a try: what is the mean density of local max
%global FIRSTTIME__
%persistent MEANLMAX__ NCALLS__

%if FIRSTTIME__ == 1
%    MEANLMAX__ = 0;
%    NCALLS__ = 0;
%    FIRSTTIME__ = 0;
%end


indxLmax = locmax1d(im);
if(isempty(indxLmax))
   ret = 0;
   return;
end;

indxLmin = locmin1d(im);
if(isempty(indxLmin))
   ret = 0;
   return;
end;

if(indxLmin(1)>indxLmax(1))
    % series starts with a loc max -> chop it
    indxLmax = indxLmax(2:end);
    if isempty(indxLmax)
       ret = 0;
       return;
    end;
end;
if(indxLmax(end)>indxLmin(end))
    % series ends with a loc max -> chop it
    indxLmax = indxLmax(1:end-1);
    if isempty(indxLmax)
       ret = 0;
       return;
    end;
end;

% there are a few pathalogic situations, where no local minima are found
% this happens for very low labeling ratios with no noise added -> the signal drops to 
% the background level without a negitive peak detected by locmin1d();
% then one has to find the positions in the list where two locmax succeed each other,
% and insert a locmin for this case
if(length(indxLmax)>=length(indxLmin))
   if length(indxLmax) == 2
      % super pathalogical case; the series consists of only two minima at the 
      % border and 2 or more local maxima in between; return 0 for this situation
      ret = 0;
      return;
   else 
      indxLminNew = [];
      iIndxLmax = 1;
      for i = 1:(length(indxLmin)-1)
         if(indxLmin(i)<indxLmax(iIndxLmax))
            indxLminNew = [ indxLminNew , indxLmin(i) ];
            iIndxLmax = iIndxLmax + 1;
         else
            indxLminNew = [ indxLminNew , round((indxLmax(iIndxLmax)+indxLmax(iIndxLmax+1))/2)];
            iIndxLmax = iIndxLmax + 2;
         end
         indxLminNew = [ indxLminNew , indxLmin(end) ];
      end;
      indxLmin = indxLminNew;   
   end;
end;

m = [];
iIndxLmin = 1;
for i = 1:length(indxLmax)
   bg = 1/2*(im(indxLmin(i))+im(indxLmin(i+1)));
   m(i) = (im(indxLmax(i)) - bg) / (im(indxLmax(i))+bg);;
end;

if ~isempty(m)
   ret = mean(m);
else
   ret = 0;
end;

% just for a try
%MEANLMAX__ = (NCALLS__ * MEANLMAX__ + length(indxLmax))/(NCALLS__ + 1);
%NCALLS__ = NCALLS__ + 1;
%MEANLMAX__
%length(indxLmax)