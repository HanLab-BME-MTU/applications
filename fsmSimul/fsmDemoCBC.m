function fsmDemoCBC(tp)
% callback functions for the fsmDemo 

if nargin == 0
   tp = 'callback';
end;

switch get(gcbo,'Tag')
case 'sliderPixSze',
   if(strcmp(tp,'create'))
      sliderPixSzeCreateCBC;
   else
      airyUnitChg;
   end;
case 'sliderGain',
   if(strcmp(tp,'create'))
      sliderGainCreateCBC;
   else
      sliderGainCBC;
   end;
case 'textGain',
   textGainCreateCBC;
case 'textPixSze',
   textPixSzeCreateCBC;
case 'listMag',
   if(strcmp(tp,'create'))
      listMagCreateCBC;
   else
      airyUnitChg;
   end;
case 'listNA',
   if(strcmp(tp,'create'))
      listNACreateCBC;
   else
      airyUnitChg;
   end;
case 'sliderLab',
   if(strcmp(tp,'create'))
      sliderLabCreateCBC;
   else
      sliderLabCBC;
   end;
case 'textLab',
   textLabCreateCBC;  
case 'checkDN',
   if(strcmp(tp,'create'))
      checkDNCreateCBC;
   else
      checkDNCBC;
   end;
case 'sliderDN',
   if(strcmp(tp,'create'))
      sliderDNCreateCBC;
   else
      sliderDNCBC;
   end;
case 'textDN',
   textDNCreateCBC;
case 'checkSN',
   if(strcmp(tp,'create'))
      checkSNCreateCBC;
   else
      checkSNCBC;
   end;
case 'textSN',
   textSNCreateCBC;   
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function airyUnitChg

global WVL__;

pixSze= get(findobj(gcbf,'Tag','sliderPixSze'),'Value');

indx = get(findobj(gcbf,'Tag','listMag'),'Value');
str = get(findobj(gcbf,'Tag','listMag'),'String');
mag = str2num(str(indx,:));

indx = get(findobj(gcbf,'Tag','listNA'),'Value');
str = get(findobj(gcbf,'Tag','listNA'),'String');
NA = str2num(str(indx,:));

airyUnits = pixSze * 1000 / mag / (1.22 * WVL__/(2*NA));

tH = findobj(gcbf,'Tag','textPixSze');
str = sprintf('pix size %4.1f um\n          %3.1f au',pixSze,airyUnits);
set(tH,'String',str);

function sliderPixSzeCreateCBC;
set(gcbo,'Value',6.4);

function textPixSzeCreateCBC
global WVL__;

airyUnits = 6.4 * 1000 / 100 / (1.22 * WVL__/(2*1.4));
str = sprintf('pix size %4.1f um\n          %3.1f au',6.4,airyUnits);
set(gcbo,'String',str);


function listMagCreateCBC
str(1,:) = ' 10';
str(2,:) = ' 20';
str(3,:) = ' 25';
str(4,:) = ' 40';
str(5,:) = ' 50';
str(6,:) = ' 63';
str(7,:) = '100';

set(gcbo,'String',str);
set(gcbo,'Value',7);

function listNACreateCBC

str(1,:) = '0.1';
str(2,:) = '0.2';
str(3,:) = '0.3';
str(4,:) = '0.4';
str(5,:) = '0.5';
str(6,:) = '0.6';
str(7,:) = '0.7';
str(8,:) = '0.8';
str(9,:) = '0.9';
str(10,:) = '1.2';
str(11,:) = '1.4';

set(gcbo,'String',str);
set(gcbo,'Value',11);


function textLabCreateCBC
str = sprintf('%6.2f %% labelled',0);
set(gcbo,'String',str);

function sliderLabCreateCBC
set(gcbo,'Value',0);  

function sliderLabCBC
tH = findobj(gcbf,'Tag','textLab');
val = 101^(get(gcbo,'Value')) - 1;   % exponetial slider
str = sprintf('%6.2f %% labelled',val);
set(tH,'String',str);

function sliderGainCreateCBC
set(gcbo,'Value',log(2)/log(51));   
% prepare for logarithmic gain slider with 
% values [0,10]

function sliderGainCBC
tH = findobj(gcbf,'Tag','textGain');
val = 51^(get(gcbo,'Value')) - 1;   % exponetial slider
str = sprintf('gain: %4.1f',val);
set(tH,'String',str);

function textGainCreateCBC
str = sprintf('gain: %4.1f',1);
set(gcbo,'String',str);

function checkDNCreateCBC
set(gcbo,'Value',0);

function checkDNCBC
if( get(gcbo,'Value') == 1)
   tH = findobj(gcbf,'Tag','textDN');
   set(tH,'Enable','on');
   sH = findobj(gcbf,'Tag','sliderDN');
   set(sH,'Enable','on');
else
   tH = findobj(gcbf,'Tag','textDN');
   set(tH,'Enable','off');
   sH = findobj(gcbf,'Tag','sliderDN');
   set(sH,'Enable','off');
end;

function sliderDNCreateCBC
set(gcbo,'Value',log(11)/log(51)); 
% set the initial dark noise value to 10%
set(gcbo,'Enable','off');

function sliderDNCBC
tH = findobj(gcbf,'Tag','textDN');
val = 51^(get(gcbo,'Value')) - 1;   % exponetial slider
str = sprintf('dark noise: %5.2f %%',val);
set(tH,'String',str);

function textDNCreateCBC
str = sprintf('dark noise: %5.2f %%',10);
set(gcbo,'String',str);
set(gcbo,'Enable','off');

function checkSNCreateCBC
set(gcbo,'Value',0);

function checkSNCBC
if( get(gcbo,'Value') == 1)
   tH = findobj(gcbf,'Tag','textSN');
   set(tH,'Enable','on');
else
   tH = findobj(gcbf,'Tag','textSN');
   set(tH,'Enable','off');
end;

function textSNCreateCBC
str = sprintf('shot noise: %5.2f %%',0);
set(gcbo,'String',str);
set(gcbo,'Enable','off');
