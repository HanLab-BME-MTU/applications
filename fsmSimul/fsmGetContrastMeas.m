function ret = fsmGetContrastMeas(funName,im);
%FSMGETCONTRASTMEAS generic routine to call a set of contrast measures
% 
% SYNOPSIS fsmGetContrastMeas(funName,im)
% 
% INPUT    funName : function name 
%                    currently supported functions: 'meanModulation' 
cmdStr = [funName,'(im)'];
ret = eval(cmdStr);


%----------------------------------------------------------------------------
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
%----------------------------------------------------------------------------