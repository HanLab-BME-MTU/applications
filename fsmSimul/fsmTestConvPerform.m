function estim = fsmTestConvPerform(shifts,nExp)
%TESTCONVPERFORM tests the performance of the convolution in estimating speckle movement
%
% SYNOPSIS estim = fsmTestConvPerform(shifts,nExp)
% 
% INPUT    shifts : can be a scalar (then the full signal shifts [pix])
%                   or a vector, where the first number must be negative;
%                   then the left part of the speckle chain is moving 
%                   to the left, and the right part to the right. The 
%                   part in between is replaced with new speckles
%          nExp   : number of experiments
%
% OUTPUT   estim :  estimated shifts
%
% NOTE     the function also displays a graphics window with the shifts
%          plots

% some global constants for the graphics
speed=0.5;

% generate a first speckle signal
ctrlData = fsmGetDflts;

% how many monomers do fall inside a pixel
nMonoPerPix = round(ctrlData.pixSize / (ctrlData.mag * ctrlData.monomSize));
    
% get the pixelated picture, the new polymer (where nMonoPerPix have been added)
% and the modulation value 
[cm,refPolym,refPixPicPolym] = fsmKernel(ctrlData,[],0);

% make sure that the number of pixel is odd; if not chop one pixel off the end
if(~rem(length(refPixPicPolym),2))
    refPixPicPolym = refPixPicPolym(1:end-1);
end;

% prepare image graphics
graphWinH = figure;
pos = get(graphWinH,'Position');
set(graphWinH,'Position',[20,20,pos(3),pos(4)]);
lineH = plot(refPixPicPolym);
set(lineH,'EraseMode','background');

wait(speed);

for i = 1 : nExp
    switch length(shifts)
    case 1, 
        nRepl = round(nMonoPerPix*shifts);
        if(shifts > 0)
            polym = [(rand(nRepl,1)<ctrlData.LR);refPolym(1:end-nRepl)];
        else
            polym = [refPolym(nRepl+1:end);(rand(nRepl,1)<ctrlData.LR)];
        end;
    case 2,
        if(shifts(1) > 0)
            shifts(1) = -1*shifts(1)
        end;
        if(shifts(2) < 0)
            shifts(2) = -1*shifts(2)
        end;
        nReplL = round(nMonoPerPix*abs(shifts(1)));
        nReplR = round(nMonoPerPix*abs(shifts(2)));
        refPolymL = refPolym(1:floor(length(refPolym)/2));
        refPolymR = refPolym(floor(length(refPolym)/2)+1:end);
        polym = [refPolymL(nReplL+1:end);(rand(nReplL,1)<ctrlData.LR);...
            (rand(nReplR,1)<ctrlData.LR);refPolymR(1:end-nReplR)];
    end;
    
    % generate the pixelated image
    [cm,polym,pixPicPolym] = fsmKernel(ctrlData,polym,0);
    
    % make sure that the number of pixel is odd; if not chop one pixel off the end
    if(~rem(length(pixPicPolym),2))
        pixPicPolym = pixPicPolym(1:end-1);
    end;
    
    % plot the new pix image
    set(lineH,'XData',1:length(pixPicPolym),'YData',pixPicPolym);
    yLim = get(get(graphWinH,'CurrentAxes'),'YLim');
    if(max(pixPicPolym(:))>yLim)
        set(get(graphWinH,'CurrentAxes'),'YLim',[0,max(pixPicPolym(:))]);
    end;
    
    % calculate the shift response according to the convolution theorem
    convOut=real(fftshift(ifft(fft(refPixPicPolym).*conj(fft(pixPicPolym)))));
    
    shiftVals = -floor(length(refPixPicPolym)/2):floor(length(refPixPicPolym)/2);
    if i == 1
        convWinH = figure;
        pos = get(convWinH,'Position');
        set(convWinH,'Position',[40+pos(3),20,pos(3),pos(4)]);
        convH = plot(shiftVals,convOut,'r-');
        set(convH,'EraseMode','background');
        set(get(convH,'Parent'),'YLim',[0,max(convOut(:))]);
    else
        set(convH,'XData',shiftVals,'YData',convOut); 
        yLim = get(get(convWinH,'CurrentAxes'),'YLim');
        if(max(convOut(:))>yLim(2))
            set(get(convWinH,'CurrentAxes'),'YLim',[0,max(convOut(:))]);
        end;
    end;
    
    rankedPeaks=sort(convOut(locmax1d(convOut)));
    
    switch length(shifts)
    case 1, 
        estim(i) = -shiftVals(find(convOut == rankedPeaks(end)));
    case 2,
        estim(i,1) = -shiftVals(find(convOut == rankedPeaks(end)));
        estim(i,2) = -shiftVals(find(convOut == rankedPeaks(end-1)));
        if(estim(i,1) > estim(i,2))
            val = estim(i,1);
            estim(i,1) = estim(i,2);
            estim(i,2) = val;
        end;
    end;
   
    [cm,refPolym,refPixPicPolym] = fsmKernel(ctrlData,polym,0);
    if(~rem(length(refPixPicPolym),2))
        refPixPicPolym = refPixPicPolym(1:end-1);
    end;
    
    wait(speed);
end;


%-------------------------------------------------------------------
function wait(speed)

t0 = clock;
j = 0;
while(etime(clock,t0)<speed)
   pause(0.1);
   j = j + 1;
end;

%-------------------------------------------------------------------
