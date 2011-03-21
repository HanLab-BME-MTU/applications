function [roughness,freqSpec]=calcEdgeRoughness(sheetEdge,r,toDoList)
% first order the curve such that it starts at the bounderies:
for frame=toDoList
    if length(sheetEdge(frame).pos(:,1))==1
        % This happens only if we have single cell data:
        roughness=NaN*zeros(1,toDoList(end));
        
        freqSpec.ssA(1,frame)=NaN;
        freqSpec.f=NaN;
        
        goOn=0;        
    end
end
if ~goOn
    display('This is single cell data, don''t calculate roguhness and freqSpec!')
    return;
end
for frame=toDoList
    % find jumps larger than one pixel:
    rawCurve=sheetEdge(frame).pos;
    steps=sqrt(sum((rawCurve(2:end,:)-rawCurve(1:end-1,:)).^2,2));
    jump=find(steps>1+eps);
    if length(jump)>1
        error('Something went wrong, found two jumps');
    elseif length(jump)==1
        % Find out which y-value is larger:
        if rawCurve(jump,1)<rawCurve(jump+1,1)
            % then it is decending
            ordCurve=vertcat(rawCurve(jump+1:end,:),rawCurve(1:jump,:));
            ordCurve=flipud(ordCurve);
        else
            % it is ascending:
            display('The ascending case should be checked by hand');
            break;
            % ordCurve=flipud(rawCurve);
            ordCurve=vertcat(flipud(rawCurve(1:jump,:)),flipud(rawCurve(jump+1:end,:)));
        end
    elseif isempty(jump)
        % the curve is ordered:
        ordCurve=rawCurve;
    end
    stepsNew=sqrt(sum((ordCurve(2:end,:)-ordCurve(1:end-1,:)).^2,2));
    jump=find(stepsNew>1+eps);
    if length(jump)>1
        error('Something went wrong, found another jump');
    end
    % plot(ordCurve(:,1),ordCurve(:,2))
    
    % smooth the curve:
    smoothCurve=ordCurve(1:r:end,:);
    % add the last point:
    if ~compPts(smoothCurve(end,:),ordCurve(end,:))
        smoothCurve(end+1,:)=ordCurve(end,:);
    end
    
    % measure the curve length:
    L0=sqrt(sum((smoothCurve(end,:)-smoothCurve(1,:)).^2,2));
    
    dL=sqrt(sum((smoothCurve(2:end,:)-smoothCurve(1:end-1,:)).^2,2));    
    Ltot=sum(dL);
    
    roughness(frame)=Ltot/L0;
    
    curve(frame).pos=ordCurve;
end
% figure(1),plot(roughness);

% figure(2)
% calculate the Fourier spectrum:
for frame=toDoList
    % Determine the shift in x
    xshift=mean([curve(frame).pos(1,2) curve(frame).pos(end,2)]);
    
    shftCurve     =curve(frame).pos;
    shftCurve(:,2)=curve(frame).pos(:,2)-xshift;
    
    % Order all points in the curve according to their y-coordinate:
    sortCurve=sortrows(shftCurve);
    
    % Find all points with the same y-coordinate:
    [~,bin] = histc(sortCurve(:,1),(0:max(sortCurve(:,1)))+0.5);
    
    % lowerCurve=sortCurve;
    for binId=1:max(bin)
        % Take always the lower value:
        lv=min(sortCurve(binId==bin,:),[],1);
        if ~isempty(lv)
            lowerCurve(binId,:)=lv;
        else
            lowerCurve(binId,:)=[NaN NaN];
        end
        if ~isnan(lowerCurve(binId,1)) && lowerCurve(binId,1)~=binId
            error('Values are not correctly sorted!')
        end
    end
    % Remove all NaNs:
    checkVec  =~isnan(lowerCurve(:,2));    
    cleanCurve=lowerCurve(checkVec,:);

    % Do also horizontal average:
    h = fspecial('average', [r,1]);
    
    avCurve=cleanCurve;
    avCurve(:,2)=imfilter(cleanCurve(:,2), h, 'replicate');
    
%     plot(shftCurve(:,1),shftCurve(:,2))
%     hold on;
%     plot(avCurve(:,1), avCurve(:,2),'k')
%     hold off;    
    
    Fsmp   = 1;                             % Sampling frequency [1/pix]
    numPts = size(avCurve,1);               % Length of signal
    y = avCurve(:,2);                       % The signal
    
    NFFT = 2^nextpow2(numPts);              % Next power of 2 from length of y
    Y    = fft(y,NFFT)/numPts;              % Calculate the FFT
    f    = (Fsmp/2*linspace(0,1,NFFT/2+1))';   % The spatial frequencies

% Plot single-sided amplitude spectrum.
%     plot(f,2*abs(Y(1:NFFT/2+1))) 
%     title('Single-Sided Amplitude Spectrum of y(x)')
%     xlabel('Frequency (1/pix)')
%     ylabel('|Y(f)|')
%     xlim([0,0.1])
%     ylim([-130,130])
%     pause(0.01)    
    
    % Store the single sided amplitude:
    freqSpec.ssA(:,frame)=2*abs(Y(1:NFFT/2+1));
    % Check if the frequency is really the same as before
    if isfield(freqSpec,'f')
        if sum(freqSpec.f-f)<eps
            freqSpec.f=f;
        else
            break;
        end
    else
        freqSpec.f=f;
    end
    
    clear lowerCurve
end

% figure(3)
% aveWindow=24;
% % show moving average:
% for frame=1:toDoList(end)-aveWindow
%     aveFreqSpeq(frame).Y    = mean(horzcat(freqSpec(frame:frame+aveWindow).Y),2);
%     aveFreqSpeq(frame).Ystd =  std(horzcat(freqSpec(frame:frame+aveWindow).Y),[],2);
%     aveFreqSpeq(frame).f    = mean(horzcat(freqSpec(frame:frame+aveWindow).f),2);
%     aveFreqSpeq(frame).fstd =  std(horzcat(freqSpec(frame:frame+aveWindow).f),[],2);
%     
%     plot(aveFreqSpeq(frame).f,2*abs(aveFreqSpeq(frame).Y(1:NFFT/2+1))) 
%     title('Single-Sided Amplitude Spectrum of y(x)')
%     xlabel('Frequency (1/pix)')
%     ylabel('|Y(f)|')
%     xlim([0,0.02])
%     ylim([0,130])
%     pause(0.01)
% end
% 
% frame=1;
% plot(aveFreqSpeq(frame).f,2*abs(aveFreqSpeq(frame).Y(1:NFFT/2+1)),'b')
% hold on
% frame=length(aveFreqSpeq);
% plot(aveFreqSpeq(frame).f,2*abs(aveFreqSpeq(frame).Y(1:NFFT/2+1)),'r')
% hold off
% title('Single-Sided Amplitude Spectrum of y(x)')
% xlabel('Frequency (1/pix)')
% ylabel('|Y(f)|')
% xlim([0,0.02])
% ylim([0,130])