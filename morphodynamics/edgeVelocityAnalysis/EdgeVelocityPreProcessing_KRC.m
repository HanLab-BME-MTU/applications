function ProcessedTimeSeries = EdgeVelocityPreProcessing_KRC(timeSeries, nIMFs)

% Number of Modes to Subtract -----------------------------------
nModes = nIMFs;
%----------------------------------------------------------------

timeSeries(:,1) = 0; % Fist Column is usually full of NaNs.
timeSeries(isnan(timeSeries)) = 0;

ProcessedTimeSeries = zeros(size(timeSeries));
nWin = size(timeSeries,1);
warning off

switch nModes
    case 0
        ProcessedTimeSeries = timeSeries; 
    case 1
        for w = 1:nWin
               TS = timeSeries(w,:);
               IMFs = emd(TS);
               if size(IMFs,1) >= 2
                   ProcessedTimeSeries(w,:) = TS - IMFs(1,:); % Subtract first IMF
               else
                   disp('Low Numer of IMFs')
                   ProcessedTimeSeries(w,:) = TS; % Subtract first IMF
               end
        end
    case 2
        for w = 1:nWin
               TS = timeSeries(w,:);
               IMFs = emd(TS);
               if isequal(size(IMFs,1), 0)
                   ProcessedTimeSeries(w,:) = TS;
                   disp('Low Numer of IMFs')
               elseif isequal(size(IMFs,1), 1)
                   ProcessedTimeSeries(w,:) = TS; 
                   disp('Low Numer of IMFs')
               elseif isequal(size(IMFs,1), 2) 
                   ProcessedTimeSeries(w,:) = TS - IMFs(1,:); % Subtract First and Second IMFs
                   disp('Low Numer of IMFs')
               else
                   ProcessedTimeSeries(w,:) = TS - IMFs(1,:) - IMFs(2,:); % Subtract First and Second IMFs
               end
        end
end


