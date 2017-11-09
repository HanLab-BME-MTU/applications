function plotWindowsColormapped(windows,matIn,cMap)

if nargin < 3 || isempty(cMap)
    cMap = jet;
end

    

nStrip = numel(windows);

minVal = min(matIn(:));
maxVal = max(matIn(:));
rng = maxVal - minVal;
nCol = size(cMap,1);

hold on

for j = 1:nStrip       
    if ~isempty(windows{j})
        nBand = numel(windows{j});            
        for k = 1:nBand
            if ~isempty(windows{j}{k})
                if ~isnan(matIn(j,k))
                    
                    %If they input an RGB specification instead of a scalar map
                    if size(matIn,3) == 3
                        winCol = matIn(j,k,:);
                    else
                        %Do lazy rounding for now. TEMP - switch to linear
                        %interpolation? Too slow?
                        iCol = max(ceil( (matIn(j,k) - minVal) / rng * (nCol-1)),1);                                            
                        winCol = cMap(iCol,:);                    
                    end
                    currWin = [windows{j}{k}{:}];
                    patch(currWin(1,:),currWin(2,:),winCol,'EdgeColor','none');
                else
                    %plotWindows(windows{j}{k},{'w','FaceAlpha',0});
                end                    
            end
        end    
    end
end

axis image

