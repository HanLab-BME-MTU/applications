function [compMat,compDate] = checkMovieArrayCompletion(MA,verbose)


processNames = {'SegmentationProcess3D',...             % 1
                'MaskGeometry3DProcess',...             % 2
                'SkeletonizationProcess',...            % 3
                'SkeletonPruningProcess',...            % 4
                'MaskObjectTrackingProcess'};           % 5


iChan =1 ;
nMov = numel(MA);
nProc = numel(processNames);

compMat = false(nMov,nProc);
compDate = zeros(nMov,nProc);


timeNow = now;

for j = 1:numel(MA)
    
    
    for k = 1:numel(processNames)
        
        tmp = MA(j).getProcessIndex(processNames{k},Inf,0);
        
        if numel(tmp) > 1
            warning(['Movie ' num2str(j) ' has ' num2str(numel(tmp)) ' ' processNames{k}])
        end
        
        if ~isempty(tmp)
            iProc = max(tmp);
            
            compMat(j,k) = MA(j).processes_{iProc}.checkChannelOutput(iChan);
            
            if compMat(j,k)                
                compDate(j,k) = datenum(MA(j).processes_{iProc}.dateTime_);
            end
            
        end
        
    end
    
end











