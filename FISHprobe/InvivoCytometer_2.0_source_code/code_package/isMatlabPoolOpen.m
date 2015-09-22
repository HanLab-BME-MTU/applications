function flagPoolOpenedAlready = isMatlabPoolOpen()

    if exist('gcp')
        flagPoolOpenedAlready = ~isempty(gcp('nocreate'));        
    else
         flagPoolOpenedAlready = matlabpool('size') > 0;
    end

end