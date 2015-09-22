function startMatlabPool()

    if exist('parpool')
        parpool;            
    else
        matlabpool open;
    end

end