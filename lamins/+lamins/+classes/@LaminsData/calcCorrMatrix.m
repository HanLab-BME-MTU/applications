function corrmatrix = calcCorrMatrix(obj,cache)
    if(nargin < 2)
        cache = true;
    end
    filename = [ obj.outputDirectory_ filesep 'corrmatrix.mat'];
    if(cache && exist(filename,'file'))
        C = load(filename);
        corrmatrix = C.corrmatrix;
    else
        % not cached
        corrmatrix = calcCorrMatrix(obj.movieData);
        save(filename,'corrmatrix');
    end
end

