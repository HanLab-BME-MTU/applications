function fprintfc( fidVec, varargin )

    for i = 1:numel(fidVec)
        fprintf(fidVec(i), varargin{:} ); 
    end
    
end