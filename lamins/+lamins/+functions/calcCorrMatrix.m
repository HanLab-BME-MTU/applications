function m = calcCorrMatrix(movieData)
% calcCorrMatrix Calculate correlation matrix which contains the values of
% corr2 between each image plane in MovieData
    if(isa(movieData,'MovieData'))
        % if we received a movieData object extract the data
        movieData = movieData.getReader();
    end
    
    if(isa(movieData,'Reader'))
        % if it is a reader, use CellReader to emulate the cell interface
        reader = CellReader ( CachedReader( movieData ) );
        nc = reader.getSizeC();
        nz = reader.getSizeZ();
%         J = reader(:);
    else
        % otherwise assume movieData is what is needed
        [nc, nz] = size(movieData);
        J = movieData(:);
    end
    % J should be a 1D cell array
    N = prod(size(J));
    m = zeros(N);
    for i=1:N
        m(i,i) = 1;
        M = J{i};
        for j=i+1:N
            m(i,j) = corr2(M,J{j});
        end
    end
    mt = m';
    m(~m) = mt(~m);
    m = reshape(m,[nc nz nc nz]);
end