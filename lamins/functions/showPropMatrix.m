function [ out ] = showPropMatrix( matrix , cmfunc, shuffle, show)
%showPropMatrix shows the property matrix
% matrix is a matrix containing property or label values
% cmfunc is a colormap function
% shuffle is a boolean value whether to shuffle or not
% show is a boolean value whether to actually create the gui
%      This is mainly in the case that one just wants the results

    bg = [0 0 0];
    
    if(isstruct(matrix))
        matrix = labelmatrix(matrix);
    end

    % make matrix a double since we will use NaN as a flag
    matrix = double(matrix);

    isNanMatrix = any(isnan(matrix(:)));
    
    if(~isNanMatrix)
        % if no NaN, then assume 0 should be black
        matrix(matrix == 0) = NaN;
    end
    
    % prop should be the unique values in matrix
    prop = unique(matrix(~isnan(matrix)));
    uProp = prop;
    nUnique = length(uProp);
    
    if(nargin < 2 || isempty(cmfunc))
        cmfunc = @jet;
    end
    if(ischar(cmfunc))
        cmfunc = eval(['@' cmfunc]);
    end
    
    if(nargin < 3 || ismempty(shuffle))
        % if there are no NaN (label matrix), then shuffle
        % otherwise, this is a property matrix, so shuffle
        shuffle = ~isNanMatrix;
    end
    
    if(nargin < 4)
        show = true;
    end

    % bgValue should be the new minimum
    bgValue = min(uProp) - (max(uProp) - min(uProp))/(nUnique-1);
    matrix(isnan(matrix)) = bgValue;
    
    cm = cmfunc(nUnique);
    
    if(shuffle)
        cm = cm(randperm(nUnique),:);
    end
    
    cm = [ bg; cm];
    
%     figure;
    if(show)
        out.him = imshow(matrix,[]);
        colormap(cm);
        if(~shuffle)
            out.hcb = colorbar;
        end
    end
    
    out.matrix = matrix;
    out.cm = cm;

end

