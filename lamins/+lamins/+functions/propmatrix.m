function [ matrix ] = propmatrix( cc, prop, show )
%propmatrix Outputs an image matrix with connected objects labeled with the
%property
%
% cc should come from bwconncomp, or a label matrix
% prop is a vector of values with the property for each corresponding label
%      prop will usually come from regionprop
%      e.g. rp = regionprop(cc,'Area'); prop = [rp.Area];
%
% output
% matrix is the same size as the original image with pixels now having the
%        value of the object property

if(nargin < 3)
    show = false;
end

if(isstruct(cc))
% lm is a label matrix the size of the original image
    lm = labelmatrix(cc);
else
    % otherwise assume that cc is a label matrix
    lm = cc;
end

% idx is a vector of linear indices of the N non-zero pixels
idx = find(lm);
% labels is a vector of size N containing the label numbers
labels = lm(idx);

% prop is a mapping between label number and the property
% so replace

% use lm as a template for matrix
matrix = lm;

% TODO: cast matrix as whatever class prop is
% or maybe not depending on how we want to flag the non-object areas
matrix = double(matrix);
% NaN is a double value, this might be a hack
matrix(~matrix) = NaN;

matrix(idx) = prop(labels);

if(show)
    % if show is a string convert to a function handle
    if(ischar(show))
        show = eval(['@' show]);
    end
    if(~isa(show,'function_handle'))
        show = @jet;
    end
    figure;
    lamins.functions.showPropMatrix(matrix);
    
%     u = unique(prop);
%     nUnique = length(u);
%     blackValue = min(u) - (max(u) - min(u))/nUnique;
%     
%     showMatrix = matrix;
%     matrix(isnan(matrix)) = blackValue;
%     
%     cm = [ [0 0 0]; show(nUnique)];
%     
%     figure;
%     imshow(matrix,[]);
%     colormap(cm);
%     colorbar;
end

end

