function [ out ] = showLabelMatrix( lm , cmfunc)
%showLabelMatrix Show a label matrix with an appropritate colormap

if(nargin < 2)
    cmfunc = [];
end

out = showPropMatrix(lm,cmfunc);

end

