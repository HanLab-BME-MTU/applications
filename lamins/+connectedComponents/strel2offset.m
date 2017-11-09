function offset = strel2offset(siz,se)
    offset = connectedComponents.mat2offset(siz,se.getnhood);
end
