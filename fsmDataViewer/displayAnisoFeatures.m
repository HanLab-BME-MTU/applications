function displayAnisoFeatures(hAxes, tag, X, layerColor)

line(X(:,[1 3])', X(:,[2 4])','LineWidth', 2, 'Color',layerColor,'Tag',tag,'Parent',hAxes);
