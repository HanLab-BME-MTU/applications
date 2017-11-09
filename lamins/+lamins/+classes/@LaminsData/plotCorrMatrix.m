function plotCorrMatrix(obj)
    import lamins.functions.*;
    order = obj.params.channels.order;
    corrmatrix = obj.calcCorrMatrix();
    corrmatrix = corrmatrix(order,:,order,:);
    plotcorrmatrix(corrmatrix,obj.movieData.getFilename());
end

