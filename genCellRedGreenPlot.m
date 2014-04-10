function genCellRedGreenPlot( imageData, cellStats, displayrange )

    hGRPlot = figure;
    hRatioPlot = figure;
    
    ratioG1 = [];
    ratioS = [];
    ratioG2M = [];
    
        for i = 1:numel(cellStats)
            for j = 1:3
                cellChannelIntensities = mat2gray( imageData{j}( cellStats(i).PixelIdxList ), displayrange(j,:) );
                cellStats(i).meanChannelIntensity{j} = mean( cellChannelIntensities );
                cellStats(i).maxChannelIntensity{j} = max( cellChannelIntensities );
                cellStats(i).medianChannelIntensity{j} = median( cellChannelIntensities );
                cellStats(i).minChannelIntensity{j} = min( cellChannelIntensities );
                cellStats(i).medianGreenRedRatio = cellStats(i).medianChannelIntensity{2} / (cellStats(i).medianChannelIntensity{3} + 0.001);
                
                cellColor = [];
                switch cellStats(i).cellPatternType
                
                    case { 'Mono_G1' }
                        
                        cellColor = [1, 0, 0];
                        ratioG1(end+1) = log( cellStats(i).medianGreenRedRatio );
                        
                    case { 'Mono_S' }
                        
                        cellColor = [0, 1, 1];
                        ratioS(end+1) = log( cellStats(i).medianGreenRedRatio );
                        
                    case { 'Mono_G2', 'Mono_Prophase' }
                        
                        cellColor = [0, 1, 0];
                        ratioG2M(end+1) = log( cellStats(i).medianGreenRedRatio );
                        
                end
                
                if ~isempty( cellColor )
                    
                    figure( hGRPlot );
                    plotstat = cellStats(i).meanChannelIntensity;     
                    hold on;
                        plot( plotstat{2}, plotstat{3}, 'x', 'Color', cellColor, 'LineWidth', 2.0 );
                    hold off;

                    figure( hRatioPlot );
                    hold on;
                        plot( log( cellStats(i).medianGreenRedRatio ), log( cellStats(i).medianGreenRedRatio ), 'x', 'Color', cellColor, 'LineWidth', 2.0 );
                    hold off;

                end
                
            end
        end    
    
        figure, nhist( {ratioG1, ratioS, ratioG2M}, 'legend', { 'G1', 'S', 'G2-M'}, 'pdf', 'smooth' );
        
end