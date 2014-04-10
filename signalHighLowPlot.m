function signalHighLowPlot( x, y, yHigh, yLow, varargin )
% signalHighLowPlot( x, y, yHigh, yLow, varargin )
%
%
% Example: 
% 
%     x = 2:100;
%     y = log2(x);
%     yLow = y - 1;
%     yHigh = y + 1;
%     
%     signalHighLowPlot( x, y, yHigh, yLow, 'lineColor', [0 0 0], 'fillColor', [0 1 0 0.2] );
%     
% Author: Deepak Roy Chittajallu
%

    p = inputParser;
    p.addParamValue( 'lineColor', zeros(1,3), @(x) ( isnumeric(x) && numel(x) == 3 ) );
    p.addParamValue( 'fillColor', [zeros(1,3), 0.2], @(x) ( isnumeric(x) && ismember(numel(x), [3, 4]) ) );
    p.parse( varargin{:} );
    parameters = p.Results;

    if numel(unique([numel(x), numel(y), numel(yHigh), numel(yLow)])) ~= 1
        error( 'ERROR: x, y, yHigh, yLow must all be of th esame size' ); 
    end
    
    if numel(parameters.fillColor) == 3
        parameters.fillColor(4) = 1.0;        
    end
    
    x = x(:);
    y = y(:);
    yHigh = yHigh(:);
    yLow = yLow(:);
    
    p = patch( [x ; x(end:-1:1)], [yLow; yHigh(end:-1:1)], parameters.fillColor(1:3), 'EdgeColor', 'None' );
    set( p, 'FaceAlpha', parameters.fillColor(4) );
    hold on;
    plot(x, y, 'Color', parameters.lineColor, 'LineWidth', 2.0 );
    hold off;

end