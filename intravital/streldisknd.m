function [diskStrel] = streldisknd( diskRadius )

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired( 'diskRadius', @(x) (isnumeric(x) && numel(x) > 1 && max(size(x)) == numel(x) && ~any(x < 0) && ~any(x-floor(x)>0)) );
    p.parse( diskRadius );
    
    indexRange = cell(1,numel(diskRadius));
    for i = 1:numel(diskRadius)
        indexRange{i} = 1:(2*diskRadius(i)+1);
    end
    
    subInd = cell(1,numel(diskRadius));
    [subInd{:}] = ndgrid( indexRange{:} );
    
    diskStrel = zeros( 2*diskRadius + 1 );
    for i = 1:numel(diskRadius)
        diskStrel = diskStrel + (subInd{i} - diskRadius(i) - 1).^2 / (eps + (diskRadius(i))^2);
    end
    diskStrel = double( (diskStrel - 1) <= 0 );
    
end