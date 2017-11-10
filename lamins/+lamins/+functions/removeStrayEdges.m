function FE = removeStrayEdges(edges,FE)
    import lamins.functions.*;
    import connectedComponents.*;
    endpts = cellfun(@(e) getEdgeEndpoints(ccFilter(edges,e)),FE,'UniformOutput',false); 
    mult = cellfun(@getMultiplicityInt,endpts,'UniformOutput',false);
    
    %% deal with edges with a single endpoint
    hasStrayEdges = cellfun(@(x) any(x(:) == 1),mult);
%     hasStrayEdgesIdx = find(hasStrayEdges);
%     uniqueEndpts = cellfun(@unique,endpts(hasStrayEdges),'UniformOutput',false);
%     singleEndpts = cellfun(@(u,m) u(m == 1), uniqueEndpts, mult(hasStrayEdges) , 'UniformOutput',false);
%     edgeLogical = cellfun(@(e,s) any(ismember(e,s),2), endpts(hasStrayEdges), singleEndpts, 'UniformOutput',false);
%     FE(hasStrayEdges) = cellfun(@(e,l) e(~l), FE(hasStrayEdges), edgeLogical, 'UniformOutput',false);
    FE(hasStrayEdges) = cellfun(@removeEdgesWithSingleEndpoints,endpts(hasStrayEdges),mult(hasStrayEdges),FE(hasStrayEdges),'UniformOutput',false);
    
    %% deal with extra loops
    hasExtraEdges = cellfun(@(x) any(x(:) == 4),mult);
    FE(hasExtraEdges) = cellfun(@removeExtraEdges,endpts(hasExtraEdges),mult(hasExtraEdges),FE(hasExtraEdges),'UniformOutput',false);
end
function FE = removeEdgesWithSingleEndpoints(endpts,mult,FE)
    uniqueEndpts = unique(endpts);
    singleEndpts = uniqueEndpts(mult == 1);
    edgeLogical = any(ismember(endpts,singleEndpts),2);
    FE = FE(~edgeLogical);
end
function FE = removeExtraEdges(endpts,mult,FE)
    assigned = assignMultiplicity(endpts,mult);
    edgeLogical = all(assigned == 4,2);
    FE = FE(~edgeLogical);
end