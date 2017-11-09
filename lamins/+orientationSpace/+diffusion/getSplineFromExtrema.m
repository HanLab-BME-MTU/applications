function sp = getSplineFromExtrema(extrema,response,K,nDerivs)
    [extrema,events] = orientationSpace.diffusion.alignExtrema(extrema);
    % Unwrap extrema from the periodic boundary to avoid discontinuities
    extrema = orientationSpace.diffusion.unwrapExtrema(extrema,events);
    extrema = fliplr(extrema);
    K = fliplr(K);
    response = fliplr(response);
    validTracks = sum(~isnan(extrema),2) > 1;
    extrema = extrema(validTracks,:);
    derivs = orientationSpace.diffusion.orientationMaximaDerivatives(response,K,nDerivs,extrema);
    
    

    
    x = joinColumns(repmat(K,nDerivs+1,1));
    sp(size(extrema,1)) = struct('form',[],'knots',[],'coefs',[],'number',[],'order',[],'dim',[]);
    for trackNum = 1:size(extrema,1)
        y = joinColumns([extrema(trackNum,:); permute(derivs(trackNum,:,:),[3 2 1])]);
        valid = ~isnan(y);
%         try
%             knots = optknt(x(valid).',nDerivs+3);
%         catch err
%             knots = aptknt(x(valid).',nDerivs+3);
%         end
        knots = aptknt(x(valid).',nDerivs+3);
        sp(trackNum) = spapi(knots,x(valid).',y(valid).');
    end
end