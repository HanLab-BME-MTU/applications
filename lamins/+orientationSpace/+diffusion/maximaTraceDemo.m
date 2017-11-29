function [q] = maximaTraceDemo()

gcp;

K = 8:-0.5:3;
threshold = pi/60;

I = orientationSpace.paper.loadLaminDemoImage;
F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
R = F*I;
maxima = R.getRidgeOrientationLocalMaxima;
maximaTrace = R.traceAllLocalMaxima(K,maxima);
parfor ki=1:length(K)
    res = real(shiftdim(R.getResponseAtOrderFT(K(ki),2).a,2));
    lm = shiftdim(maximaTrace(:,:,:,ki),2)*2;
    lmd = orientationSpace.diffusion.orientationMaximaDerivatives(res,K(ki),1,lm);
    maximaTraceDeriv(:,:,:,ki) = shiftdim(lmd,1);
end
[~,ind] = max(cumsum(abs(maximaTraceDeriv) < threshold,4),[],4);
[q.X,q.Y,q.Z] = meshgrid(1:1024,1:1024,1:7);
q.maximaTraceSelect = maximaTrace(sub2ind(size(maximaTrace),q.Y,q.X,q.Z,ind));
R3 = R.getResponseAtOrderFT(3,2);
q.nlms = R3.nonLocalMaximaSuppressionPrecise(q.maximaTraceSelect);
q.maximaTrace = maximaTrace;
q.maximaTraceDeriv = maximaTraceDeriv;

end