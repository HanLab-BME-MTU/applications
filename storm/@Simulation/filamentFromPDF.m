function t = filamentFromPDF(obj,controlPoints,sigmaNoise,pdf)

len = lengthBezier(controlPoints);
nT = 2*round(len/sum(pdf(:,1).*pdf(:,2)));

val = valueFromPDF(pdf,nT);
pos = [0;cumsum(val)];
pos = pos(pos<len);
t = [pos;len]/len;
t = arcLengthToNativeBezierParametrization(controlPoints,t);

% Generate normalized 3D noise
nSamples = numel(t);
noiseGaussianNorm = randn(nSamples,3);

% Compute the tangent with unitary length in the normalized data set
cPNorm = controlPoints./repmat(sigmaNoise,size(controlPoints,1),1);
[~,tangent] = tangentBezier(cPNorm,t);

% Projection the normalized 3D noise onto the orthogonal plane (tangent is normalized)
% noiseGaussianNormOrtho = tangent x (noiseGaussianNorm x tangent/norm(tangent))/norm(tangent)
% noiseGaussianNormOrtho = tangent x (noiseGaussianNorm x tangent)
noiseOrthoNorm = cross(tangent,cross(noiseGaussianNorm,tangent,2),2);

% Rescale the noise
noiseOrtho = noiseOrthoNorm.*repmat(sigmaNoise,nSamples,1);
points = renderBezier(controlPoints,t);

% % Add the noise to the points
points = points + noiseOrtho;
obj.data.points = [obj.data.points;points];

% figure(44)
% bar(pdf(:,1),pdf(:,2));

% figure(55)
% counts = hist(val,pdf(:,1));
% bar(pdf(:,1),counts/sum(counts));

    function val = valueFromPDF(pdf,n)
        pdf = pdf(pdf(:,2)~=0,:); % Remove empty bins
        cdf = cumsum(pdf(:,2)); % Compute the cdf
        
        binSize = pdf(2,1)-pdf(1,1);
        
        if cdf(1) ~= 0
            cdf = [0;cdf];
            pdf = [[pdf(1,1)-binSize/2,0];pdf];
        end
        if cdf(end) ~= 1
            cdf = [cdf;1];
            pdf = [pdf;[pdf(end,1)+binSize/2,0]];
        end
        
        prob = rand(n,1);
        val = interp1(cdf,pdf(:,1),prob,'cubic');
    end
end
