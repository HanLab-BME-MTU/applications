function [ theta, residual ] = localExtrema( orientationMatrix, ind )
%orientationSpace.localExtrema finds the absolute maxima
    import orientationSpace.*;
    s = size(orientationMatrix);
    M = real(reshape(orientationMatrix,s(1)*s(2),s(3)));
    n = (s(3))/2;
    
    % new
    % estimate zeros in derivative
    deriv = orientationSpace.derivative(real(orientationMatrix),pi/s(3)/8);
    derivDiff = -diff(deriv(:,:,[1:end 1]),1,3);
    offset = deriv./derivDiff;
    zeroEstimate = bsxfun(@plus,offset,shiftdim(0:s(3)*8-1,-1));
    zeroEstimate(offset < 0 | offset > 1) = NaN;
    zeroEstimate(derivDiff < 0 ) = NaN;

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(2*n-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = 2*exp(-xx.^2/2);
    
    if(nargin < 2)
        [~,ind] = max(M,[],2);
    elseif(isscalar(ind))
        ind = ones(s(1)*s(2),1)*ind;
    end
    
    xind = x(ind)';

    M = M/A;
    nM = size(M,1);
    
    % Initialize Batch Structures for Optimization
    batchSize = s(1);
    nBatches  = nM/batchSize;
    Mcell     = cell(1,nBatches);
    xindc     = cell(1,nBatches);
    theta     = cell(1,nBatches);
    residual     = cell(1,nBatches);

    
    % group indices into batches
    for i=1:nBatches
        idx = (1:batchSize)+(i-1)*batchSize;
        Mcell{i} = M(idx,:);
        xindc{i} = xind(idx);
    end

    options = optimoptions('lsqnonlin','Jacobian','on','Display','off','TolFun',1e-9);
%     progressText(0,'Solving');
    
    parfor i=1:nBatches
        offset = xindc{i};
        f = @(k) orientationDerivWrap(k,x,Mcell{i},n);
        [theta{i},~,residual{i}] = lsqnonlin(f,offset,offset-0.6,offset+0.6,options);
        %         progressText(i/nM*1024);
    end

    % assemble image
    theta = horzcat(theta{:});
    % wrap around one more time, needed in the even case
%     theta = wraparoundN(theta,-n,n);
    
%     if(~isreal(orientationMatrix))
%         theta_i = orientationSpace.localExtrema(cat(3,imag(orientationMatrix),-imag(orientationMatrix)),ind);
%         theta = theta + 2j*theta_i;
%     end
%     
%     theta = real(theta)*pi/s(3) + 1j*imag(theta);
    
    if(nargout > 1)
        residual = horzcat(residual{:});
    end

end
function [F,J] = orientationDerivWrap(k,x,M,n)
    % k is the angle of the candidate maximum
    % x is the angle space
    % M represents the coefficients, the deconvolved response
    % n represents the wraparound limits
    t = wraparoundN(bsxfun(@minus,k,x),-n,n);
    et = exp(-t.^2/2);
    F = -sum(M.*(t.*et),2);
    if(nargout > 1)
        J = sum(M.*(t.^2.*et - et),2);
        % along the main diagonal is of interest
        J = spdiags(J(:),0,length(F),length(k));
    end
end

