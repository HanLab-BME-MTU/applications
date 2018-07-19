function [ func ] = generateBifurcationFunction( response, Korg, freq )
%generateBifurcationFunction Generates a function that outputs the value of
%the first and second partial derivative of the orientation response
%function with respect to theta as well as the Jacobian

if(nargin < 3)
    freq = false;
end

if(freq)
    response_hat = response;
else
    response_hat = fft(response);
end

func = @bifurcationFunction;

    function [F,J] = bifurcationFunction(x)
        if(nargout > 1)
            derivs = 1:4;
        else
            derivs = 1:2;
        end
        D = 2*pi^2;
        theta = x(1,:);
        time = x(2,:);
        K = (1./sqrt(time)-1)/2;
        response_hat_at_K = orientationSpace.getResponseAtOrderVecHat(response_hat,Korg, K);
        values = interpft1_derivatives(response_hat_at_K,theta,derivs,2*pi,true);
        F = values(:,:,1:2);
        F = permute(F,[ 3 2 1]);
%         F(1) = values(1);
%         F(2) = values(2);
        if(nargout > 1)
            J = values(:,:,[2 3 3 4]);
            J(:,:,3:4) = D*J(:,:,3:4);
            J = reshape(J,[],2,2);
            J = permute(J,[2 3 1]);
            J = num2cell(J,[1 2]);
            J = blkdiag(J{:});
        end
%         J(1,1) = values(2);
%         J(1,2) = D*values(3);
%         J(2,1) = values(3);
%         J(2,2) = D*values(4);
    end
    


end

