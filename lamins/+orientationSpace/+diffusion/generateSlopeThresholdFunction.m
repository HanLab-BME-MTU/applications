function [ func ] = generateSlopeThresholdFunction( response, Korg, T, freq )
%generateSlopeThresholdFunction Generates a function that outputs the value of
%the first and second partial derivative of the orientation response
%function with respect to theta as well as the Jacobian

if(nargin < 4)
    freq = false;
end

if(freq)
    response_hat = response;
else
    response_hat = fft(response);
end

func = @slopeThresholdFunction;

    function [F,J] = slopeThresholdFunction(x)
        if(nargout > 1)
            derivs = 1:5;
        else
            derivs = 1:3;
        end
        D = 2*pi^2;
        theta = x(1,:);
        time = x(2,:);
        K = (1./sqrt(time)-1)/2;
        dt_dK = -4*time.^(1.5);
        p_dt_dK_dt = -6*time.^(0.5);
        disp(time);
        response_hat_at_K = orientationSpace.getResponseAtOrderVecHat(response_hat,Korg, K);
        values = interpft1_derivatives(response_hat_at_K,theta,derivs,2*pi,true);
        F = values(:,:,1:2);
        F(:,:,2) = dt_dK.*values(:,:,3) - T.*values(:,:,2);
%         F(:,:,2) = F(:,:,2).*(abs(dt_dK.*values(:,:,3)./values(:,:,2)) >= abs(T));
        F = permute(F,[ 3 2 1]);
%         F(1) = values(1);
%         F(2) = values(2);
        if(nargout > 1)
            J = values(:,:,[2 3 3 4]);
            J(:,:,2) = dt_dK.*values(:,:,4) - T.*values(:,:,3);
            J(:,:,4) = p_dt_dK_dt.*values(:,:,3) + D*dt_dK.*values(:,:,5) - T.*D.*values(:,:,4);
%             J(:,:,4) = dt_dK.*values(:,:,5) - T.*values(:,:,4);
%             J(:,:,4) = J(:,:,4)/2;
            J(:,:,3) = D*J(:,:,3);
%             J(:,:,3:4) = D*J(:,:,3:4);
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

