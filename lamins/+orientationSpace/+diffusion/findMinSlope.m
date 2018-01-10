function [ x_out ] = findMinSlope( x0, response, K_response, K_min , freq)
%findMinSlope Find minimum abs( dtheta_dK)

    if(nargin < 5)
        freq = false;
    end

    if(freq)
        response_hat = response;
    else
        response_hat = fft(response);
    end
    
    lb = [-Inf; K_min];
    ub = [Inf; K_response];
    
    lb = repmat(lb,1,size(x0,2));
    ub = repmat(ub,1,size(x0,2));
    
    options = optimoptions('fmincon');
    options.SpecifyObjectiveGradient = true;
    options.SpecifyConstraintGradient = true;
    options.CheckGradients = true;
    
    x_out = fmincon(@objectiveFunc,x0, [], [], [], [], lb, ub, @constraintFunc,options);


    function [F,J] = objectiveFunc(x)
        % (dtheta_dK).^2
        if(nargout > 1)
            derivs = 1:5;
        else
            derivs = 1:3;
        end
        D = 2*pi^2;
        theta = x(1,:);
        K = x(2,:);
        dt_dK   = -4./(2*K+1).^3;
        d2t_dK2 = 24./(2*K+1).^4;
        response_hat_at_K = orientationSpace.getResponseAtOrderVecHat(response_hat,K_response, K);
        values = interpft1_derivatives(response_hat_at_K,theta,derivs,2*pi,true);
        sqrtF = -D .* values(:,:,3) ./ values(:,:,2) .* dt_dK;
        F = sqrtF.^2;
%         F = D.^2 .* values(:,:,3).^2 ./ values(:,:,2).^2 .* (dt_dK).^2;
        % Convert into proper format
%         F = permute(F,[ 3 2 1]);
        if(nargout > 1) 
            % partial with respect to theta
            J(:,:,1) = values(:,:,4)    ./ values(:,:,2)    - ...
                       values(:,:,3).^2 ./ values(:,:,2).^2;
            J(:,:,1) = -D .* dt_dK .* J(:,:,1);
            J(:,:,1) = 2 .* sqrtF(:,:,1) .* J(:,:,1);
            % partial with respect to K
            J(:,:,2) = (values(:,:,5) .* dt_dK.^2 .* D + values(:,:,3).* d2t_dK2)./values(:,:,2) - ...
                       (values(:,:,3) .* dt_dK.^2 .* D .* values(:,4))./values(:,:,2).^2;
            J(:,:,2) = -D .* J(:,:,2);
            J(:,:,2) = 2 .* sqrtF(:,:,1) .* J(:,:,2);
            % Convert into proper format
            J = permute(J,[3 2 1]);
%             J = reshape(J,[],2);
%             J = sum(J);
%             J = reshape(J,[],2,2);
%             J = permute(J,[2 3 1]);
%             J = num2cell(J,[1 2]);
%             J = blkdiag(J{:});
        end
        F = sum(F(:));
    end

    function [F1,F2,J1,J2] = constraintFunc(x)
        % p2_rho_p_theta2 (should be less than zero), p_rho_p_theta (should be zero)
        if(nargout > 1)
            derivs = 1:4;
        else
            derivs = 1:2;
        end
        D = 2*pi^2;
        theta = x(1,:);
        K = x(2,:);
        dt_dK   = -4./(2*K+1).^3;
        response_hat_at_K = orientationSpace.getResponseAtOrderVecHat(response_hat,K_response, K);
        values = interpft1_derivatives(response_hat_at_K,theta,derivs,2*pi,true);
        F1 = values(:,:,2);
        F1 = permute(F1,[ 3 2 1]);
        if(nargout > 2)
            J1(:,:,1) = values(:,:,3);
            J1(:,:,2) = D.* values(:,:,4) .* dt_dK;
            % Convert into proper format
            J1 = reshape(J1,[],2,1);
            J1 = permute(J1,[2 3 1]);
            J1 = num2cell(J1,[1 2]);
            J1 = blkdiag(J1{:});
%             J1 = J1(:);
        end

%         D = 2*pi^2;
%         theta = x(1,:);
%         K = x(2,:);
%         dt_dK   = -4./(2*K+1).^3;
%         response_hat_at_K = orientationSpace.getResponseAtOrderVecHat(response_hat,K_response, K);
        values = interpft1_derivatives(response_hat_at_K,theta,derivs,2*pi,true);
        F2 = values(:,:,1);
        F2 = permute(F2,[ 3 2 1]);
        if(nargout > 3)
            J2(:,:,1) = values(:,:,2);
            J2(:,:,2) = D.* values(:,:,3) .* dt_dK;
            % Convert into proper format
            J2 = reshape(J2,[],2,1);
            J2 = permute(J2,[2 3 1]);
            J2 = num2cell(J2,[1 2]);
            J2 = blkdiag(J2{:});
%             J2 = J2(:);
        end
    end


end

