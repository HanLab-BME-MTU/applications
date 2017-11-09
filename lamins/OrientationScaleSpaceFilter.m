classdef OrientationScaleSpaceFilter < OrientationSpaceFilter
    %OrientationScaleSpaceFilter implements a polar separable filter in
    %frequency space using an angular Gaussian for orientation and a
    %raised cos log filter (a la Simoncelli) for scale
    
    properties
        components
    end
    
    methods
        function obj = OrientationScaleSpaceFilter(f_c,b_f,K,normEnergy)
            if(nargin < 3)
                K = b_f;
                b_f = NaN;
            end
            if(nargin < 4 && ischar(K))
                normEnergy = K;
                K = b_f;
                b_f = NaN;
            end
            if(nargin < 4 || isempty(normEnergy))
                normEnergy = 'none';
            else
                if(iscell(normEnergy))
                    normEnergy = normEnergy{1};
                end
            end
            obj = obj@OrientationSpaceFilter(f_c, b_f, K, normEnergy);
            for ii=1:numel(obj)
                angles = 0:ceil(obj(ii).K); % is K an integer?
                angles = angles/obj(ii).n*pi;
                scales = abs(obj(ii).f_c./cos(angles));
                bandwidths = abs(obj(ii).b_f./cos(angles));
                obj(ii).components = OrientationSpaceFilter(scales,bandwidths,obj(ii).K,normEnergy);
                obj(ii).components = obj(ii).components(logical(eye(length(obj(ii).components))));
                % Preshift the components
%                 for jj=2:length(obj(ii).components)
%                     obj(ii).components(jj).circshiftAngles(jj-1);
%                 end
            end
        end
        function setupFilter(obj,siz)
            if(isscalar(siz))
                siz = siz([1 1]);
            end
%             coords = orientationSpace.getFrequencySpaceCoordinates(siz);
            notSetup = ~cellfun(@(x) isequal(siz,x),{obj.size});
            notSetup = notSetup | cellfun('isempty',{obj.F});
            obj = obj(notSetup);
            if(isempty(obj))
                return;
            end
            [obj.size] = deal(siz);
            
            % Setup all the filters at once since they probably share the
            % same angular order K
%             allComponents = [obj.components];
%             allComponents.setupFilter(siz);
            
            for ii=1:numel(obj)
                obj(ii).components.setupFilter(siz);
%                 obj(ii).components(1).setupFilter(siz);
                obj(ii).F = obj(ii).components(1).F;
%                 obj(ii).components(1).clearTransients();
                for jj=2:length(obj(ii).components)
%                     obj(ii).components(jj).setupFilter(siz)
                    obj(ii).F = obj(ii).F + circshift(obj(ii).components(jj).F, jj-1,3) ...
                                          + circshift(obj(ii).components(jj).F,-jj+1,3);
                    % We preshifted the first half already
%                     obj(ii).F = obj(ii).F +  ...
%                                 obj(ii).components(jj).F + ...
%                                 circshift(obj(ii).components(jj).F,-2*jj+2,3);
%                     obj(ii).components(jj).clearTransients();
                end
                obj(ii).components.clearTransients();
            end
            
%             allComponents.clearTransients();
            
%             if( all(obj(1).K == [obj.K]) )
%                 % angular component is all the same
%                 
%                 A = orientationSpace.angularKernel(obj(1).K,obj(1).angles,coords);
%                 R = orientationSpace.radialKernel([obj.f_c], [obj.b_f],coords);
%                 for o=1:numel(obj)
%                     obj(o).F = bsxfun(@times,A, R(:,:,o));
%                 end
%             else
%                 for o=1:numel(obj)
%                     obj(o).F = orientationSpace.kernel(obj(o).f_c, obj(o).b_f, obj(o).K, obj(o).angles, coords);
%                 end
%             end
%             for o=1:numel(obj)
%                 if(isempty(obj(o).normEnergy))
%                     break;
%                 end
%                 switch(obj(o).normEnergy)
%                     case 'energy'
%                         % E is complex
%                         E = shiftdim(obj(o).getEnergy(),-1);
%                         F = obj(o).F;
%                         obj(o).F = bsxfun(@rdivide,real(F),real(E)) +1j*bsxfun(@rdivide,imag(F),imag(E));
%                     case 'peak'
%                         F = obj(o).F;
%                         sumF = sum(F(:))./numel(F);
%                         obj(o).F = real(F)./real(sumF) + 1j*imag(F)./imag(sumF);
%                     case 'scale'
%                         obj(o).F = obj(o).F ./ obj(o).f_c ./ sqrt(siz(1)*siz(2));
%                     case 'sqrtscale'
%                         obj(o).F = obj(o).F ./ sqrt(obj(o).f_c) ./ sqrt(siz(1)*siz(2))
%                     case 'none'
%                     otherwise
%                         error('OrientationSpaceFilter:setupFilterNormEnergy', ...
%                             'Invalid normEnergy property');
%                 end
%             end
        end
    end
    methods (Static)
        function F = constructByRadialOrder(f_c, K_f, K, normEnergy, constructor)
            if(nargin < 4)
                normEnergy = [];
            end
            if(nargin < 5)
                constructor = @OrientationScaleSpaceFilter;
            end
            F = OrientationSpaceFilter.constructByRadialOrder(f_c, K_f, K, normEnergy, constructor);
        end
    end
    
end