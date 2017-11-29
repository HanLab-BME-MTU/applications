classdef OrientationSpaceWaveletFilter < OrientationSpaceFilter
    %ORIENTATIONSPACEWAVELETFILTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = OrientationSpaceWaveletFilter(f_c,b_f,K,normEnergy)
%             if(nargin == 0)
%                 return;
%             end
%             if(~isscalar(f_c) || ~isscalar(b_f) || ~isscalar(K))
%                 s = [length(f_c) length(b_f) length(K)];
%                 s(2) = max(1,s(2));
%                 f_c = repmat(f_c(:),1,s(2),s(3));
%                 if(isempty(b_f))
%                     b_f = 1/sqrt(2) * f_c;
%                 else
%                     b_f = repmat(b_f(:)',s(1),1,s(3));
%                 end
%                 K   = repmat(shiftdim(K(:),-2),s(1),s(2),1);
%                 if(nargin < 4)
%                         normEnergy = [];
%                 end
%                 normEnergy = repmat({normEnergy},size(K));
%                 constructor = str2func(class(obj));
%                 obj = arrayfun(constructor,f_c,b_f,K,normEnergy,'UniformOutput',false);
%                 obj = reshape([obj{:}],size(obj));
%                 return;
%             end
%             if(isempty(b_f))
%                 % Set the bandwidth to be 0.8 of the central frequency by
%                 % default
%                 b_f = 1/sqrt(2) * f_c;
%             end
            if(nargin < 4 || isempty(normEnergy))
                normEnergy = 'none';
            else
                if(iscell(normEnergy))
                    normEnergy = normEnergy{1};
                end
            end
%             obj.f_c = f_c;
%             obj.b_f = b_f;
%             obj.K = K;
%             obj.normEnergy = normEnergy;
%             
% %             obj.n = 2*ceil(K) + 1;
%             obj.n = 2*obj.sampleFactor*ceil(K) + 1;
% %             obj.angles = (0:obj.n-1)/obj.n*pi;
% %             obj.angles = 0:pi/obj.n:pi-pi/obj.n;
%             
            obj = obj@OrientationSpaceFilter(f_c,b_f,K,normEnergy);
        end
    end
    
    methods
        function setupFilter(obj,siz)
            if(isscalar(siz))
                siz = siz([1 1]);
            end
            coords = orientationSpace.getFrequencySpaceCoordinates(siz);
            notSetup = ~cellfun(@(x) isequal(siz,x),{obj.size});
            notSetup = notSetup | cellfun('isempty',{obj.F});
            obj = obj(notSetup);
            if(isempty(obj))
                return;
            end
            [obj.size] = deal(siz);
            if( all(obj(1).K == [obj.K]) )
                % angular component is all the same
                
                A = orientationSpace.angularKernel(obj(1).K,obj(1).angles,coords);
                R = orientationSpace.radialWaveletKernel([obj.f_c], [obj.b_f],coords);
                for o=1:numel(obj)
                    obj(o).F = bsxfun(@times,A, R(:,:,o));
                end
            else
                for o=1:numel(obj)
                    obj(o).F = orientationSpace.kernel(obj(o).f_c, obj(o).b_f, obj(o).K, obj(o).angles, coords);
                end
            end
            for o=1:numel(obj)
                if(isempty(obj(o).normEnergy))
                    break;
                end
                switch(obj(o).normEnergy)
                    case 'energy'
                        % E is complex
                        E = shiftdim(obj(o).getEnergy(),-1);
                        F = obj(o).F;
                        obj(o).F = bsxfun(@rdivide,real(F),real(E)) +1j*bsxfun(@rdivide,imag(F),imag(E));
                    case 'peak'
                        F = obj(o).F;
                        sumF = sum(F(:))./numel(F);
                        obj(o).F = real(F)./real(sumF) + 1j*imag(F)./imag(sumF);
                    case 'scale'
                        obj(o).F = obj(o).F ./ obj(o).f_c ./ sqrt(siz(1)*siz(2));
                    case 'sqrtscale'
                        obj(o).F = obj(o).F ./ sqrt(obj(o).f_c) ./ sqrt(siz(1)*siz(2));
                    case 'none'
                    otherwise
                        error('OrientationSpaceFilter:setupFilterNormEnergy', ...
                            'Invalid normEnergy property');
                end
            end
        end
    end
    methods (Static)
        function F = constructEqualLengthFilters(f_c, b_f, K, normEnergy, constructor)
            if(nargin < 5)
                constructor = @OrientationSpaceWaveletFilter;
            end
            %% Approximate cone by height of triangle
            % Largest central frequency or smallest scale
            f_c_max = max(f_c(:));
%             height = sin(pi/(2*K+1))*f_c_max;
            arcLength = pi/(2*K+1)*f_c_max;
            
            assert(isscalar(K));
            
            %% Normal constructor
            s = [length(f_c) length(b_f) 1];
            s(2) = max(1,s(2));
            f_c = repmat(f_c(:),1,s(2),s(3));
            if(isempty(b_f))
                b_f = 1/sqrt(2) * f_c;
            else
                b_f = repmat(b_f(:)',s(1),1,s(3));
            end
            
            if(nargin < 4)
                    normEnergy = [];
            end
            normEnergy = repmat({normEnergy},size(f_c));
            
%             K = (pi./asin(height ./ f_c) - 1)./2;
            K = (pi/arcLength*f_c-1)/2;
            
            F = arrayfun(constructor,f_c,b_f,K,normEnergy,'UniformOutput',false);
            F = reshape([F{:}],size(F));

        end
        function F = constructByRadialOrder(f_c, K_f, K, normEnergy, constructor)
            if(nargin < 4)
                normEnergy = [];
            end
            if(nargin < 5)
                constructor = @OrientationSpaceWaveletFilter;
            end
            b_f = f_c ./ sqrt(K_f);
            F = constructor(f_c, b_f, K, normEnergy);
            F = F(logical(eye(length(f_c))));
        end
    end

    
end

