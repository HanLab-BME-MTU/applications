function [ maxima ] = traceAllLocalMaxima( obj, K , maxima_org)
%traceAllLocalMaxima Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    maxima_org = obj.getRidgeOrientationLocalMaxima;
end

    maxima_org = shiftdim(maxima_org,2)*2;
    
    maxima = cell(1,length(K));
    
    [K,sort_idx] = sort(K,'descend');
    
    Kidx = 1:length(K);
    Kidx(obj.filter.K ==K) = [];

    pool = gcp('nocreate');
    if(isempty(pool))
        prev_maxima = maxima_org;
        for ki=Kidx
            maxima{ki} = halleyft(shiftdim(obj.getResponseAtOrderFT(K(ki),2).a,2),prev_maxima,false,1);
            prev_maxima = maxima{ki};
        end
%         maxima{obj.filter.K ==K} = maxima_org;
    else
        maxima_org_sz = size(maxima_org);
        prev_maxima = maxima_org(:,:);

        pool_size = pool.NumWorkers;
        block_size = ceil(size(prev_maxima,2)/pool_size);
        partition = [repmat(block_size,1,floor(size(prev_maxima,2)/block_size)) mod(size(prev_maxima,2),block_size)];
        partition = partition(partition ~= 0);
        
        prev_maxima = mat2cell(prev_maxima,maxima_org_sz(1),partition);
        
        response = cell(length(K),length(partition));
        for ki=Kidx
            response{ki,1} = shiftdim(obj.getResponseAtOrderFT(K(ki),2).a,2);
            response{ki,1} = response{ki,1}(:,:);
            response(ki,:) = mat2cell(response{ki,1},size(response{ki},1),partition);
        end
        response = num2cell(response,1);
        
        maxima = parcellfun_progress(@(r,m) halleyftDescent(r,m,Kidx),response,prev_maxima,'UniformOutput',false);
        maxima = vertcat(maxima{:});
        maxima = num2cell(maxima,1);
        
%         for ki=Kidx
%             response = shiftdim(obj.getResponseAtOrderFT(K(ki),2).a,2);
%             response = response(:,:);
%             response = mat2cell(response,size(response,1),partition);
%             maxima{ki} = parcellfun_progress(@(r,m) halleyft(r,m,false,1),response,prev_maxima,'UniformOutput',false);
%             prev_maxima = maxima{ki};
%         end
        for ki=Kidx
            maxima{ki} = [maxima{ki}{:}];
            maxima{ki} = reshape(maxima{ki},maxima_org_sz);
        end
%         maxima{obj.filter.K ==K} = maxima_org;
    end
    maxima{obj.filter.K ==K} = maxima_org;
    
    maxima(sort_idx) = maxima;
    maxima = cat(4,maxima{:});
    maxima = maxima/2;
    maxima = permute(maxima,[2 3 1 4]);

end
function maxima = halleyftDescent(response,prev_maxima,Kidx)
    maxima = cell(1,length(response));
    for ki=Kidx
        maxima{ki} = halleyft(response{ki},prev_maxima,false,1);
        prev_maxima = maxima{ki};
    end
end