function tmp=cleantmp(tmp)

c=1;
for i=1:size(tmp,1)
        
        start=tmp(:,1:2);
        stop=start; %tmp(:,3:4);
        
        extr=[]; rem=[];
            
        % Find repetitions
        t=start(i,1)==stop(:,1);
        u=start(i,2)==stop(:,2);
        y=find(t & u);
        
        if length(y)==1
            
            % Okay. This is the only occurrence
            
        elseif length(y)>1
            extr=[];
            extr=tmp(y,:);
            
            rem=find(tmp(y,3)==0);
            
            if length(rem)==length(y)
                rem=rem(1);
            end
           
            toBeRem(c:c+length(rem)-1)=y(rem)';
            c=c+length(rem);
            
        end
end

if ~isempty(toBeRem)
    tmp(toBeRem,:)=[];
end

            
