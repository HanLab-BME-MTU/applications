
%%
% Crossing comparison

[idx_cell, dist_cell] = KDTreeBallQuery([Y1 X1],[Y1 X1],0.1);

crossing_flag_1 = zeros(1, length(X1));
angle_crossing_1 = nan(1, length(X1));

for iQ = 1 : length(Y1)
    idx_this = idx_cell{iQ};
    if(length(idx_this)>1)
        
        AA = O1(idx_this(1))-O1(idx_this(2));
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        
        if(abs(AA)>0.2)
            crossing_flag_1(iQ)=1;
            angle_crossing_1(iQ)=O1(idx_this(1))-O1(idx_this(2));
            ind_pair(iQ,1)=idx_this(1);
            ind_pair(iQ,2)=idx_this(2);
        end
    end
end
angle_crossing_1(angle_crossing_1>pi/2) = angle_crossing_1(angle_crossing_1>pi/2)-pi;
angle_crossing_1(angle_crossing_1<-pi/2) = angle_crossing_1(angle_crossing_1<-pi/2)+pi;
angle_crossing_1(angle_crossing_1>pi/2) = angle_crossing_1(angle_crossing_1>pi/2)-pi;
angle_crossing_1(angle_crossing_1<-pi/2) = angle_crossing_1(angle_crossing_1<-pi/2)+pi;



[idx_cell, dist_cell] = KDTreeBallQuery([Y2 X2],[Y2 X2],0.2);

crossing_flag_2 = zeros(1, length(X2));
angle_crossing_2 = nan(1, length(X2));

for iQ = 1 : length(Y2)
    idx_this = idx_cell{iQ};
    if(length(idx_this)>1)
        
        
         AA = O2(idx_this(1))-O2(idx_this(2));
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        if(AA>pi/2)
            AA=AA-pi;
        end
        if(AA<-pi/2)
            AA=AA+pi;
        end
        
        if(abs(AA)>0.2)
            crossing_flag_2(iQ)=1;
            angle_crossing_2(iQ)=O2(idx_this(2))-O2(idx_this(2));     
        end
    end    
end
angle_crossing_2(angle_crossing_2>pi/2) = angle_crossing_2(angle_crossing_2>pi/2)-pi;
angle_crossing_2(angle_crossing_2<-pi/2) = angle_crossing_2(angle_crossing_2<-pi/2)+pi;
angle_crossing_2(angle_crossing_2>pi/2) = angle_crossing_2(angle_crossing_2>pi/2)-pi;
angle_crossing_2(angle_crossing_2<-pi/2) = angle_crossing_2(angle_crossing_2<-pi/2)+pi;



[idx, dist] = KDTreeClosestPoint([Y2(crossing_flag_2>0) X2(crossing_flag_2>0)],...
    [Y1(crossing_flag_1>0) X1(crossing_flag_1>0) ]);

dist_pool_for_crossing = [dist_pool_for_crossing;dist];
ang_pool_for_crossing = [ang_pool_for_crossing;(angle_crossing_1(crossing_flag_1>0))'];


[idx, dist] = KDTreeClosestPoint([Y1(crossing_flag_1>0) X1(crossing_flag_1>0)],...
[Y2(crossing_flag_2>0) X2(crossing_flag_2>0)]);

dist_pool_for_crossing = [dist_pool_for_crossing;dist];
ang_pool_for_crossing = [ang_pool_for_crossing;(angle_crossing_2(crossing_flag_2>0))'];
