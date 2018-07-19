function matchOK=acceptTrackedSpotDist(match,model)
% sub-function of trackeframe
%test estimated distance between target spots for relevance above noise

val=ones(length(model),1);
T_TEST_PROB=0.05;
for i=1:length(model)
    for j=(i+1):length(model)
        Qxx=blkdiag(model(i).QMatrix, model(j).QMatrix);
        dif=((match(i).center+match(i).parms)-(match(j).center+match(j).parms));
        distance=sqrt(sum(dif.^2));
        H=1/distance*[dif -dif];
        Qdd=H*Qxx*H';
        %unlike in the slist, the chi is NOT squared here -> square
        chiT=(model(i).sigma0^2+model(j).sigma0^2)/2;
        sigD(i,j)=sqrt(Qdd*chiT);
        tvalue=distance/sigD(i,j);
        %keep only coords that pass as nonzero distance
        if tvalue<tinv(1-(T_TEST_PROB/2),1);
            val(i)=0;
            val(j)=0;
        end;
    end;
end;
matchOK=val;