function ncord = testdistance(cord,Q,chi,snr,dataProperties)
% TESTDISTANCE test the distance if signifcantly larger than zero
%
% SYNOPSIS ncord = testdistance(cord,Q,chi,snr)
%
% INPUT   cord   : list of all coordinates
%              Q       : Q matrix of fitted model
%              chi     : chi squared of fitted model 
%              snr    :  snr of spots
%
% OUTPUT ncord : significant coordinates 


%c: 04/07/01 dT

%CONST DEFINITIONS
T_TEST_PROB=dataProperties.T_TEST_PROB;

% crosscheck every pair for zero distance
% statistical t-test 
nspots=length(cord.sp);
valid=[1:nspots];
for i=1:(nspots-1)
    for j=(i+1):nspots
        % distinguish between evtl. overlapped spots and separated spots
        if j==i+1
            Qxx=Q((i-1)*3+1:(i-1)*3+6,(i-1)*3+1:(i-1)*3+6);
        else
            Qxx=blkdiag(Q((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3), Q((j-1)*3+1:(j-1)*3+3,(j-1)*3+1:(j-1)*3+3));
        end;
       % dif=(cord.sp(i).cord-cord.sp(j).cord).*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
        dif=(cord.sp(i).cord-cord.sp(j).cord);
        distance=sqrt(sum(dif.^2));
        H=1/distance*[dif -dif];
        Qdd=H*Qxx*H';
        chiT=(chi(i)+chi(j))/2;
        sigD(i,j)=sqrt(Qdd*chiT);
        tvalue=distance/sigD(i,j);
        %keep only coords that pass as nonzero distance
        if tvalue<tinv(1-(T_TEST_PROB/2),1);
            % spot i from multi check, remove it
            if cord.sp(i).mult    
                rElem=i;
                ind=find(valid==rElem);
            elseif cord.sp(j).mult
                rElem=j;
                ind=find(valid==rElem);
            else
                ind=[];
            end;
            
            if ~isempty(ind)
                valid=[valid(1:ind-1) valid((ind+1):end)];
            end;
        end;
    end;
end;

sigD(nspots,nspots)=0;
% copy back accepted coords
ncord.sp=[];
ct=1;
Qtemp=[];
for i=valid
    ncord.sp(ct).cord=cord.sp(i).cord;
    ncord.sp(ct).amp=cord.sp(i).amp;
    ncord.sp(ct).bg=cord.sp(i).bg;
    ncord.sp(ct).mnint=cord.sp(i).mnint;
    ncord.sp(ct).mult=cord.sp(i).mult;
    %ncord.sp(ct).parms=cord.sp(i).parms;
    ct=ct+1;
    Qtemp=blkdiag(Qtemp,Q((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3));
end;
ncord.mnint=cord.mnint;
ncord.statistics.sigD=sigD(valid,valid);
ncord.statistics.chi=chi(valid);
ncord.statistics.snr=snr(valid);
ncord.statistics.Q=Qtemp;
ncord.parms=cord.parms;
