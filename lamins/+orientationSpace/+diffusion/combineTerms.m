function [ ss ] = combineTerms( s )
%combineTerms Combine terms from calcThetaDeriv or calcTimeDeriv

    max_rhod = max(arrayfun(@(s) size(s.rhod,2),s));
    has_timed = isfield(s,'timed');
    has_thetad = isfield(s,'thetad');
    if(has_timed)
        max_timed = max(arrayfun(@(s) size(s.timed,2),s));
    end
    if(has_thetad)
        max_thetad = max(arrayfun(@(s) size(s.thetad,2),s));
    end
    for i=1:length(s)
        s(i).rhod(end+1:max_rhod) = 0;
        if(has_timed)
            s(i).timed(end+1:max_timed) = 0;
        end
        if(has_thetad)
            s(i).thetad(end+1:max_thetad) = 0;
        end
    end
    ss.rhod = cat(1,s.rhod);
    ss.coeff = cat(1,s.coeff);
    if(has_timed)
        ss.timed = cat(1,s.timed);
    end
    if(has_thetad)
        ss.thetad = cat(1,s.thetad);
    end
    ss.D = cat(1,s.D);

end

