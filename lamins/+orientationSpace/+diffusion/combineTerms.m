function [ ss ] = combineTerms( s )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    max_rhod = max(arrayfun(@(s) size(s.rhod,2),s));
    max_timed = max(arrayfun(@(s) size(s.timed,2),s));
    for i=1:length(s)
        s(i).rhod(end+1:max_rhod) = 0;
        s(i).timed(end+1:max_timed) = 0;
    end
    ss.rhod = cat(1,s.rhod);
    ss.coeff = cat(1,s.coeff);
    ss.timed = cat(1,s.timed);
    ss.D = cat(1,s.D);

end

