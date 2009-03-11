function [] = main283AUTO_cluster(iclean,oriImageName, oriImagePath);
% wrapper function for cluster

if iscell(oriImageName)
    for k=1:length(oriImageName)
        currImageName = oriImageName{k};
        currImagePath = oriImagePath{k};
        main283AUTO(iclean, currImageName, currImagePath);
    end
else
    main283AUTO(iclean, oriImageName, oriImagePath);
end


end % of function