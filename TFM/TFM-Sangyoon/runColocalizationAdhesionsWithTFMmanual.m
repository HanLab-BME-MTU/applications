% runColocalizationAdhesionsWithTFMmanual.m for figure 2o in Han Nature
% Methods paper
% [forceNA_fttc,forceFC_fttc,forceFA_fttc,forceBG_fttc] = colocalizationAdhesionsWithTFMmanual( '/project/cellbiology/gdanuser/adhesion/Sangyoon/Youbean/130429 cell11 1frame/ROI','L2 Lcorner',16,1600);
% [forceNA_opt,forceFC_opt,forceFA_opt,forceBG_opt] = colocalizationAdhesionsWithTFMmanual( '/project/cellbiology/gdanuser/adhesion/Sangyoon/Youbean/130429 cell11 1frame/ROI','L2 optimal',16,1600);
% [forceNA_L1,forceFC_L1,forceFA_L1,forceBG_L1] = colocalizationAdhesionsWithTFMmanual( '/project/cellbiology/gdanuser/adhesion/Sangyoon/Youbean/130429 cell11 1frame/ROI','L1new',16,1600);

pathd = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft Nat Methods/Revision1/Figures/Fig2/cell1/ROI';
[forceNA_fttc,forceFC_fttc,forceFA_fttc,forceBG_fttc] =...
    colocalizationAdhesionsWithTFMmanual( pathd,'L2 Lcorner',16,1600);

[forceNA_opt,forceFC_opt,forceFA_opt,forceBG_opt] = ...
    colocalizationAdhesionsWithTFMmanual( pathd,'L2 optimal',16,1600);

[forceNA_L1,forceFC_L1,forceFA_L1,forceBG_L1] = ...
    colocalizationAdhesionsWithTFMmanual( pathd,'L1new',16,1600);
save([pathd filesep 'forceBGNAFCFA.mat'])