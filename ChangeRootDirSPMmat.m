newrootdir = '/media/4A26DE8B26DE7781/data';
SPMdir = '/media/4A26DE8B26DE7781/data/MeynielMoreno_NACONF_2016/MRI_data/analyzed_data/subj02/first_level_estimates/Model01_LocalGlobal';
SPMname = 'SPM_laptopPath';

load(SPMname)
nScan = size(SPM.xY.P,1);
P = [];
for iScan = 1:nScan
    P(iScan,:) = [newrootdir, SPM.xY.P(iScan,33:end)];
end
SPM.xY.P = char(P);
SPM.swd = char([newrootdir, SPM.swd(1,33:end)]);
save(SPMname, 'SPM')