%% Paths
% Add gannet to path
addpath('..')

% Add spm to path
% addpath('spmlocation here')

% Initilise spm jobman
spm('defaults','fmri')
spm_jobman('initcfg')

%% Data locations
% List batched files in cell array. Each column is a subject/session/voi
% etc. 

mainFiles = {'MPRESSExamples/Sub01/meas_MID01388_FID104333_mega_press_TR1500_LOC_Acq_CORR.dat',...
             'MPRESSExamples/Sub02/meas_MID01067_FID107301_mega_press_TR1500_LOC_Acq_1_CORR.dat'};
waterRefs = {'MPRESSExamples/Sub01/meas_MID01392_FID104337_mega_press_wref3_LOC.dat',...
             'MPRESSExamples/Sub02/meas_MID01071_FID107305_mega_press_wref3_LOC.dat'};
niftiFiles = {'MPRESSExamples/Sub01/images_003_t1mprax1mmiso64chv21001.nii',...
              'MPRESSExamples/Sub02/images_002_t1mprax1mmiso64chv21001.nii'};

%% Call gannet processing
% loading and fitting
MRS_struct = GannetLoad(mainFiles,waterRefs);
MRS_struct = GannetFit(MRS_struct);

% Quantification
MRS_struct = GannetCoRegister(MRS_struct,niftiFiles);
MRS_struct = GannetSegment(MRS_struct);

MRS_struct.p.csv = true; % Output to cvs
MRS_struct = GannetQuantify(MRS_struct);

%% Output
tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
for iDx = 1:numel(mainFiles)
    fprintf('Subject %02.0f:\n',iDx)
    fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_struct.out.vox1.GABA.ConcIU_TissCorr(iDx));
    fprintf('%s Glx/Water: %.2f i.u.\n', tmp1, MRS_struct.out.vox1.Glx.ConcIU_TissCorr(iDx));
    fprintf('The estimated GABA fit error relative to water is %0.2f%%.\n\n',MRS_struct.out.vox1.GABA.FitError_W(iDx))
end