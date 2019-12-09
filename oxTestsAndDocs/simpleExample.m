%% Paths
% Add gannet to path
addpath('..')

% Add spm to path
% addpath('spmlocation here')

% Initilise spm jobman
spm('defaults','fmri')
spm_jobman('initcfg')

%% Data locations
mainFile = 'MPRESSExamples/Sub01/meas_MID01388_FID104333_mega_press_TR1500_LOC_Acq_CORR.dat';
waterRef = 'MPRESSExamples/Sub01/meas_MID01392_FID104337_mega_press_wref3_LOC.dat';
niftiFile = 'MPRESSExamples/Sub01/images_003_t1mprax1mmiso64chv21001.nii';

%% Call gannet processing
% loading and fitting
MRS_struct = GannetLoad({mainFile},{waterRef});
MRS_struct = GannetFit(MRS_struct);

% Quantification
MRS_struct = GannetCoRegister(MRS_struct,{niftiFile});
MRS_struct = GannetSegment(MRS_struct);
MRS_struct = GannetQuantify(MRS_struct);

%% Output
% Fit output can be found in MRS_struct.out.vox1 and more specifically
% MRS_struct.out.vox1.GABA and MRS_struct.out.vox1.Glx

% The values printed on the GannetQuantify output are
clc
tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_struct.out.vox1.GABA.ConcIU_TissCorr);
fprintf('%s Glx/Water: %.2f i.u.\n\n', tmp1, MRS_struct.out.vox1.Glx.ConcIU_TissCorr);

tmp1 = 'Relaxation-, tissue-, alpha-corrected (Harris et al. method)';
fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_struct.out.vox1.GABA.ConcIU_AlphaTissCorr);
fprintf('%s Glx/Water: %.2f i.u.\n\n', tmp1, MRS_struct.out.vox1.Glx.ConcIU_AlphaTissCorr);

tmp1 = 'Relaxation-, tissue-, alpha-corrected; average-voxel-normalized (Harris et al. method)';
fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_struct.out.vox1.GABA.ConcIU_AlphaTissCorr_GrpNorm);
fprintf('%s Glx/Water: %.2f i.u.\n\n', tmp1, MRS_struct.out.vox1.Glx.ConcIU_AlphaTissCorr_GrpNorm);

% The estimated fitting error can be found here:
fprintf('The estimated GABA fit error relative to water is %0.2f%%.\n',MRS_struct.out.vox1.GABA.FitError_W)

% I would recommed using the first output (Gasparovic) as this is what is
% used in the Big GABA II paper.

%% Output to cvs
% Optionally output to cvs by setting the following line befroe calling
% GannetQuantify. Note that the output directory is the current directory,
% so you might want to cd somewhere first (see 2nd following commented line)
% Note the Gasparovic corrected value is ConcIU_TissCorr in the csv file.

MRS_struct.p.csv = true;
% oldDir = cd('Outputdir'); % Location you want the cvs file output
MRS_struct = GannetQuantify(MRS_struct); % will now output cvs file
% cd(oldDir); % Return to original dir