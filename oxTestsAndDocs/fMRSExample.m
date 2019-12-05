%% Paths
% Add gannet to path
addpath('..')

% Add spm to path
% addpath('spmlocation here')

% Initilise spm jobman
spm('defaults','fmri')
spm_jobman('initcfg')

%% Data locations
mainFile = '/Users/wclarke/Documents/Data/tmp/BI_MEGAPRESS_Comparison/2013_40_525/meas_MID00241_FID57300_mega_press_TR1500_2.dat';
waterRef = '/Users/wclarke/Documents/Data/tmp/BI_MEGAPRESS_Comparison/2013_40_525/meas_MID00245_FID57304_mega_press_wref3.dat';
niftiFile = '/Users/wclarke/Documents/Data/tmp/BI_MEGAPRESS_Comparison/2013_40_525/F3T_2013_40_525/images_015_t1mprax1mmisowithNose64ch1001.nii';

%% Output directory
mkdir('fMRSTestOutput')

%% Call gannet load
% loading and fitting
originalDir = cd('fMRSTestOutput');
MRS_struct = GannetLoad({mainFile},{waterRef});

%% Split the MRS_struct across blocks
% In this example I'm using 4 equal size blocks.
blocks = 4;
averages = MRS_struct.p.Navg;
blockSize = MRS_struct.p.Navg/blocks;
blockStep = blockSize; % no overlap
fprintf('There are %0.0f averages, block size will be %0.0f with a step size of %0.0f.\n',averages,blockSize,blockStep);

MRS_structs = GannetChop(MRS_struct,blockSize,blockStep);
%% Show subdivided spectra
figure(1000)
clf
subplot(1,1+blocks,1)
GannetPlotPrePostAlign(MRS_struct,MRS_struct.p.Vox, 1, 1);
title('Combined')
for iDx = 1:blocks
    subplot(1,1+blocks,iDx+1)
    GannetPlotPrePostAlign(MRS_structs{iDx},MRS_structs{iDx}.p.Vox, 1, 1);
    title(sprintf('Block %0.0f',iDx))
end
%% Process the blocks
clc
for iDx = 1:blocks
    % Output folder
    mkdir(sprintf('Block_%03.0f',iDx))
    blockDir = cd(sprintf('Block_%03.0f',iDx));
    
    MRS_structs{iDx} = GannetFit(MRS_structs{iDx});

    % Quantification
    MRS_structs{iDx} = GannetCoRegister(MRS_structs{iDx},{niftiFile});
    MRS_structs{iDx} = GannetSegment(MRS_structs{iDx});
    
    MRS_structs{iDx}.p.csv = true;
    MRS_structs{iDx} = GannetQuantify(MRS_structs{iDx});
    
    tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
    fprintf('Results of block %0.0f.\n',iDx)
    fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_structs{iDx}.out.vox1.GABA.ConcIU_TissCorr);
    fprintf('%s Glx/Water: %.2f i.u.\n', tmp1, MRS_structs{iDx}.out.vox1.Glx.ConcIU_TissCorr);
    fprintf('The estimated GABA fit error relative to water is %0.2f%%.\n',MRS_structs{iDx}.out.vox1.GABA.FitError_W)

    cd(blockDir)
end
%%
cd(originalDir)