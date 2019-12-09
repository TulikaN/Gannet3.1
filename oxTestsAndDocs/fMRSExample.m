%% Paths
cd(fileparts(mfilename('fullpath')))
% Add gannet to path
addpath('..')

% Add spm to path
% addpath('spmlocation here')

% Initilise spm jobman
spm('defaults','fmri')
spm_jobman('initcfg')

%% Data locations
mainFiles = {'MPRESSExamples/Sub01/meas_MID01388_FID104333_mega_press_TR1500_LOC_Acq_CORR.dat',...
             'MPRESSExamples/Sub02/meas_MID01067_FID107301_mega_press_TR1500_LOC_Acq_1_CORR.dat'};
waterRefs = {'MPRESSExamples/Sub01/meas_MID01392_FID104337_mega_press_wref3_LOC.dat',...
             'MPRESSExamples/Sub02/meas_MID01071_FID107305_mega_press_wref3_LOC.dat'};
niftiFiles = {'MPRESSExamples/Sub01/images_003_t1mprax1mmiso64chv21001.nii',...
              'MPRESSExamples/Sub02/images_002_t1mprax1mmiso64chv21001.nii'};

% Convert ot absolute paths
mainFiles= fullfile(pwd,mainFiles);
waterRefs= fullfile(pwd,waterRefs);
niftiFiles= fullfile(pwd,niftiFiles);

%% Output directory
mkdir('fMRSTestOutput')

%% Call gannet load
% loading and fitting
originalDir = cd('fMRSTestOutput');
MRS_struct = GannetLoad(mainFiles,waterRefs);

%% Split the MRS_struct across blocks
% In this example I'm using 4 equal size blocks.
blocks = 4;
averages = MRS_struct.p.Navg;
blockSize = MRS_struct.p.Navg/blocks;
blockStep = blockSize; % no overlap
for iDx = 1:numel(averages)
    fprintf('There are %0.0f averages, block size will be %0.0f with a step size of %0.0f.\n',averages(iDx),blockSize(iDx),blockStep(iDx));
end

MRS_structs = GannetChop(MRS_struct,blockSize(1),blockStep(1));
%% Show subdivided spectra
for jDx = 1:MRS_struct.p.numscans
    figure(1000+jDx)
    clf
    subplot(1,1+blocks,1)
    GannetPlotPrePostAlign(MRS_struct,MRS_struct.p.Vox, 1, 1);
    title('Combined')
    for iDx = 1:blocks
        subplot(1,1+blocks,iDx+1)
        GannetPlotPrePostAlign(MRS_structs{iDx},MRS_structs{iDx}.p.Vox, jDx, 1);
        title(sprintf('Block %0.0f',iDx))
    end
end

%% Plot the mean of the blocks
% There should only be significant differences in regions containing
% lipids/water which have been removed per block.
meanSpec = [];
for iDx = 1:blocks
    meanSpec(:,:,iDx) = MRS_structs{iDx}.spec.vox1.GABAGlx.diff;
end
meanSpec = mean(meanSpec,3);

figure(1500)
clf
subplot(1,2,1)
plot(MRS_struct.spec.freq,real(MRS_struct.spec.vox1.GABAGlx.diff(1,:)))
hold on
plot(MRS_struct.spec.freq,real(meanSpec(1,:)))
xlim([0 5]);set(gca,'xdir','reverse');legend('Original','combined blocks'); title('Subject 1')
subplot(1,2,2)
plot(MRS_struct.spec.freq,real(MRS_struct.spec.vox1.GABAGlx.diff(2,:)))
hold on
plot(MRS_struct.spec.freq,real(meanSpec(2,:)))
xlim([0 5]);set(gca,'xdir','reverse');legend('Original','combined blocks'); title('Subject 2')

%% Process the blocks
clc
for iDx = 1:blocks
    % Output folder
    mkdir(sprintf('Block_%03.0f',iDx))
    blockDir = cd(sprintf('Block_%03.0f',iDx));
    
    MRS_structs{iDx} = GannetFit(MRS_structs{iDx});

    % Quantification
    MRS_structs{iDx} = GannetCoRegister(MRS_structs{iDx},niftiFiles);
    MRS_structs{iDx} = GannetSegment(MRS_structs{iDx});
    
    MRS_structs{iDx}.p.csv = true;
    MRS_structs{iDx} = GannetQuantify(MRS_structs{iDx});
    
    tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
    fprintf('Results of block %0.0f.\n',iDx)
    for sDx = 1:numel(mainFiles)
        fprintf('Subject %02.0f:\n',sDx)
        fprintf('%s GABA+/Water: %.2f i.u.\n', tmp1, MRS_structs{iDx}.out.vox1.GABA.ConcIU_TissCorr(sDx));
        fprintf('%s Glx/Water: %.2f i.u.\n', tmp1, MRS_structs{iDx}.out.vox1.Glx.ConcIU_TissCorr(sDx));
        fprintf('The estimated GABA fit error relative to water is %0.2f%%.\n\n',MRS_structs{iDx}.out.vox1.GABA.FitError_W(sDx))
    end

    cd(blockDir)
end
%%
cd(originalDir)

%% Plot
for iDx = 1:numel(MRS_structs)
    GABAConc(iDx,:) = MRS_structs{iDx}.out.vox1.GABA.ConcIU_TissCorr;
    GABASD(iDx,:) = MRS_structs{iDx}.out.vox1.GABA.FitError_W;
end

figure(2000)
clf
hold on
for iDx = 1:size(GABAConc,2)
    currentSD = GABAConc(:,iDx).*GABASD(:,iDx)/100;
    errorbar(1:blocks,GABAConc(:,iDx),currentSD,currentSD,'color',[0.5 0.5 0.5])
end
errorbar(1:blocks,mean(GABAConc,2),std(GABAConc,[],2),std(GABAConc,[],2),'k')

