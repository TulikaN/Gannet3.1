function outStructs = GannetChop(MRS_struct,blocksize,blockstep)
% Function to chop up a gannet MRS_struct loaded by GannetLoad into 
% multiple MRS_structs with different numbers of averages.
% Input:
% MRS_struct - MRS_struct from GannetLoad
% blocksize - size of blocks, must be even.
% blockstep - number of averages between blocks.

% Will Clarke, University of Oxford, 2019.

% Currently this ONLY works with GABAGlx as a target. HERMES isn't
% supported either.
if numel(MRS_struct.p.target)>1 || ~strcmp(MRS_struct.p.target,'GABAGlx')
    error('Only GABAGlx is supported by GannetChop currently')
end

if MRS_struct.p.HERMES || MRS_struct.p.HERCULES
    error('HERMES and HERCULES not supported')
end
%% Block size maths.
if mod(blocksize,2)
    error('This is an interleaved MEGA edited sequence so blocksize must be even.')
end

if mod(blockstep,2)
    error('This is an interleaved MEGA edited sequence so blockstep must be even.')
end

if (blocksize<0) || (blockstep<0)
    error('Both blocksize and blockstep must be > 0');
end

% Calculate indicies.
currentAverage = 1;
numberOfAverages =  MRS_struct.p.Navg;
blockIndicies = logical([]); % Need the explicet cast otherwise tmp gets cast back to a double on assignment.
while (currentAverage+blocksize-1) <= numberOfAverages
    tmp = false(numberOfAverages,1);
    tmp(currentAverage:(currentAverage+blocksize-1)) = true;
    blockIndicies(:,end+1) =  tmp;
    currentAverage = currentAverage + blockstep;
end
numBlocks = size(blockIndicies,2);

%% Make the new structures

for iDx = 1:numBlocks
    for ii = 1:MRS_struct.p.numscans % Loop of batched files
        outStructs{iDx} = MRS_struct;

        % Go through each field and touch up the values
        % p field
        outStructs{iDx}.p.nrows = blocksize;
        outStructs{iDx}.p.Navg = blocksize;

        % spec field
        outStructs{iDx}.spec.F0freq = MRS_struct.spec.F0freq(ii,blockIndicies(:,iDx));
        outStructs{iDx}.spec.F0freq2 = MRS_struct.spec.F0freq2(ii,blockIndicies(:,iDx));
        outStructs{iDx}.spec.AllFramesFTrealign = MRS_struct.spec.AllFramesFTrealign(:,blockIndicies(:,iDx));

        % fids field
        outStructs{iDx}.fids.data = MRS_struct.fids.data(:,blockIndicies(:,iDx));
        outStructs{iDx}.fids.ON_OFF = MRS_struct.fids.ON_OFF(blockIndicies(:,iDx));
        outStructs{iDx}.fids.data_align = MRS_struct.fids.data_align(:,blockIndicies(:,iDx));
    end
    
    % fid and spec VoxX fields
    [outStructs{iDx}.spec.vox1.GABAGlx,currReject] = redoAveraging(outStructs{iDx});
    outStructs{iDx} = redoWaterRemoval(outStructs{iDx});
    
    % out field
    outStructs{iDx}.out.SpecReg = structfun(@(x) x(blockIndicies(:,iDx)),MRS_struct.out.SpecReg,'UniformOutput',false);
    outStructs{iDx}.out.reject = currReject;
    
end

end % End of function

% Function which replicates lines 601-651 of GannetLoad
% require to populate the voxX fields of MRS_struct.spec
function [out,reject] = redoAveraging(MRS_struct)

AllFramesFTrealign = MRS_struct.spec.AllFramesFTrealign;

for ii = 1:MRS_struct.p.numscans % Loop of batched files
    freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
    if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg') && ~strcmp(MRS_struct.p.vendor,'Siemens_rda') % if .rda data, use conventional averaging

        % Determine weights for weighted averaging (MM: 190423)
        ON_inds  = find(MRS_struct.fids.ON_OFF == 1);
        OFF_inds = find(MRS_struct.fids.ON_OFF == 0);
        DIFFs = zeros(size(AllFramesFTrealign,1), size(AllFramesFTrealign,2)/2);

        for ll = 1:size(AllFramesFTrealign,2)/2
            DIFFs(:,ll) = AllFramesFTrealign(:,ON_inds(ll)) - AllFramesFTrealign(:,OFF_inds(ll));
        end

        DIFFs = ifft(ifftshift(DIFFs,1),[],1);
        DIFFs = fftshift(fft(DIFFs(1:MRS_struct.p.npoints(ii),:),[],1),1);

        freq = (size(DIFFs,1) + 1 - (1:size(DIFFs,1))) / size(DIFFs,1) * freqRange + 4.68 - freqRange/2;
        freqLim = freq <= 4.25 & freq >= 1.8;
        D = zeros(size(AllFramesFTrealign,2)/2);

        for ll = 1:size(AllFramesFTrealign,2)/2
            for mm = 1:size(AllFramesFTrealign,2)/2
                tmp = sum((real(DIFFs(freqLim,ll)) - real(DIFFs(freqLim,mm))).^2) / sum(freqLim);
                if tmp == 0
                    D(ll,mm) = NaN;
                else
                    D(ll,mm) = tmp;
                end
            end
        end

        d = nanmedian(D);
        w = 1./d.^2;
        w = w/sum(w);
        w = repmat(w, [size(AllFramesFTrealign,1) 1]);

        out.off(ii,:) = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF==0),2);
        out.on(ii,:)  = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF==1),2);

        reject(:,ii) = zeros(1,size(AllFramesFTrealign,2));

    else

        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            reject(:,ii) = zeros(size(AllFramesFTrealign,2),1);
        end
        out.off(ii,:) = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)' & MRS_struct.out.reject(:,ii)==0), 2);
        out.on(ii,:)  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)' & MRS_struct.out.reject(:,ii)==0), 2);

    end

    out.diff(ii,:) = (out.on(ii,:) - out.off(ii,:))/2;

    if isfield(MRS_struct.spec,'AllFramesFT') % Catch case where unmodified version of GannetLoad might have been used.
        AllFramesFT = MRS_struct.spec.AllFramesFT;
        out.diff_noalign(ii,:) = (mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==1)),2) - mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==0)),2))/2;
    else
        warning('AllFramesFT hasn''t been stored, diff_noalign will be the same as diff.')
        out.diff_noalign(ii,:) = out.diff(ii,:);
    end
end

end

function MRS_struct = redoWaterRemoval(MRS_struct)
% Redo the water removal as per lines657 682 of Gannet Load
% Force 'vox 1'
vox = MRS_struct.p.Vox(1);
kk = 1;

if MRS_struct.p.water_removal
    for ii = 1:MRS_struct.p.numscans % Loop of batched files
        for jj = 1:length(MRS_struct.p.target)
            if jj == 1
                fprintf('\nFiltering out residual water signal...\n');
            end

            % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
            MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:).')), ...
                MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:)));

            MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:).')), ...
                MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:)));

            % MM (170703): Need to perform baseline correction on filtered data
            freqbounds = MRS_struct.spec.freq <= 10 & MRS_struct.spec.freq >= 9;
            baseMean_diff = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,freqbounds)));
            baseMean_diffnoalign = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,freqbounds)));

            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) - baseMean_diff;
            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) - baseMean_diffnoalign;
        end
    end
    
end
end