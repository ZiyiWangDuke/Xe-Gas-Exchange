function demo_dixon_recon_ventGas(varargin)
% reconstruct gas with the same parameter 
if(nargin < 1 | ~exist(varargin{1}))
    disp('Select Dixon pfile');
    dixon_pfile = filepath('/home/scott/');
else
    dixon_pfile = varargin{1};
end

% Prepare to override values
pfileOverride = GE.Pfile.Pfile();

clear startDissolved startGas;

% For VERY OLD scan format (46)
% output_image_sizeg = 64*[1 1 1];
% output_image_sized = 64*[1 1 1];
% overgrid_factor = 3;
% gasKernel.sharpness = 0.14;
% gasKernel.extent = 9*gasKernel.sharpness;
% dissolvedKernel.sharpness = 0.14;
% dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
% verbose = 1;
% nPipeIter = 25;
% pfileOverride.rdb.rdb_hdr_user1  = 0.252; % pw_gxwa
% pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
% pfileOverride.rdb.rdb_hdr_user44 = 1.024; % pw_gxw/1000
% pfileOverride.rdb.rdb_hdr_user22 = 0.1325; %toff
% % pfileOverride.rdb.rdb_hdr_user22 = 0.0625; %toff
% pfileOverride.rdb.rdb_hdr_user32 = 1; % Archimedian spiral
% downstream_magFrames = 50;
% rmBline = 0;
% rmFirstGas = 1;
% rmFirstDis = 1;

% For old format (46A,54,55,39B,49C,49D)
% output_image_sizeg = 64*[1 1 1];
% % output_image_sized = 64*[1 1 1];
% overgrid_factor = 3;
% gasKernel.sharpness = 0.32;
% gasKernel.extent = 9*gasKernel.sharpness;
% % dissolvedKernel.sharpness = 0.14;%0.14;
% % dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
% verbose = 1;
% nPipeIter = 25;
% pfileOverride.rdb.rdb_hdr_user1  = 0.252; % pw_gxwa
% pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
% pfileOverride.rdb.rdb_hdr_user44 = 1.024; % pw_gxw/1000
% pfileOverride.rdb.rdb_hdr_user22 = 0.1325; %toff
% pfileOverride.rdb.rdb_hdr_user32 = 0;  % Golden Means
% % downstream_magFrames = 50;
% rmBline = 1;
% rmFirstGas = 0;
% rmFirstDis = 0;

% For 64
% output_image_sizeg = 64*[1 1 1];
% output_image_sized = 64*[1 1 1];
% overgrid_factor = 3;
% gasKernel.sharpness = 0.32;
% gasKernel.extent = 9*gasKernel.sharpness;
% dissolvedKernel.sharpness = 0.14;
% dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
% verbose = 1;
% nPipeIter = 25;
% pfileOverride.rdb.rdb_hdr_user1  = 0.508; % pw_gxwa
% pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
% pfileOverride.rdb.rdb_hdr_user44 = 2.048; % pw_gxw/1000
% pfileOverride.rdb.rdb_hdr_user22 = 0.125; %toff
% % pfileOverride.rdb.rdb_hdr_user32 = 0;  % Golden Means
% downstream_magFrames = 50;
% rmBline = 0;
% rmFirstGas = 1;
% rmFirstDis = 0;
% startDissolved = 2;
% startGas = 1;

% For most recent format (65)
output_image_sizeg = 64*[1 1 1];
output_image_sized = 64*[1 1 1];
overgrid_factor = 3;
gasKernel.sharpness = 0.14;
gasKernel.extent = 9*gasKernel.sharpness;
dissolvedKernel.sharpness = 0.14;
dissolvedKernel.extent = 9*dissolvedKernel.sharpness;
verbose = 1;
nPipeIter = 25;
pfileOverride.rdb.rdb_hdr_user1  = 0.512; % pw_gxwa
pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
pfileOverride.rdb.rdb_hdr_user44 = 1.536; % pw_gxw/1000
pfileOverride.rdb.rdb_hdr_user22 = 0.125; %toff
% pfileOverride.rdb.rdb_hdr_user32 = 0;  % Golden Means
downstream_magFrames = 50;
rmBline = 0;
rmFirstGas = 1;
rmFirstDis = 0;

deltaf_gas = 0;

% Gradient delays
delays.x_delay = 0.000;
delays.y_delay = 0.00;
delays.z_delay = 0.000;

%Optional parameters
% revision_override = [];  %Optional override if it can't be automatically read from the pfile

%% Read Raw Pfile and process pfile
pfile = GE.Pfile.read(dixon_pfile);

% Convert from Pfile format
pfile = convertLegacyPfile(pfile);

% Override header values (optional)
displayPfileHeaderInfo(pfile);
pfile = overridePfile(pfile, pfileOverride);
if(verbose)
    % Display key header info
    displayPfileHeaderInfo(pfile);
end

% Check for overranging
MRI.DataProcessing.checkForOverranging(pfile);

% Remove baselines
if(rmBline)
    pfile = MRI.DataProcessing.removeBaselineViews(pfile);
end

%% Split pfile into 2: disolved and gas
if(~                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   exist('startDissolved') & ~exist('startGas'))
    if(pfile.rdb.rdb_hdr_ps_mps_freq == 176604450)
        % dissolved first
%         startDissolved = 2;
        startGas = 1;
    else
        % gas first
%         startDissolved = 1;
        startGas = 2;
    end
end

gas_pfile = pfile;
gas_pfile.data = gas_pfile.data(:,startGas:2:end);
if(rmFirstGas)
    gas_pfile.data = gas_pfile.data(:,2:end);
end
     gas_pfile.rdb.rdb_hdr_user20 = size(gas_pfile.data,2);
%% Calculate trajectories (will be same for both pfiles)
% Calculate trajectory for a single radial ray
radialDistanceg = MRI.Trajectories.Centric.Radial.calcRadialRay(gas_pfile, delays, output_image_sizeg);

% % 	% Only keep data during gradients on
% [junk, dissolved_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistanced, dissolved_pfile);
% [radialDistance, gas_pfile] = MRI.DataProcessing.removeNonReadoutSamples(radialDistanceg, gas_pfile);

% Distribute rays onto 3d sphere
trajg = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistanceg, gas_pfile);

% Undo loopfactor
[trajg, gas_pfile] = MRI.DataProcessing.undoloopfactor(trajg, gas_pfile);

% Demodulate signal
bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
dwell_time = 1/(2*bw*1000);                                             % Time between each sample
dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
nPts = pfile.rdb.rdb_hdr_frame_size;           % Number of sample points per frame/ray

tMatg = repmat(dwell_time*(1:nPts)',[1 size(gas_pfile.data,2)]);

gas_pfile.data = gas_pfile.data.*exp(1i*2*pi*deltaf_gas*tMatg);

% Calculate Maximum volume size for Nyquist

MRI.DataProcessing.calculateNyquistMatrixSize(radialDistanceg, gas_pfile);


% Vectorize data and traj for recon

[trajg, gas_pfile] = MRI.DataProcessing.vectorizeDataAndTraj(trajg, gas_pfile);

% Enforce Nyquist limits

[trajg, gas_pfile] = MRI.DataProcessing.enforceNyquistBounds(trajg, gas_pfile);

%% Reconstruct Gas data
% Choose kernel
gasKernelObj = Recon.SysModel.Kernel.Gaussian(gasKernel.sharpness, gasKernel.extent, verbose);

% Choose Proximity object
gasProxObj = Recon.SysModel.Proximity.L2Proximity(gasKernelObj, verbose);
clear gasKernelObj;

% Create System model
gasSystemObj = Recon.SysModel.MatrixSystemModel(trajg, overgrid_factor, ...
    output_image_sizeg, gasProxObj, verbose);

% Choose density compensation function (DCF)
gasDcfObj = Recon.DCF.Iterative(gasSystemObj, nPipeIter, verbose);

% Choose Reconstruction Model
gasReconObj = Recon.ReconModel.LSQGridded(gasSystemObj, gasDcfObj, verbose);
gasReconObj.crop = 1;
gasReconObj.deapodize = 1;

gasVol = gasReconObj.reconstruct(gas_pfile.data, trajg);


%% save the result
[pathstr,name,ext] = fileparts(dixon_pfile);
save([pathstr filesep() name '_gas_recon.mat'],'gasVol');

[pathstr,name,ext] = fileparts(dixon_pfile);
niiname = [pathstr filesep() name '_gas_recon.nii'];
nii = make_nii(abs(gasVol));
save_nii(nii,niiname);

end