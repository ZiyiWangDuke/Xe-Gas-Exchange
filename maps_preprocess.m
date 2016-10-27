function [ gas,dissolvedVol,rbc_vol,barrier_vol,lungMask,mrbc2bar ] = maps_preprocess(parafile)
%Packed the pre-process part before calculating map  
load(parafile);
mrbc2bar = meanRbc2barrier;
RBC_barrier_spect = meanRbc2barrier;
% varargin{1} = meanRbc2barrier;
% varargin{2} = bhute_pfile;
% varargin{3} = dixon_pfile;

% if(nargin < 3 | ~isnumeric(varargin{1}) | ~exist(varargin{2}) | ~exist(varargin{3}))
%     RBC_barrier_spect = input('What was the RBC:barrier ratio from spectroscopy?:');
%     disp('Select BHUTE pfile');
%     bhute_pfile = filepath('/home/scott/Desktop/');
%     [pathstr,name,ext] = fileparts(bhute_pfile);
%     disp('Select Dixon pfile');
%     dixon_pfile = filepath(pathstr);
% else
%     RBC_barrier_spect = varargin{1};
%     bhute_pfile = varargin{2};
%     dixon_pfile = varargin{3};
% end 

[pathstr,name,ext] = fileparts(bhute_pfile);
lungMask = load_nii([pathstr filesep() name '_lungMask.nii']);
lungMask = boolean(lungMask.img);
[pathstr,name,ext] = fileparts(dixon_pfile); 
load([pathstr filesep() name '_gas_recon.mat']);
load([pathstr filesep() name '_dissolved_recon.mat']);   
enforceInLungOnly = 1;

%% Finds the phase offset to give the correct spectroscopy ratio
desired_angle = atan2(RBC_barrier_spect,1);
% Ziyi: Here may need RBC_gas phase offset to give the correct spect ratio
if(enforceInLungOnly)
    netVec = sum(dissolvedVol(lungMask));
else
    netVec = sum(dissolvedVol(:));
end
current_angle_complete = atan2(imag(netVec),real(netVec));
delta_angle_complete = desired_angle-current_angle_complete;

% % Check that angle gives correct RBC:barrier ratio
rotVol = dissolvedVol.*exp(1i*delta_angle_complete);
if(enforceInLungOnly)
    finalVec = sum(rotVol(lungMask));
else
    finalVec = sum(rotVol(:));
end
rbc_barrier_ratio_check1 = imag(finalVec)/real(finalVec);

%% Correct for B0 inhomogeneities
% Itterate until mean phase is zero
gas2 = gasVol;
iterCount = 0;
meanphase = inf;
while((abs(meanphase) > 1E-7))
    if(iterCount > 20)
        warning('Could not get zero mean phase within 10 iterations...');
    end
    if(iterCount > 100)
        error('Could not get zero mean phase within 10 iterations...');
    end
    iterCount = iterCount + 1;
    diffphase = angle(gas2);
    meanphase = mean(diffphase(lungMask(:)));
    gas2 = gas2*exp(-1i*meanphase);
end
diffphase = angle(gas2);

% % remove B0 inhomogeneities
dissolvedVol = dissolvedVol.*exp(1i*(delta_angle_complete-diffphase));

if(enforceInLungOnly)
    finalVec2 = sum(dissolvedVol(lungMask));
else
    finalVec2 = sum(dissolvedVol(:));
end
rbc_barrier_ratio_check2 = imag(finalVec2)/real(finalVec2);

rbc_vol = imag(dissolvedVol);
barrier_vol = real(dissolvedVol);

% Force positiveness and scale 
rbc_vol(rbc_vol<0) = 1e-5;
barrier_vol(barrier_vol<0) = 1e-5;
% rbc_vol = abs(rbc_vol);
% barrier_vol = abs(barrier_vol);
% rbc_vol(rbc_vol<0) = 0;
% barrier_vol(barrier_vol<0) = 0;

if(enforceInLungOnly)
    rbc_barrier_ratio_check3 = sum(rbc_vol(lungMask))/sum(barrier_vol(lungMask));
else
    rbc_barrier_ratio_check3 = sum(rbc_vol(:))/sum(barrier_vol(:));
end

gas = abs(gas2);
end

