% this file calculates RBC:Gas map and conducts binning based on histogram
% Ratiomap based on pixel-wise division
%% set output path
path_root = 'C:\Git\Data_component\';
cohortfolder = 'IPF';

pathcom = [path_root,cohortfolder];
if ~exist(pathcom,'dir') 
    mkdir(pathcom);
end

direct = dir();
nfolder = length(direct)-2;

for k = 1:1:nfolder
    subfolder = direct(k+2).name;
    cd(subfolder);
    % check if folder have been created
    path_sub = [pathcom,'\',subfolder];
    if (exist(path_sub,'dir') == 0)
        mkdir(path_sub);
    end
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    
    niiname = [path_sub,'\',subfolder,'_dissolved.nii'];
    nii = make_nii(abs(dissolved));
    save_nii(nii,niiname);
    
    niiname = [path_sub,'\',subfolder,'_barrier.nii'];
    nii = make_nii(abs(barrier));
    save_nii(nii,niiname);
    
    niiname = [path_sub,'\',subfolder,'_rbc.nii'];
    nii = make_nii(abs(rbc));
    save_nii(nii,niiname);
    
    niiname = [path_sub,'\',subfolder,'_gas_hSNR.nii'];
    nii = make_nii(abs(gas));
    save_nii(nii,niiname);
    
    load('parameters');
    [pathstr,name,ext] = fileparts(dixon_pfile);
    load([pathstr filesep() name '_gas_recon_old.mat']);
    
    niiname = [path_sub,'\',subfolder,'_gas_lSNR.nii'];
    nii = make_nii(abs(gasVol));
    save_nii(nii,niiname);
    
    cd('../');
end

