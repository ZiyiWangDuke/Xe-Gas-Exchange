% this file register BHUTE to high-resolution gas 
% 2016.7.14
tic;
direct = dir();
nfolder = length(direct)-2;
ratio = struct();
flag_patient = 0; % 1 for patient, 0 for healthy

for k = 1:1:nfolder
    folder = direct(k+2).name;
    ratio(k).name = folder;
    cd(folder);
    
    % load high resolution gas, lung mask and UTE data
    load('parameters');
    
    [pathstr,name,~] = fileparts(dixon_pfile);
    gas_file = [pathstr filesep() name '_gas_recon_old.nii'];
    if(exist(gas_file,'file') ~=2 )
        demo_dixon_recon_ventGas(dixon_pfile); % recon high resolution gas
    end
    gasVol = load_nii(gas_file);
    gasVol = abs(gasVol.img);
    gas_nom = gasVol/max(gasVol(:));
     
    [pathstr,name,~] = fileparts(bhute_pfile);
%     lung_mask = load_nii([pathstr filesep() name '_lungMask.nii']);
%     lung_mask = abs(lung_mask.img);
    lung_mask = load_nii('LungMask_grow.nii');
    lung_mask = abs(lung_mask.img);
    uteVol = load_nii([pathstr filesep() name '_bhute_recon.nii']);
    uteVol = abs(uteVol.img);
          
    % define registration transform
    [optimizer,metric] = imregconfig('multimodal');
    optimizer = registration.optimizer.RegularStepGradientDescent();
    optimizer.MaximumIterations = 20;
    transformtype = 'affine';
    % apply transformation
    [lung_reg,scaleinfo]  = imregister(lung_mask,gas_nom,transformtype,optimizer,metric,'DisplayOptimization',true);
    % apply the same transformation on UTE to obtain a registered UTE
    tform = imregtform(lung_mask,gas_nom,transformtype,optimizer,metric);
    ute_reg = imwarp(uteVol,tform,'OutputView',scaleinfo);
    % generate new binary mask from the registered lung mask
    lung_reg(lung_reg<0.5) = 0;
    lung_reg = boolean(lung_reg);
    
    % save the files
    lung_nii = make_nii(double(lung_reg));
    save_nii(lung_nii,'Reg_lungMask_grow.nii');
    ute_nii = make_nii(double(ute_reg));
    save_nii(ute_nii,'Reg_ute.nii');
    
    cd('../');
end
toc;