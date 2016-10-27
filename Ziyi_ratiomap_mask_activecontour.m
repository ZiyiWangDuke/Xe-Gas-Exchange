%% This file generates a intact mask from original UTE mask with active contour to include the wedge area of the lung 
close all;clear;
%% read in
direct = dir();
nfolder = length(direct)-2;
ratio = struct();
%%
maxrate = 0;
for k = 1:1:nfolder
    folder = direct(k+2).name;
    ratio(k).name = folder;
    cd(folder);
    % load path data for data saving 
    load('parameters');
    [pathstr,name,ext] = fileparts(dixon_pfile);
    gas_old = load_nii([pathstr filesep() name '_gas_recon_old.nii']);
    gas_old = abs(gas_old.img);
    
    % load ute before registration and original mask
    [pathstr,name,ext] = fileparts(bhute_pfile);
    img_ute = load_nii([pathstr filesep() name '_bhute_recon.nii']);
    img_ute = abs(img_ute.img);
    mask_ute = load_nii([pathstr filesep() name '_lungMask.nii']);
    mask_ute = boolean(mask_ute.img);
    
%     img_ute = load_nii('Reg_ute.nii');
%     img_ute = abs(img_ute.img);
%     
%     if(exist('Reg_lungMask.nii','file') == 2)
%         mask_ute = load_nii('Reg_lungMask.nii');
%     else
%         [pathstr,name,ext] = fileparts(bhute_pfile);
%         mask_ute = load_nii([pathstr filesep() name '_lungMask.nii']);
%     end
%     mask_ute = boolean(mask_ute.img);
    
    %% Active contour for wedge recruitment 
    % output of the section: mask_grow
    % function: find lung slices in original mask
    mask_line = squeeze(sum(sum(mask_ute,2),1));
    index = find(mask_line>0);
    
    mask_grow = zeros(64,64,64);

    for kk = index(1):1:index(end)
        tag_ute = img_ute(:,:,kk);
        tag_msk = mask_ute(:,:,kk);
        
        %% this active contour still need to be improved
        bw_agre = activecontour(tag_ute,tag_msk,'Chan-Vese','SmoothFactor',0.8);
%         bw_cons = activecontour(tag_ute,tag_msk,'Chan-Vese','SmoothFactor',0.8);
%         
%         dif_second = bw_agre - bw_cons;
%         dif_second = sum(dif_second(:));
%         dif_primary = bw_cons - tag_msk;
%         dif_primary = sum(dif_primary(:));
%         sum_msk = sum(tag_msk(:));
%         
%         thre_primary = (sum_msk*0.25>120)*sum_msk*0.25 + (sum_msk*0.25<=120)*120;
%         thre_second = (sum_msk*0.1>40)*sum_msk*0.1 + (sum_msk*0.1<=40)*40;
%         
%         if dif_primary>thre_primary
%             mask_grow(:,:,kk) = tag_msk;
%         else if dif_second>thre_second;
%                 mask_grow(:,:,kk) = bw_cons;
%             else
%                 mask_grow(:,:,kk) = bw_agre;
%             end
%         end
        
        mask_grow(:,:,kk) = bw_agre;
%         if dif_primary>120
%             % original mask
%             mask_grow(:,:,kk) = tag_msk;
%         else 
%             mask_grow(:,:,kk) = (dif_second>50).*bw_cons+(dif_second<=50).*bw_agre;
%         end
    end
    mask_grow_nii = make_nii(abs(mask_grow));
    save_nii(mask_grow_nii,'LungMask_grow.nii');
    %% Ventilation mask
%     % derive mean and std from gas_old and generate threshold from them
%     [row,col,height] = size(gas_old);
%     reference_slice = mask_ute(:,:,round(height/2));
%     reference_gas = gas_old(:,:,round(height/2));
%     % Define the region of noise: the regions for noise calculaton are the left, right and bottom side
%     % of the mask
%     noise_region_mask = zeros(row,col);
%     
%     updown = sum(reference_slice,1);
%     nonzero = find(updown>0);
%     left = nonzero(1)-5;left = (left>0)*left;
%     right = nonzero(end)+5;right = (right>col)*col+(right<col)*right;
%     noise_region_mask(:,1:left) = 1;
%     noise_region_mask(:,right:end) = 1;
% 
%     leftright = sum(reference_slice,2);
%     nonzero = find(leftright>0);
%     bottom = nonzero(end)+5;bottom = (bottom>row)*row+(bottom<row)*bottom;
%     noise_region_mask(bottom:end,:) = 1;
%     noise_region_mask = boolean(noise_region_mask);
%     
%     % calculate ventilation threshold
%     noise_mean = mean(reference_gas(noise_region_mask));
%     noise_std = std(reference_gas(noise_region_mask));
%     ven_threshold = noise_mean+2*noise_std;
%     
%     mask_ven = (gas_old>ven_threshold);
%     mask_combined = mask_ven & mask_grow;
%        
%     mask_combined_nii = make_nii(abs(mask_combined));
%     save_nii(mask_combined_nii,'Reg_lungMask_vent.nii');
       
    % imopen and save ventilation mask
%     [xgrid, ygrid, zgrid] = meshgrid(-2:2);
%     ball = (sqrt(xgrid.^2 + ygrid.^2 + zgrid.^2) <= 1);
%     mask_ven = imopen(mask_ven,ball);
%     mask_ven_nii = make_nii(abs(mask_ven));
%     save_nii(mask_ven_nii,'Reg_lungMask_vent_pure.nii');
    
%     figure;
%     subplot(1,2,1);
%     imshow(reference_gas,[]);
%     title(name);
%     subplot(1,2,2);
%     imshowpair(mask_combined(:,:,32),mask_ute(:,:,32));
%     dif = abs(mask_combined - mask_ute);
%     rate = sum(dif(:))/sum(mask_ute(:));
%     title(['Percentage: ',num2str(rate)]);
%     maxrate = (maxrate>rate)*maxrate + (maxrate<rate)*rate;
    display([folder,' is completed,',num2str(k/nfolder),' is completed']);
    cd('../');
end
