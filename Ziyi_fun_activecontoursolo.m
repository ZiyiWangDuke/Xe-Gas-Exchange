% % registration test with subject, read in lung_mask, gasVol and uteVol;
% % close all;
% gasVol = abs(gasVol);
% lung_mask= abs(lung_mask);
% uteVol = abs(uteVol);
% lung_nom = lung_mask/max(lung_mask(:));
% gas_comp = gasVol/max(gasVol(:));
% % gas_comp = ones(size(gasVol))-gas_comp;
% [optimizer,metric] = imregconfig('multimodal');
% % optimizer.GrowthFactor = 1.001;
% % optimizer.InitialRadius = 0.0009;
% % optimizer.MaximumIterations = 500;
% % optimizer.Epsilon = 1.5e-06;'
% optimizer = registration.optimizer.RegularStepGradientDescent();
% optimizer.MaximumIterations = 20;
% transformtype = 'affine';
% [lung_reg,scaleinfo]  = imregister(lung_nom,gas_comp,transformtype,optimizer,metric,'DisplayOptimization',true);
% % apply the same transformation on UTE to obtain a registered UTE
% tform = imregtform(lung_nom,gas_comp,transformtype,optimizer,metric);
% ute_reg = imwarp(uteVol,tform,'OutputView',scaleinfo);
% % generate new binary mask from the old mask
% lung_reg(lung_reg<0.5) = 0;
% lung_reg = boolean(lung_reg);
% % figure;
% % imslice(uteVol);
% % figure;
% % imslice(ute_reg);
% % figure;
% % imslice(gasVol);
% %%
% figure;
% slice = 37;
% subplot(1,3,1);
% imagesc(uteVol(:,:,slice));
% axis image; colormap('gray');
% subplot(1,3,2);
% imagesc(ute_reg(:,:,slice));axis image;
% subplot(1,3,3);
% imshowpair(double(lung_reg(:,:,slice)),lung_nom(:,:,slice))
% % imagesc(gas_comp(:,:,slice));axis image;
% %%
% % ute_reg(ute_reg<0.5) = 0;
% % ute_reg = boolean(ute_reg);
% lung_nii = make_nii(double(lung_reg));
% save_nii(lung_nii,'lungMask_reg.nii');
% ute_nii = make_nii(double(ute_reg));
% save_nii(ute_nii,'ute_reg.nii');
%% active contour lung segmentation
img_ute = load_nii('Reg_ute.nii');
img_ute = abs(img_ute.img);
img_msk = load_nii('Reg_lungMask.nii');
img_msk = abs(img_msk.img);
img_vent = load_nii('Reg_lungMask_vent_pure.nii');
img_vent = abs(img_vent.img);

mask_line = squeeze(sum(sum(img_msk,2),1));
index = find(mask_line>0);

grow_msk = zeros(64,64,64);
filename = 'comparison.gif';

for k = index(1):1:index(end)
tag_ute = img_ute(:,:,k);
tag_msk = img_msk(:,:,k);
subplot(1,4,1);
imshowpair(tag_ute,tag_msk);
title('Original Mask');
bw_agre = activecontour(tag_ute,tag_msk,'Chan-Vese','SmoothFactor',0.1);
subplot(1,4,2);
imshowpair(tag_ute,bw_agre);
title('Intact Mask Aggresive');
subplot(1,4,3);
bw_cons = activecontour(tag_ute,tag_msk,'Chan-Vese','SmoothFactor',0.8);
dif_second = bw_agre - bw_cons;
dif_second = sum(dif_second(:));
dif_primary = bw_cons - tag_msk;
dif_primary = sum(dif_primary(:));
sum_msk = sum(tag_msk(:));
thre_primary = (sum_msk*0.25>120)*sum_msk*0.25 + (sum_msk*0.25<=120)*120;
thre_second = (sum_msk*0.1>40)*sum_msk*0.1 + (sum_msk*0.1<=40)*40;

imshowpair(tag_ute,bw_cons);
title(['Intact Mask Conservative',num2str(dif_second),'  ',num2str(thre_second)]);
subplot(1,4,4);


if dif_primary>thre_primary
    imshowpair(tag_ute,tag_msk);
    title('Intact Mask with Original');
else if dif_second>50;
        imshowpair(tag_ute,bw_cons);
        title('Intact Mask with Mild');
    else
        imshowpair(tag_ute,bw_agre);
        title('Intact Mask with Aggresive');
    end
end

pause(0.1);
% grow_msk(:,:,k) = bw;

% save into gif
% frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if k == index(1);
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
end
% grow_msk = make_nii(grow_msk);
% save_nii(grow_msk,'intactmask.nii');
% subplot(1,3,1);
% imshow(tag_ute,[]);
% subplot(1,3,2);
% bw = activecontour(tag_ute,tag_msk);
% imshow(bw,[]);
% subplot(1,3,3);
% imshowpair(tag_ute,bw);
%%
