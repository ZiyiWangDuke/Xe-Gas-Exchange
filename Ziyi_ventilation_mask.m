% input gasVol;
% input lung_mask;
%% determine ventilation threshold;
factor = 3; % the higher, the smaller excluded area
absgas = abs(gasVol);
absgas = absgas.*lung_mask; % constrain with thoracic cavity
line = absgas(:);
line = line(line>0);
histogram(line,200);
thre_98 = prctile(line,98);
% line_95 = line(line<thre_95);
% [x,y] = hist(line_95,200);
% mean_ven = mean(line_95);
% std_ven = std(line_95);
% [fitresult, gof] = f_gaussianfit(y, x);
% peak_x = fitresult.b1;
% vent_thre = peak_x/2;
vent_thre = thre_98/factor;
%% generate ventilation mask;
gas_mask = absgas>vent_thre;
mask_vent  = gas_mask & lung_mask;
nii = make_nii(abs(mask_vent));
save_nii(nii,'ventilation_mask.nii');