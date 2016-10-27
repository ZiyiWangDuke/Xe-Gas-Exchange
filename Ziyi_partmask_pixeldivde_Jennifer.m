     % this file calculates RBC:Bar, RBC:Gas, Bar:Gas maps and the average
% ratios. The three inputs are used for comparison
% Ratio maps are calculated by divide them on a pixel base.
% close all;clear all;
%% read in
direct = dir();
nfolder = length(direct)-2;
ratio = struct();
num_msk = 8;
msk_dir = 'Z:\shared\jennifer\Final_partialmasks_092816\IPF_gasexchange';
% msk_dir = 'Z:\shared\jennifer\Final_partialmasks_092816\healthy_vent';
for k = 1:1:nfolder
    folder = direct(k+2).name;
    ratio(k).name = folder;
    cd(folder);
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    % load path data for data saving 
    load('parameters');

    TE90 = exp((TE90/1000000)/(2/1000))/exp((TE90/1000000)/0.05);
    flipoff = sin(flipoff/(2.5714*180)*pi)/sin(22/180*pi)*100;
    
    [pathstr,name,ext] = fileparts(bhute_pfile);
    % mask direction
    msk_sub_dir = [msk_dir,'\',strrep(folder, 'S', 's')];
    direct_msk = dir(msk_sub_dir);
    %% for ventilation 
    [pathstr,name,ext] = fileparts(dixon_pfile);
    gas = load_nii([pathstr filesep() name '_gas_recon_old.nii']);
    gas = gas.img;
    maskall = load_nii('Reg_lungMask_grow.nii');
    maskall = boolean(maskall.img);
    
    gas = gas.*maskall;
    
    % binning for gas
    gasline = gas(gas>0);
    gasline = sort(gasline,'ascend');
    sumgas = cumsum(gasline./sum(gasline));
    % find the 99% sumup
    difarray = abs(sumgas-0.99);
    index = find(difarray == min(difarray));
    threshold = (sumgas(index(1)+1) - sumgas(index(1)))*sum(gasline);
    % gas normalization
    gasline(gasline>threshold) = threshold;
    gasline = gasline/threshold;
    gas(gas>threshold) = threshold;
    gas = gas/threshold;
    
    %%
    for m = 1:1:num_msk
        msk_name = direct_msk(m+2).name;
        jf_mask = load_nii([msk_sub_dir,'\',msk_name]);
        mask = boolean(jf_mask.img);

        rbc2bar = zeros(size(mask));
        
        mask(barrier<1) = 0;
%         rbc2bar(mask) = rbc(mask)./barrier(mask); 
%         rbc2bar(mask) = barrier(mask)*TE90*flipoff./gas(mask); 
        rbc2bar(mask) = gas(mask);
        rbc2bar(rbc2bar == inf) = 0;
        rbc2bar(isnan(rbc2bar)) = 0;
        
        finalmask = (rbc2bar>0);
        maskall = sum(finalmask(:));
        switch m
            case 1
                ratio(k).wholelung = sum(rbc2bar(:))/maskall;
            case 2
                ratio(k).anterior = sum(rbc2bar(:))/maskall;
            case 3
                ratio(k).apex = sum(rbc2bar(:))/maskall;
            case 4
                ratio(k).base = sum(rbc2bar(:))/maskall;
            case 5
                ratio(k).central = sum(rbc2bar(:))/maskall;
            case 6
                ratio(k).final = sum(rbc2bar(:))/maskall;
            case 7
                ratio(k).peripheral = sum(rbc2bar(:))/maskall;
            case 8
                ratio(k).posterior = sum(rbc2bar(:))/maskall;
            otherwise
                display('There is an error');
        end
%         ratio(k).rbc2gas = sum(rbc2gas(:))/maskall;
%         ratio(k).bar2gas = sum(bar2gas(:))/maskall;

%         nii = make_nii(bar2gas);
%         save_nii(nii,[pathstr filesep() 'rmp_bar2gas.nii']);
%         nii = make_nii(rbc2gas);
%         save_nii(nii,[pathstr filesep() 'rmp_rbc2gas.nii']);
%         nii = make_nii(rbc2bar);
%         save_nii(nii,[pathstr filesep() 'rmp_rbc2bar.nii']);
%         nii = make_nii(dis2gas);
%         save_nii(nii,[pathstr filesep() 'dis2gas.nii']);
    end
    cd('../');
end
