% this file normalizas the gas recon file and conducts binning based on histogram
% Ratiomap based on pixel-wise division
%% read in
% clear;
direct = dir();
nfolder = length(direct)-2;
cluster = struct();
clustermax = zeros(nfolder,6);

% define color map for ventilation 
mmap= zeros(7,3);  
mmap(1,:) =[0 0 0]; % background
mmap(2,:) = [1 0 0]; % defect 1 
mmap(3,:) = [1 0.7143 0];% yellow 2 
mmap(4,:) = [0.4 0.7 0.4];% yellow - GREEN 3 
mmap(5,:) = [0 1 0];%  GREEN 4 
mmap(6,:) = [0 0.57 0.71]; % green 5 
mmap(7,:) = [0 0 1]; % high-intense 6 

uteflag = 1; %1 for overlap with UTE, 0 for not
imjflag = 1; %1 for imagej opeclustern, 0 for not
histflag = 1 ; %1 for ploting histogram, 0 for not

mean_vent = 0.5098;
std_vent = 0.1899;
thre = [mean_vent-2*std_vent,mean_vent-std_vent,mean_vent,mean_vent+std_vent,mean_vent+2*std_vent]; % threshold for bar:gas

% lauch ImageJ, for montage
montagedir = 'C:\Users\zw73\Desktop\ratiocolors_reg\ventilations\';
stackdir = 'C:\Users\zw73\Desktop\stacks_reg\ventilations\';
path_colormap = 'C:\Users\zw73\Desktop\colorhist_reg\ventilations\';
ImageJ;
sumpixel = [];
for k = 1:1:2%nfolder
    folder = direct(k+2).name;
    cluster(k).name = folder;
    cd(folder);
    % check if folder have been created
    folder_tif = 'ventilation';
    if (exist(folder_tif,'dir') == 0)
        mkdir(folder_tif);
    end
    
    % read in and preprocess data
    load('parameters');
    [pathstr,name,ext] = fileparts(dixon_pfile);
    gas = load_nii([pathstr filesep() name '_gas_recon_old.nii']);
    gas = gas.img;
    mask = load_nii('Reg_lungMask_grow.nii');
    mask = boolean(mask.img);
    
    gas = gas.*mask;
    
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

    if(histflag)
        f_ratiomap_colorhistogram(gasline,folder,2,path_colormap);
    end
    %% collect all the voxels
    sumpixel = [sumpixel;gasline];
    %%
    binning = gas;
    binning(binning > thre(5)) = 7;
    binning(binning > thre(4) & binning <= thre(5)) = 6;
    binning(binning > thre(3) & binning <= thre(4)) = 5;
    binning(binning > thre(2) & binning <= thre(3)) = 4;
    binning(binning > thre(1) & binning <= thre(2)) = 3;
    binning(binning > 0 & binning <= thre(1)) = 2;
    binning(binning == 0) = 1;
    
    %% save ventilation mask from the ventilation map
    mask_vent = binning;
    mask_vent(mask_vent<3) = 0;
    mask_vent(mask_vent>0)= 1;
    save_nii(make_nii(mask_vent),'Reg_lungMask_growvent.nii');
    %%
    %save color binning
    if(imjflag)
       dimension = size(binning);
       rmdir('ventilation','s');
       mkdir('ventilation');
        for i= 1:dimension(3)
%            mkdir('binning_bar');
           filename =['ventilation/IMG',int2str(i),'.tif'];
           im = binning(:,:,i);
           rgb = ind2rgb(im,mmap);   
           imwrite(rgb,filename);
        end
    end
    %% Quantitative analysis
   psum = sum(mask(:));
   a = binning == 2;
   clustermax(k,1) = sum(a(:));
   cluster(k).block1 = clustermax(k,1)/psum;
   a = binning == 3;
   clustermax(k,2) = sum(a(:));
   cluster(k).block2 = clustermax(k,2)/psum;
   a = binning == 4;
   clustermax(k,3) = sum(a(:));
   cluster(k).block3 = clustermax(k,3)/psum;
   a = binning == 5;
   clustermax(k,4) = sum(a(:));
   cluster(k).block4 = clustermax(k,4)/psum;
   a = binning == 6;
   clustermax(k,5) = sum(a(:));
   cluster(k).block5 = clustermax(k,5)/psum;
   a = binning == 7;
   clustermax(k,6) = sum(a(:));
   cluster(k).block6 = clustermax(k,6)/psum;
   cluster(k).ventilation = mean(gasline);
   %% ImageJ section
   if(imjflag)
        if(uteflag)
             % clear the thoracic space
             [pathstr,name,ext] = fileparts(bhute_pfile);
%              bhute = load_nii([pathstr filesep() name '_bhute_recon.nii']);
             bhute = load_nii('Reg_ute.nii');
             bhute = double(bhute.img);
             mbhute = bhute.*~mask;
             niiname = [pathstr filesep() name '_bhute_clean.nii'];
             nii = make_nii(mbhute);
             save_nii(nii,niiname);
            % open BHUTE image
            fname = f_findnamewithpart( 'bhute_clean.nii' ); % find ute recon file
            MIJ.run('Open...',['path=',pwd,'\',fname]);
            MIJ.run('Rotate 90 Degrees Right');
            MIJ.run('Flip Horizontally', 'stack');
            MIJ.run('Set Slice...','slice=35');
            MIJ.run('Enhance Contrast', 'saturated=0.35');  
            MIJ.run('RGB Color');
            MIJ.run('8-bit Color', 'number=256');
            % input the binning map
            MIJ.run('Image Sequence...', ['open=',pwd,'\ventilation\IMG1.tif sort']);
            MIJ.run('8-bit Color', 'number=256');
            % merge channel
            ims = MIJ.getListImages;
            MIJ.run('Merge Channels...', ['c1=',char(ims(1)),' c4=',char(ims(2)),' create keep']);
            MIJ.run('RGB Color', 'slices keep');
        else
            MIJ.run('Image Sequence...', ['open=',pwd,'\ventilation\IMG1.tif sort']);
        end
        % save stack
        MIJ.run('NIfTI-1', ['save=',stackdir,'ven_',folder,'.nii']);
        % montage
        MIJ.run('Make Montage...', 'columns=8 rows=2 scale=1 first=18 last=64 increment=2 border=0 font=12');
        mon_im = MIJ.getImage('Montage');
        MIJ.run('Close All');
        imwrite(uint8(mon_im),[montagedir,'ven_',folder,'.png']);
   end
   cd('../');
end
MIJ.exit;
cluster(k+1).name = 'All';
psum = sum(clustermax(:));
cluster(k+1).block1 = sum(clustermax(:,1))/psum;
cluster(k+1).block2 = sum(clustermax(:,2))/psum;
cluster(k+1).block3 = sum(clustermax(:,3))/psum;
cluster(k+1).block4 = sum(clustermax(:,4))/psum;
cluster(k+1).block5 = sum(clustermax(:,5))/psum;
cluster(k+1).block6 = sum(clustermax(:,6))/psum;
% subplot(2,1,1);
%%
% figure;
% map = brewermap(4,'Blues');
% histogram(sumpixel,150,'facecolor',map(4,:),'facealpha',.6,'Normalization','probability');
% xlim([0,1]);ylim([0,0.02])
% ax = gca;
% ax.LineWidth = 3;
% ax.FontSize = 25;
% hold on;
% plot(thre,0.01*ones(1,5),'r*');
% hold on;
% 
% x = 0:1/150:1;
% a1 = 0.01359;
% b1 = 0.4939;
% c1 = 0.278 ;
% normal = a1*exp(-((x-b1)/c1).^2);
% area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',4);