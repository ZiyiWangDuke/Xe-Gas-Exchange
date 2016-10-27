% this file calculates RBC:Gas map and conducts binning based on histogram
% Ratiomap based on pixel-wise division
%% read in
direct = dir();
nfolder = length(direct)-2;
cluster = struct();
clustermax = zeros(nfolder,6);

mmap= zeros(7,3);  
mmap(1,:) =[0 0 0]; % background
mmap(2,:) = [1 0 0]; % defect 1 
mmap(3,:) = [1 0.7143 0];% yellow 2 
mmap(4,:) = [0.4 0.7 0.4];% yellow - GREEN 3 
mmap(5,:) = [0 1 0];%  GREEN 4 
mmap(6,:) = [0 0.57 0.71]; % green 5 
mmap(7,:) = [0 0 1]; % high-intense 6 

uteflag = 1; %1 for overlap with UTE, 0 for not
imjflag = 1; %1 for imagej operation, 0 for not
histflag = 1 ; %1 for ploting histogram, 0 for not

mean_rbc = 0.2592;
std_rbc = 0.1038;%0.1037;
thre = [mean_rbc-2*std_rbc,mean_rbc-std_rbc,mean_rbc,mean_rbc+std_rbc,mean_rbc+2*std_rbc]; % threshold for rbc:gas

% figure;
sumpixel = [];
% lauch ImageJ, for montage
% montagedir = 'C:\Users\zw73\Desktop\ratipcolors\';
% stackdir = 'C:\Users\zw73\Desktop\stacks\RBC2gas_binning\';
montagedir = 'C:\Users\zw73\Desktop\ratiocolors_reg\RBC2gas\';
stackdir = 'C:\Users\zw73\Desktop\stacks_reg\RBC2gas\';
path_colormap = 'C:\Users\zw73\Desktop\colorhist_reg\RBC2gas\';
ImageJ;
for k = 1:1:nfolder
    folder = direct(k+2).name;
    cluster(k).  name = folder;
    cd(folder);
    % check if folder have been created
    folder_tif = 'binning_rbc';
    if (exist(folder_tif,'dir') == 0)
        mkdir(folder_tif);
    end
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    
    mask = load_nii('Reg_lungMask_growvent.nii');
    mask = boolean(mask.img);
    
    % load path data for data saving 
    load('parameters');
    TE90 = exp((TE90/1000000)/(2/1000))/exp((TE90/1000000)/0.05);
    flipoff = sin(flipoff/(2.5714*180)*pi)/sin(22/180*pi)*100;
    [pathstr,name,ext] = fileparts(bhute_pfile);
    
    % mask
    rbc = rbc.*mask;
    gas = gas .*mask;
    
    % initial ratio maps
    rbc2gas = zeros(size(mask));
    
    % small value rule out
    % 0.01 of the average value of the biggest 5% of the data is taken as a
    % threshold
    gmask = mask;
    
    % average ratio 
    rbc2gas(gmask) = rbc(gmask)./gas(gmask);
 
%%     binning for RBC
    binning = rbc2gas;
    binning = binning*TE90*flipoff;
    
    line = binning(:);
    line = line(line>0);
    if(histflag)
        f_ratiomap_colorhistogram(line,folder,0,path_colormap);
    end
    sumpixel = [sumpixel,line'];   
 
    binning(binning > thre(5)) = 7;
    binning(binning > thre(4) & binning <= thre(5)) = 6;
    binning(binning > thre(3) & binning <= thre(4)) = 5;
    binning(binning > thre(2) & binning <= thre(3)) = 4;
    binning(binning > thre(1) & binning <= thre(2)) = 3;
    binning(binning > 0 & binning <= thre(1)) = 2;
    binning(binning == 0) = 1;
    
    RBC_lo = (length(binning(binning==2))+length(binning(binning==3)))/length(line);
%    % save color binning
   if(imjflag)
        dimension = size(binning);
        for i= 1:dimension(3)
%            mkdir('binning_rbc');
           filename =['binning_rbc/IMG',int2str(i),'.tif'];
           im = binning(:,:,i);
           rgb = ind2rgb(im,mmap);   
%            rgb = binning(:,:,i);
           imwrite(rgb,filename);
        end
   end
    
   %% Quantitative analysis
   psum = sum(gmask(:));
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
    %% ImageJ section
    if(imjflag)
         if(uteflag)
             % clear the thoracic space
             
%              bhute = load_nii([pathstr filesep() name '_bhute_recon.nii']);
%              bhute = double(bhute.img);
             
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
            MIJ.run('Image Sequence...', ['open=',pwd,'\binning_rbc\IMG1.tif sort']);
            MIJ.run('8-bit Color', 'number=256');
            % merge channel
            ims = MIJ.getListImages;
            MIJ.run('Merge Channels...', ['c1=',char(ims(1)),' c4=',char(ims(2)),' create keep']);
            MIJ.run('RGB Color', 'slices keep');
        else
            MIJ.run('Image Sequence...', ['open=',pwd,'\binning_rbc\IMG1.tif sort']);
        end
        % save stack
        MIJ.run('NIfTI-1', ['save=',stackdir,'RBC2gas_',folder,'rr.nii']);
        % montage
        MIJ.run('Make Montage...', 'columns=8 rows=2 scale=1 first=18 last=64 increment=2 border=0 font=12');
        mon_im = MIJ.getImage('Montage');
        MIJ.run('Close All');
        imwrite(uint8(mon_im),[montagedir,'RBC2gas_',folder,'rr.png']);
    end
      
    cd('../');
end
MIJ.exit;
cluster(k+1).  name = 'All';
psum = sum(clustermax(:));
cluster(k+1).block1 = sum(clustermax(:,1))/psum;
cluster(k+1).block2 = sum(clustermax(:,2))/psum;
cluster(k+1).block3 = sum(clustermax(:,3))/psum;
cluster(k+1).block4 = sum(clustermax(:,4))/psum;
cluster(k+1).block5 = sum(clustermax(:,5))/psum;
cluster(k+1).block6 = sum(clustermax(:,6))/psum;
% %%
% % subplot(2,1,1);
% % figure;
% upbound = 0.8;
% sumpixel = [sumpixel,upbound];
% sumpixel(sumpixel>upbound) = upbound;
% nbins = 200;
% 
% % [disY,disX] = hist(sumpixel,nbins,'Normalization','probability');
% % hist(sumpixel,500,'Normalization','probability');
% map = brewermap(4,'Set1');
% histogram(sumpixel,nbins,'facecolor',map(1,:),'facealpha',.6,'Normalization','probability');
% xlim([0,0.8]);ylim([0,0.02]);
% ax = gca;
% ax.LineWidth = 3;
% ax.FontSize = 25;
% hold on;
% plot(thre,0.01*ones(1,5),'r*');
% 
% x = 0:1/nbins:1;
% a1 = 0.01529;
% b1 = 0.2604;
% c1 = 0.1481 ;
% normal = a1*exp(-((x-b1)/c1).^2);
% area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',4);
