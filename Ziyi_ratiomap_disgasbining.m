% this file calculates Barrier:Gas map and conducts binning based on histogram
% Ratiomap based on pixel-wise division
%% read in
% clear all;
direct = dir();
nfolder = length(direct)-2;
cluster = struct();
clustermax = zeros(nfolder,6);

% define color map for bar2gas (brown color from colorbrewer)
mmap= zeros(9,3);  
mmap(1,:) =[0 0 0]; % background
mmap(9,:) =[197,27,125]/255; % highest value
mmap(8,:) =[225,129,162]/255;
mmap(7,:) =[243,205,213]/255;
mmap(6,:) =[184,226,145]/255;
mmap(2,:) = [1 0 0]; % lowest value,red
mmap(3,:) = [1 0.7143 0];% yellow 2 
mmap(4,:) = [0.4 0.7 0.4];% yellow - GREEN 3 
mmap(5,:) = [0 1 0];%  GREEN 4 


uteflag = 1; %1 for overlap with UTE, 0 for not
imjflag = 1; %1 for imagej operation, 0 for not
histflag = 0; %1 for ploting histogram, 0 for not

mean = 0.4447;
std = 0.1419;
thre = [mean-2*std,mean-std,mean,mean+1*std,mean+2*std,mean+3*std,mean+4*std]; % 8 bins for map
% thre = [mean-1.5*std,mean,mean+1.5*std,mean+3*std,mean+4.5*std]; % threshold for bar:gas

% figure;
sumpixel = [];

% lauch ImageJ, for montage
montagedir = 'C:\Users\zw73\Desktop\ratipcolors\Dis2gas\';
stackdir = 'C:\Users\zw73\Desktop\stacks\dis2gas_binning\';
ImageJ;
for k = 1:1:nfolder
    folder = direct(k+2).name;
    cluster(k).name = folder;
    cd(folder);
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    
    % load path data for data saving ;
    load('parameters');
%     TE90 = exp((TE90/1000-0.932)/2);
    TE90 = exp((TE90/1000000)/(2/1000))/exp((TE90/1000000)/(0.0106));
    flipoff = sin(flipoff/(2.5714*180)*pi)/sin(22/180*pi)*100;
    [pathstr,name,ext] = fileparts(bhute_pfile);
    
    % mask
    dissolved = abs(dissolved).*mask;
    gas = gas .*mask;
    
    % initial ratio maps
    dis2gas = zeros(size(mask));
    
    % small value rule out
    % 0.01 of the average value of the biggest 5% of the data is taken as a
    % threshold
%     gmask = mask;
    [ gmask ] = f_smallruleout( gas,0.1,0.05);

    % average ratio 
    dis2gas(gmask) = dissolved(gmask)./gas(gmask);
    dis2gas(dis2gas<0) = 1e-5;
    
%     ratio(k).rbc2bar = sum(rbc2bar(:))/maskall;
%     ratio(k).rbc2gas = sum(rbc2gas(:))/maskall;
%     ratio(k).bar2gas = sum(bar2gas(:))/maskall;
    
    

%     subplot(4,4,k);
%     [disY,disX] = hist(line,100);
%     hist(line,100);
%     xlim([0,1]);
%     ylim([0,1500]);
%     title([folder,' , ',num2str(sum(gmask(:))),' , ',num2str(sum(bar2gas(:))/sum(gmask(:)))]);
    
    % binning for barrier
    binning = dis2gas;
    binning = binning*TE90*flipoff; % correction for TE and gas flip off
    
    line = binning(:);
    line = line(line>0);
    if(histflag)
        f_ratiomap_colorhistogram(line,folder,1);
    end
    sumpixel = [sumpixel,line'];

    binning(binning > thre(7)) = 9;
    binning(binning > thre(6) & binning <= thre(7)) = 8;
    binning(binning > thre(5) & binning <= thre(6)) = 7;
    binning(binning > thre(4) & binning <= thre(5)) = 6;
    binning(binning > thre(3) & binning <= thre(4)) = 5;
    binning(binning > thre(2) & binning <= thre(3)) = 4;
    binning(binning > thre(1) & binning <= thre(2)) = 3;
    binning(binning > 0 & binning <= thre(1)) = 2;
    binning(binning == 0) = 1;
    
%     bar_hi = length(binning(binning==9|binning==8|binning==7))/length(line);
    
    %save color binning
    if(imjflag)
       dimension = size(binning);
        for i= 1:dimension(3)
%            mkdir('binning_bar');
           filename =['binning_dis/IMG',int2str(i),'.tif'];
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
   a = binning == 8;
   clustermax(k,7) = sum(a(:));
   cluster(k).block7 = clustermax(k,7)/psum;
   a = binning == 9;
   clustermax(k,8) = sum(a(:));
   cluster(k).block8 = clustermax(k,8)/psum;
   %% ImageJ section
   if(imjflag)
        if(uteflag)
             % clear the thoracic space
             bhute = load_nii([pathstr filesep() name '_bhute_recon.nii']);
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
            MIJ.run('Image Sequence...', ['open=',pwd,'\binning_dis\IMG1.tif sort']);
            MIJ.run('8-bit Color', 'number=256');
            % merge channel
            ims = MIJ.getListImages;
            MIJ.run('Merge Channels...', ['c1=',char(ims(1)),' c4=',char(ims(2)),' create keep']);
            MIJ.run('RGB Color', 'slices keep');
        else
            MIJ.run('Image Sequence...', ['open=',pwd,'\binning_dis\IMG1.tif sort']);
        end
        % save stack
        MIJ.run('NIfTI-1', ['save=',stackdir,'Dis2gas_',folder,'.nii']);
        % montage
        MIJ.run('Make Montage...', 'columns=8 rows=2 scale=1 first=18 last=64 increment=2 border=0 font=12');
        mon_im = MIJ.getImage('Montage');
        MIJ.run('Close All');
        imwrite(uint8(mon_im),[montagedir,'Dis2gas_',folder,'.png']);
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
cluster(k+1).block7 = sum(clustermax(:,7))/psum;
cluster(k+1).block8 = sum(clustermax(:,8))/psum;
 %%
% subplot(2,1,1);
% map = brewermap(4,'Set2');
%  histogram(sumpixel,300,'facecolor',map(1,:),'facealpha',.6,'Normalization','probability');
%  xlim([0,1.8]);ylim([0,0.03]);
%  ax = gca;
%  ax.LineWidth = 3;
%  ax.FontSize = 25;
%  hold on;
%  plot(thre,0.01*ones(1,7),'r*');
