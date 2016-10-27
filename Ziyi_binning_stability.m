% this file calculates Barrier:Gas map and conducts binning based on histogram
% Ratiomap based on pixel-wise division
%% read in
% clear all;
direct = dir();
nfolder = length(direct)-2;
cluster = struct();
clustermax = zeros(nfolder,6);

% TEcorrect for IPF patients
% {'Sub46';'Sub46A';'Sub46B';'Sub46C';'Sub47';'Sub50';'Sub54';'Sub54A';'Sub56';'Sub67';'Sub67A';'Sub68';'Sub68A';'Sub69';'Sub69A';'Sub69B';'Sub73';'Sub76';'Sub79';'Sub80'}
TEc = [1030,855,738,900,740,912,816,1000,740,1300,750,1200,900,1300,900,890,849,900,914,844,1300];
flipoff = [1,1.44,1.5,1.44,1,1,1.45,1.44,1.5,1.5000,1.4400,1.4400,1.4400,1.4380,1.4400,1.5060,1.44,1.44,1.44,1.4370,1.4620];

% TEcorrect for healthy subjects (include 003002)
% for {'Sub002003C';'Sub003';'Sub003002_health';'Sub003003_health';'Sub39B_health';'Sub48';'Sub49C';'Sub49D';'Sub55_health';'Sub64_health';'Sub65_health'}
% TEc = [985,947,819,1044,1027,1011,1005,1007,1070,1000,971];
% flipoff = [1.5,1.4400,1.44,1.4400,1.3600,1.3700,1,1.51,1.5000,1.5000,1.4400];

% TEcorrect for healthy subjects (exclude 003002)
% for {'Sub002003C';'Sub003';'Sub003003_health';'Sub39B_health';'Sub48';'Sub49C';'Sub49D';'Sub55_health';'Sub64_health';'Sub65_health'}
% TEc = [985,947,1044,1027,1011,1005,1007,1070,1000,971];
% flipoff = [1.5,1.4400,1.4400,1.3600,1.3700,1,1.51,1.5000,1.5000,1.44];

% TEcorrect for other disease
% TEc = [1000,910,794,914,815];
% flipoff = [1.5000,1.4270,1.4400,1.4400,1.4400];

TEc = exp((TEc/1000-0.932)/2);
% flipoff = sin(flipoff/180*pi)/sin(1.44/180*pi);
flipoff = sin(flipoff/(2.5714*180)*pi)/sin(22/180*pi)*100;
% define color map for bar2gas (brown color from colorbrewer)
mmap= zeros(7,3);  
mmap(1,:) =[0 0 0];
mmap(7,:) =[197,27,125]/255;
mmap(6,:) =[233,163,201]/255;
mmap(5,:) =[253,224,239]/255;
mmap(4,:) =[230,245,208]/255;
mmap(3,:) =[161,215,106]/255;
mmap(2,:) =[77,146,33]/255;

mean = 0.3408;
% std = 0.1179;
m = 1;

for std = 0.05:0.002:0.3;
thre = [mean-1.5*std,mean,mean+1.5*std,mean+3*std,mean+4.5*std]; % threshold for bar:gas

% figure;
sumpixel = [];

for k = 1:1:nfolder
    folder = direct(k+2).name;
    cd(folder);
    % read in and preprocess data
    [ gas,dissolved,rbc,barrier,mask,meanRbc2barrier ] = maps_preprocess('parameters');
    
    % load path data for data saving ;
    load('parameters');
    [pathstr,name,ext] = fileparts(bhute_pfile);
    
    % mask
    barrier = barrier.*mask;
    gas = gas .*mask;
    
    % initial ratio maps
    bar2gas = zeros(size(mask));
    
    % small value rule out
    % 0.01 of the average value of the biggest 5% of the data is taken as a
    % threshold
%     gmask = mask;
    [ gmask ] = f_smallruleout( gas,0.1,0.05);

    % average ratio 
    bar2gas(gmask) = barrier(gmask)./gas(gmask);

    
    % binning for barrier
    binning = bar2gas;
    binning = binning*TEc(k)*flipoff(k); % correction for TE and gas flip off
    
    binning(binning > thre(5)) = 7;
    binning(binning > thre(4) & binning <= thre(5)) = 6;
    binning(binning > thre(3) & binning <= thre(4)) = 5;
    binning(binning > thre(2) & binning <= thre(3)) = 4;
    binning(binning > thre(1) & binning <= thre(2)) = 3;
    binning(binning > 0 & binning <= thre(1)) = 2;
    binning(binning == 0) = 1;
    %save color binning
    %% Quantitative analysis
   a = binning == 2;
   clustermax(k,1) = sum(a(:));
   a = binning == 3;
   clustermax(k,2) = sum(a(:));
   a = binning == 4;
   clustermax(k,3) = sum(a(:));
   a = binning == 5;
   clustermax(k,4) = sum(a(:));
   a = binning == 6;
   clustermax(k,5) = sum(a(:));
   a = binning == 7;
   clustermax(k,6) = sum(a(:));

    cd('../');
end
cluster(m).name = ['std ',num2str(m)];
psum = sum(clustermax(:));
cluster(m).block1 = sum(clustermax(:,1))/psum;
cluster(m).block2 = sum(clustermax(:,2))/psum;
cluster(m).block3 = sum(clustermax(:,3))/psum;
cluster(m).block4 = sum(clustermax(:,4))/psum;
cluster(m).block5 = sum(clustermax(:,5))/psum;
cluster(m).block6 = sum(clustermax(:,6))/psum;
m = m+1;
end
%% plot;
std = 0.1:0.001:0.2;
a = zeros(6,101);
% for k =1:1:101
%     a(1,k) = cluster(k).block1;
%     a(2,k) = cluster(k).block2;
%     a(3,k) = cluster(k).block3;
%     a(4,k) = cluster(k).block4;
%     a(5,k) = cluster(k).block5;
%     a(6,k) = cluster(k).block6;
% end
for k =1:1:101
    a(1,k) = bar2gas_c(k).block1-bar2gas_ht(k).block1;
    a(2,k) = bar2gas_c(k).block2-bar2gas_ht(k).block2;
    a(3,k) = bar2gas_c(k).block3-bar2gas_ht(k).block3;
    a(4,k) = bar2gas_c(k).block4-bar2gas_ht(k).block4;
    a(5,k) = bar2gas_c(k).block5-bar2gas_ht(k).block5;
    a(6,k) = bar2gas_c(k).block6-bar2gas_ht(k).block6;
end
figure;
std = std/1.5;
plot(std,a(1,:),'r');hold on;
plot(std,a(2,:),'b'); hold on;
plot(std,a(3,:),'g'); hold on;
plot(std,a(4,:),'c');plot(std,a(5,:),'m');plot(std,a(6,:),'k');
legend('cluster1','cluster2','cluster3','cluster4','cluster5','cluster6');
grid on;
title('cluster IPF - cluster Healthy');
xlim([0.1/1.5,0.185/1.5]);