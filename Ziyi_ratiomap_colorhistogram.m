%% This file plot the colorhistogram for each single 
% sumpixel = sort(sumpixel,'descend');
%% fit healthy reference
key = 2; % 2 for normalized ventilation, 1 for barrier and 0 for RBC
% range for rbc:gas is 0 - 0.8;
% range for bar:gas is 0 - 1.8;
% range for normalzied gas is 0 - 1;
% sumpixel = sumpixel(sumpixel<=1.8);
sumpixel = [sumpixel,1.8];
sumpixel(sumpixel>1.8) = 1.8;
numsum = length(sumpixel);
[disy,disx] = hist(sumpixel,70,'Normalization','probability');
disy = disy/numsum;
%% plot the histogram for RBC:gas
switch key
    case 0
    dis_thre = 0.5;
    nbins = 50;
    inter = dis_thre/nbins;
    map = brewermap(4,'Set2');
    figure;
%     line = [line;dis_thre]; % standerized the dynamic range
%     line(line>dis_thre) = dis_thre;
    histogram(line,nbins,'facecolor',map(1,:));%,'facealpha',.6);
%     xlim([0,0.5]);
    ylim([0,1500]);
    % hold on;
    
   
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',26)
    ylabel('Fraction of total pixels', 'FontSize',26)
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 25;
    
    hold on;
    
    x = 0:inter:dis_thre;
    a1 = 0.08056;
    b1 = 0.2367;
    c1 = 0.1369;
    % a2 = 162.8;
    % b2 =0.3389;
    % c2 =0.1066;
    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',2);
    case 1
    %% plot the histogram for barrier:gas
    dis_thre = 1.4;
    nbins = 70;
    inter = dis_thre/nbins;
    map = brewermap(4,'Set1');
    figure;
    line = [line;dis_thre];
    line(line>dis_thre) = dis_thre;
    % plot
    histogram(line,nbins,'facecolor',map(1,:));%,'facealpha',.6)
    xlim([0,dis_thre]);
    ylim([0,4000]);
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',30)
    ylabel('Fraction of total pixels', 'FontSize',30)
    
    hold on;
    
    x = 0:inter:dis_thre;
    a1 = 0.1186;
    b1 = 0.4106;
    c1 = 0.1638 ;
    % a2 = 0.0849;
    % b2 =0.4105;
    % c2 =0.1633;
    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',2);
    case 2
    %% plot the histogram for normalzied gas
    nbins = 150;
    inter = 1/nbins;
%     map = brewermap('list');
    map = brewermap(4,'Blues');
    figure;
    % plot
    histogram(sumpixel,nbins,'facecolor',map(3,:),'facealpha',.6,'Normalization','probability');
    xlim([0,1]);
%     ylim([0,4000]);
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',30)
    ylabel('Fraction of total pixels', 'FontSize',30)
    
    hold on;
    
    x = 0:inter:1;
    a1 = 0.01359;
    b1 = 0.4939;
    c1 = 0.278 ;
    % a2 = 0.0849;
    % b2 =0.4105;
    % c2 =0.1633;
    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',2);
end