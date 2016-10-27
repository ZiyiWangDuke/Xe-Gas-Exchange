function f_ratiomap_colorhistogram(line,sub,key,pathhist)
% function plot histogram with the reference
% line, the aggregated vonxel intensities
% key, 1 for barrier and 0 for RBC
% sub, subject name
if ~exist('pathhist','var')
   pathhist = 'C:\Users\zw73\Desktop\ratipcolors\histograms\';
end
switch key 
    %% plot the histogram for RBC:gas
    case 0
    dis_thre = 0.8;
    nbins = 50;
    inter = dis_thre/nbins;
    map = brewermap(4,'Set1');
    figure;
    line = [line;dis_thre;0]; % standerized the dynamic range
    line = line(line<=dis_thre);
    filename = [pathhist,'hist_rbc_',sub];
    histogram(line,nbins,'facecolor',map(1,:),'facealpha',.6,'Normalization','probability');
    xlim([0,dis_thre]);
    ylim([0,0.1]);
%     ylim([0,0.15]);
    % hold on;
    
   
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',26)
    ylabel('Fraction of total pixels', 'FontSize',26)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 25;
    
    hold on;
    
    x = 0:inter:dis_thre;
    a1 = 0.06106;
    b1 = 0.2604;
    c1 = 0.1481 ;

    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',4);
    print(filename,'-dpng');close;
    %% plot the histogram for barrier:gas
    case 1 
    dis_thre = 1.8;
    nbins = 70;
    inter = dis_thre/nbins;
    map = brewermap(4,'Set2');
    figure;
    line = line(line<=dis_thre);
    line = [line;dis_thre;0];
    % plot
    filename = [pathhist,'hist_bar_',sub];
    histogram(line,nbins,'facecolor',map(1,:),'facealpha',.6,'Normalization','probability');
    xlim([0,dis_thre]);
    ylim([0,0.18]);
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',26)
    ylabel('Fraction of total pixels', 'FontSize',26)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 25;
    
    hold on;
    
    x = 0:inter:dis_thre;
    a1 = 0.07006;
    b1 = 0.4632;
    c1 = 0.1953 ;
    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',4);
    print(filename,'-dpng');close;
    case 2
    nbins = 50;
    inter = 1/nbins;
    map = brewermap(4,'Blues');
    figure;
    filename = [pathhist,'hist_vent_',sub];
    histogram(line,nbins,'facecolor',map(3,:),'facealpha',.6,'Normalization','probability');
    xlim([0,1]);
    ylim([0,0.07]);
    % title('Histogram after Scaling', 'FontSize',12);
    xlabel('Pixel Intensity', 'FontSize',26)
    ylabel('Fraction of total pixels', 'FontSize',26)
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 25;
    
    hold on;
    
    x = 0:inter:1;
    a1 = 0.04074;
    b1 = 0.494;
    c1 = 0.278 ;
    normal = a1*exp(-((x-b1)/c1).^2);
    area(x,normal,'FaceColor','none','LineStyle','--','LineWidth',4);
    print(filename,'-dpng');close; 
end
end