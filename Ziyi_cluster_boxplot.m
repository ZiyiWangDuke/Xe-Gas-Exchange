%% double tail t test
clear all;
load('cluster.mat');
block_ht_bar = [extractfield(cluster_ht_bar,'block1')',extractfield(cluster_ht_bar,'block2')',extractfield(cluster_ht_bar,'block3')',extractfield(cluster_ht_bar,'block4')',extractfield(cluster_ht_bar,'block5')',extractfield(cluster_ht_bar,'block6')',extractfield(cluster_ht_bar,'block7')',extractfield(cluster_ht_bar,'block8')'];
% block_ht_dis = [extractfield(cluster_ht_dis,'block1')',extractfield(cluster_ht_dis,'block2')',extractfield(cluster_ht_dis,'block3')',extractfield(cluster_ht_dis,'block4')',extractfield(cluster_ht_dis,'block5')',extractfield(cluster_ht_dis,'block6')'];
block_ht_rbc = [extractfield(cluster_ht_rbc,'block1')',extractfield(cluster_ht_rbc,'block2')',extractfield(cluster_ht_rbc,'block3')',extractfield(cluster_ht_rbc,'block4')',extractfield(cluster_ht_rbc,'block5')',extractfield(cluster_ht_rbc,'block6')'];
block_ipf_bar = [extractfield(cluster_ipf_bar,'block1')',extractfield(cluster_ipf_bar,'block2')',extractfield(cluster_ipf_bar,'block3')',extractfield(cluster_ipf_bar,'block4')',extractfield(cluster_ipf_bar,'block5')',extractfield(cluster_ipf_bar,'block6')',extractfield(cluster_ipf_bar,'block7')',extractfield(cluster_ipf_bar,'block8')'];
% block_ipf_dis = [extractfield(cluster_ipf_dis,'block1')',extractfield(cluster_ipf_dis,'block2')',extractfield(cluster_ipf_dis,'block3')',extractfield(cluster_ipf_dis,'block4')',extractfield(cluster_ipf_dis,'block5')',extractfield(cluster_ipf_dis,'block6')'];
block_ipf_rbc = [extractfield(cluster_ipf_rbc,'block1')',extractfield(cluster_ipf_rbc,'block2')',extractfield(cluster_ipf_rbc,'block3')',extractfield(cluster_ipf_rbc,'block4')',extractfield(cluster_ipf_rbc,'block5')',extractfield(cluster_ipf_rbc,'block6')'];

% cluster barrier
hp_bar = zeros(2,8);
for k = 1:1:8
    [hp_bar(1,k),hp_bar(2,k)] = ttest2(block_ht_bar(:,k),block_ipf_bar(1:20,k),'Tail','both','Vartype','unequal');
end
% cluster dissolved
% hp_dis = zeros(2,6);
% for k = 1:1:6
%     [hp_dis(1,k),hp_dis(2,k)] = ttest2(block_ht_dis(:,k),block_ipf_dis(:,k),'Tail','both','Vartype','unequal');
% end
% cluster RBC
hp_rbc = zeros(2,6);
for k = 1:1:6
    [hp_rbc(1,k),hp_rbc(2,k)] = ttest2(block_ht_rbc(:,k),block_ipf_rbc(1:20,k),'Tail','both','Vartype','unequal');
end
% ratio 
% hp_ratios = zeros(1,4);
% for k = 1:1:4
%     [hp_ratios(1,k),hp_ratios(2,k)] = ttest2(ratio_ipf(:,k),ratio_ht(:,k),'Tail','both','Vartype','unequal');
% end

%% box plot
num1 = 21; % number of ipf subjects
num2 = 10; % number of healthy subjects
cluster_ipf = cluster_ipf_rbc;
cluster_healthy = cluster_ht_rbc;
cur_p = hp_rbc;
n_clust = 6;

% cluster_ipf = cluster_ipf_rbc;
% cluster_healthy = cluster_ht_rbc;
% cur_p = hp_rbc;
% n_clust = 6;

names_ipf = fieldnames(cluster_ipf);
names_ipf = names_ipf(2:end,:);
names_ht = fieldnames(cluster_healthy);
names_ht = names_ht(2:end,:);
figure;
for k = 1:1:n_clust
    subplot(2,n_clust/2,k);
    ipf = extractfield(cluster_ipf,char(names_ipf(k)));
    ipf = ipf(1:num1-1);
    ht = extractfield(cluster_healthy,char(names_ht(k)));
    ht = ht(1:num2-1);
    sam = [ht ipf];
    grp = [zeros(1,num2-1),ones(1,num1-1)];
    boxplot(sam,grp,'labels',{'Healthy','IPF'});
    title(['Cluster',num2str(k),'  P: ',num2str(cur_p(2,k),'%1.4f')]);
%     title(['cluster',num2str(k),'  P Value:',num2str(hp_rbc(2,k))]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    ax = gca;
    ax.LineWidth = 2;
    ax.FontSize = 25;
%     boxplot(sam,grp,'notch','on','labels',{'Healthy Subject','IPF'},'whisker',1)
end
