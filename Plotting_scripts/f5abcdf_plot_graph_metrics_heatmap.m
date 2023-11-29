%% Setup
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init

%% Load graphs

load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','grpahStatsTable_08_26'),'graphStatsTable')

results = struct();
corrType = 3;  % 1 - Similarity, 2 - All corr, 3 - Sig corr

for row =  1:size(graphStatsTable,1)

    switch corrType
        case 1
            c_rest = graphStatsTable{row,4}{:};
            c_run = graphStatsTable{row,5}{:};
            tt_clust = {'Clustering coefficient','Cosine similarity'};
            tt_close = {'Closeness centrality','cosine similarity'};
        case 2
            c_rest = graphStatsTable{row,8}{:};
            c_run = graphStatsTable{row,9}{:};
            tt_clust = {'Clustering coefficient','Pearson correlation'};
            tt_close = {'Closeness centrality','Pearson correlation'};
        case 3
            c_rest = graphStatsTable{row,12}{:};
            c_run = graphStatsTable{row,13}{:};
            tt_clust = {'Clustering coefficient','Pearson correlation sig'};
            tt_close = {'Closeness centrality','Pearson correlation sig'};
    end

    d_rest = sqrt(-log(c_rest+eye(size(c_rest))));
    d_rest = d_rest-diag(diag(d_rest));
    g_rest = graph(d_rest,'omitselfloops');
%     inactive_cells = graphStatsTable{row,16}{:};
%     g_rest = rmnode(g_rest,inactive_cells);

    g_rest = rmedge(g_rest,find(isinf(g_rest.Edges{:,2})));
    cl_coeff_rest = getClusteringCoeff(g_rest);
    wcc_rest = size(c_rest,1)*centrality(g_rest,'closeness','Cost',g_rest.Edges.Weight);

    clust_rest = mean((cl_coeff_rest));
    close_rest = mean((wcc_rest));


    d_run = sqrt(-log(c_run+eye(size(c_run))));
    d_run = d_run-diag(diag(d_run));
    g_run = graph(d_run,'omitselfloops');
    % Add centroid here
%     inactive_cells = graphStatsTable{row,17}{:};
%     g_run = rmnode(g_run,inactive_cells);
    g_run = rmedge(g_run,find(isinf(g_run.Edges{:,2})));
    wcc_run = size(c_run,1)*centrality(g_run,'closeness','Cost',g_run.Edges.Weight);
    cl_coeff_run = getClusteringCoeff(g_run);
    %     [~ ,i_run]= min((cl_coeff_run)*100);

    clust_run = mean((cl_coeff_run));
    close_run = mean((wcc_run));



    % Assign session details
    results(row).animal = graphStatsTable{row,1}{:};
    results(row).day = graphStatsTable{row,2};
    results(row).condition = graphStatsTable{row,3}{:};

    results(row).g_rest = g_rest;
    results(row).g_run = g_run;
    results(row).clust_rest = clust_rest;
    results(row).clust_run = clust_run;
    results(row).close_rest = close_rest;
    results(row).close_run = close_run;
    results(row).heat_rest = wcc_rest;
    results(row).heat_run = wcc_run;
end


%% Plot boxplots
close all

% 
% clust_WT_rest = [results(21:end).clust_rest];
% clust_KO_rest = [results(1:20).clust_rest];
% clust_WT_run = [results(21:end).clust_run];
% clust_KO_run = [results(1:20).clust_run];

close_WT_rest = [results(21:end).close_rest];
close_KO_rest = [results(1:20).close_rest];
close_WT_run = [results(21:end).close_run];
close_KO_run = [results(1:20).close_run];

% % Boxplot - Clustering coefficient
% 
% figure('units','normalized','position',[0.2,0.4,0.5,0.4])
% b1 = boxchart(1*ones(1,numel(clust_WT_rest)),...
%     clust_WT_rest,...
%     'BoxFaceAlpha',0,...
%     'MarkerStyle','none',...
%     'BoxFaceColor',colors(1,:),...
%     'WhiskerLineColor', colors(1,:),'LineWidth',1);
% hold on
% s1 = scatter(1,clust_WT_rest,48,colors(1,:),'filled');
% 
% b2 = boxchart(2*ones(1,numel(clust_WT_run)),...
%     clust_WT_run,...
%     'BoxFaceAlpha',0,...
%     'MarkerStyle','none',...
%     'BoxFaceColor',colors(2,:),...
%     'WhiskerLineColor', colors(2,:),'LineWidth',1);
% hold on
% s2 = scatter(2,clust_WT_run,48,colors(2,:),'filled');
% 
% b3 = boxchart(3*ones(1,numel(clust_KO_rest)),...
%     clust_KO_rest,...
%     'BoxFaceAlpha',0,...
%     'MarkerStyle','none',...
%     'BoxFaceColor',colors(3,:),...
%     'WhiskerLineColor', colors(3,:),'LineWidth',1);
% hold on
% s3 = scatter(3,clust_KO_rest,48,colors(3,:),'filled');
% 
% b4 = boxchart(4*ones(1,numel(clust_KO_run)),...
%     clust_KO_run,...
%     'BoxFaceAlpha',0,...
%     'MarkerStyle','none',...
%     'BoxFaceColor',colors(4,:),...
%     'WhiskerLineColor', colors(4,:),'LineWidth',1);
% hold on
% s4 = scatter(4,clust_KO_run,48,colors(4,:),'filled');
% 
% xticks([1,2,3,4])
% xticklabels({'WT Rest','WT Run','KO Rest','KO Run'})
% [hrest, prest] = ttest2(clust_WT_rest,clust_KO_rest);
% [hrun, prun] = ttest2(clust_WT_run,clust_KO_run);
% [hWT, pWT] = ttest(clust_WT_rest,clust_WT_run);
% [hKO, pKO] = ttest(clust_KO_rest,clust_KO_run);
% ttests = [pWT,pKO,prest,prun]';
% means = [mean(clust_WT_rest),mean(clust_WT_run),mean(clust_KO_rest),mean(clust_KO_run)]';
% stds = [std(clust_WT_rest),std(clust_WT_run),std(clust_KO_rest),std(clust_KO_run)]';
% 
% yy = ylim;
% text(1.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pWT),'HorizontalAlignment','center')
% % text(2, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prest),'HorizontalAlignment','center')
% % text(3, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prun),'HorizontalAlignment','center')
% text(3.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pKO),'HorizontalAlignment','center')
% title(tt_clust')
% axis square
% set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
% box off
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('clustering_coefficient_%s.png',strjoin(strsplit(tt_clust{2}),'_'))));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('clustering_coefficient_%s.epsc',strjoin(strsplit(tt_clust{2}),'_'))));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('clustering_coefficient_%s.fig',strjoin(strsplit(tt_clust{2}),'_'))));

%% 4 Groups - Normality test
sessionwise_stats =struct();
 
yData = [close_WT_rest';close_WT_run';close_KO_rest';close_KO_run'];
gType = [ones(size(close_WT_rest)),ones(size(close_WT_run)),zeros(size(close_KO_rest)),zeros(size(close_KO_run))]';
movSp = [zeros(size(close_WT_rest)),ones(size(close_WT_run)),zeros(size(close_KO_rest)),ones(size(close_KO_run))]';

figure;

names = {'WTrest';'WTrun';'KOrest';'KOrun'};
mat = [1,0;1,1;0,0;0,1];
isNormal = [];
for pp = 1:4
    yData_group = yData(gType == mat(pp,1) & movSp == mat(pp,2));
    [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
    subplot(2,2,pp)
    qqplot(yData_group);
    title({names{pp}, sprintf('SW test p: %.4f',p_SW)})
    isNormal = [isNormal;[~h_SW,p_SW]];
end
sessionwise_stats(1).type = '4_groups';
sessionwise_stats(1).isNormal= isNormal';

sessionwise_stats(1).means = [mean(close_WT_rest),mean(close_WT_run),mean(close_KO_rest),mean(close_KO_run)]';
sessionwise_stats(1).stds = [std(close_WT_rest),std(close_WT_run),std(close_KO_rest),std(close_KO_run)]';


% Do two way anova
[sessionwise_stats(1).ANOVA.p...
    sessionwise_stats(1).ANOVA.tbl...
    sessionwise_stats(1).ANOVA.stats...
    sessionwise_stats(1).ANOVA.terms] = anovan(yData,{gType,movSp},"model","interaction",'display','off');


% FIT GLM
sessionwise_stats(1).GLM.mdl = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl);
sessionwise_stats(1).GLM.fit = deviance_test{2,4};


% Boxplot - Closeness centrality 
figure('units','normalized','position',[0.2,0.4,0.5,0.4])
b1 = boxchart(1*ones(1,numel(close_WT_rest)),...
    close_WT_rest,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'WhiskerLineColor', colors(1,:),'LineWidth',1);
hold on
s1 = scatter(1,close_WT_rest,48,colors(1,:),'filled');

b2 = boxchart(2*ones(1,numel(close_WT_run)),...
    close_WT_run,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'WhiskerLineColor', colors(2,:),'LineWidth',1);
hold on
s2 = scatter(2,close_WT_run,48,colors(2,:),'filled');

b3 = boxchart(3*ones(1,numel(close_KO_rest)),...
    close_KO_rest,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'WhiskerLineColor', colors(3,:),'LineWidth',1);
hold on
s3 = scatter(3,close_KO_rest,48,colors(3,:),'filled');

b4 = boxchart(4*ones(1,numel(close_KO_run)),...
    close_KO_run,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'WhiskerLineColor', colors(4,:),'LineWidth',1);
hold on
s4 = scatter(4,close_KO_run,48,colors(4,:),'filled');

xticks([1,2,3,4])
xticklabels({'WT Rest','WT Run','KO Rest','KO Run'})

% t tests
% [hrest, prest] = ttest2(close_WT_rest,close_KO_rest);
% [hrun, prun] = ttest2(close_WT_run,close_KO_run);
% [hWT, pWT] = ttest(close_WT_rest,close_WT_run);
% [hKO, pKO] = ttest(close_KO_rest,close_KO_run);
% ttests = [pWT,pKO,prest,prun]';

% yy = ylim;
% text(1.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pWT),'HorizontalAlignment','center')
% text(2, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prest),'HorizontalAlignment','center')
% text(3, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prun),'HorizontalAlignment','center')
% text(3.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pKO),'HorizontalAlignment','center')
title(tt_close')

axis square
set(gca,'Units','inches','InnerPosition',[.8 .75 3 4], 'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
box off
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'closeness_centrality.png'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'closeness_centrality.epsc'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'closeness_centrality.fig'));

%% 2 Groups - T test

% Normality test
yData = [close_WT_rest'-close_WT_run';close_KO_rest'-close_KO_run'];
gType = [ones(size(close_WT_rest)),zeros(size(close_KO_rest))]';
% movSp = [zeros(size(close_WT_rest)),ones(size(close_WT_run)),zeros(size(close_KO_rest)),ones(size(close_KO_run))]';

figure;
names = {'WT';'KO'};
mat = [1;0];
isNormal = [];
for pp = 1:2
    yData_group = yData(gType == mat(pp));
    [h_SW, p_SW, W] = swtest(yData_group, 0.05,1);
    subplot(1,2,pp)
    qqplot(yData_group);
    title({names{pp}, sprintf('SW test p: %.4f',p_SW)})
    isNormal = [isNormal;[~h_SW,p_SW]];
end

sessionwise_stats(2).type = '2_groups';
sessionwise_stats(2).isNormal= isNormal';
sessionwise_stats(2).means = [mean(close_WT_rest-close_WT_run),mean(close_KO_rest-close_KO_run)]';
sessionwise_stats(2).stds = [std(close_WT_rest-close_WT_run),std(close_KO_rest-close_KO_run)]';

% t test
[sessionwise_stats(2).Ttest.h_t, sessionwise_stats(2).Ttest.p_t] = ttest2(close_WT_rest'-close_WT_run',close_KO_rest'-close_KO_run');
% ranksum
[sessionwise_stats(2).ranksum.h_rs, sessionwise_stats(2).ranksum.p_rs] = ranksum(close_WT_rest'-close_WT_run',close_KO_rest'-close_KO_run');
% % ranksum
% [sessionwise_stats(2).ranksum.h_WT, sessionwise_stats(2).ranksum.p_WT] = ranksum(close_WT_rest'-close_WT_run',close_KO_rest'-close_KO_run');
% signrank
[sessionwise_stats(2).signrank.h_WT, sessionwise_stats(2).signrank.p_WT] = signrank(close_WT_rest'-close_WT_run');
[sessionwise_stats(2).signrank.h_KO, sessionwise_stats(2).signrank.p_KO] = signrank(close_KO_rest'-close_KO_run');
sessionwise_stats(2).signrank.median_WT = median(close_WT_rest'-close_WT_run');
sessionwise_stats(2).signrank.median_KO = median(close_KO_rest'-close_KO_run');

% Boxplot - Closeness centrality difference
figure('units','normalized','position',[0.2,0.4,0.5,0.4])
b1 = boxchart(1*ones(1,numel(close_WT_rest - close_WT_run)),...
    close_WT_rest-close_WT_run,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'WhiskerLineColor', colors(2,:),'LineWidth',1);
hold on
s1 = scatter(1,close_WT_rest -close_WT_run,48,colors(2,:),'filled');

b2 = boxchart(2*ones(1,numel(close_KO_rest-close_KO_run)),...
    close_KO_rest-close_KO_run,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'WhiskerLineColor', colors(4,:),'LineWidth',1);
hold on
s2 = scatter(2,close_KO_rest-close_KO_run,48,colors(4,:),'filled');
axis square
set(gca,'Units','inches','InnerPosition',[.8 .75 3 2], 'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
xticks([1,2])
xticklabels({'WT','KO'})
box off
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'diff_closeness_centrality.png'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'diff_closeness_centrality.epsc'));
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs',...
    'diff_closeness_centrality.fig'));

save('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\j_plot_graph_metrics_heatmap.mat','sessionwise_stats')

%% Plot networks
% close all


for row = 1:32
    % decide color 
    if row>=21 % WT
        edge_col_run = colors(2,:);
        edge_col_rest = colors(1,:);
    else
        edge_col_run = colors(4,:);
        edge_col_rest = colors(3,:);
    end

    g_rest = results(row).g_rest;
    g_run =  results(row).g_run;
    % Get centroid and edge weights 
    switch corrType
        case 1
            cen_rest = graphStatsTable{row,6}{:};
            cen_run = graphStatsTable{row,7}{:};
            c_rest = graphStatsTable{row,4}{:};
            c_run = graphStatsTable{row,5}{:};
            tt = {'Cosine similarity'};
        case 2
            cen_rest = graphStatsTable{row,10}{:};
            cen_run = graphStatsTable{row,11}{:};
            c_rest = graphStatsTable{row,8}{:};
            c_run = graphStatsTable{row,9}{:};
            tt = {'Pearson correlation'};

        case 3
            cen_rest = graphStatsTable{row,14}{:};
            cen_run = graphStatsTable{row,15}{:};
            c_rest = graphStatsTable{row,12}{:};
            c_run = graphStatsTable{row,13}{:};
            tt = {'Pearson correlation sig'};
    end
    
    % Get edge weight 
    e_rest = graph(c_rest,'omitselfloops');
    LWidths_rest = 0.75*e_rest.Edges.Weight/max(e_rest.Edges.Weight);

    e_run = graph(c_run,'omitselfloops');
    LWidths_run = 0.75*e_run.Edges.Weight/max(e_run.Edges.Weight);


    % plot network with physical roi location
    figure('Units','Normalized','Position',[0 0 1 1])
    subplot(1,2,1)
    g1 = plot(g_rest,'NodeColor',[.3185,.3185,.3185],'EdgeColor',[0.4196 0.4196 0.4196],'XData',cen_rest(:,1),'YData',cen_rest(:,2),'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_rest,'EdgeAlpha',1); %,'NodeLabel',cen_rest(:,3),'EdgeLabel',round(g_rest.Edges.Weight,1));
    g1.NodeCData = results(row).heat_rest;
    colormap jet
    axis square
    colorbar
    caxis([0.2 0.35])
    title([tt, sprintf('%s day %i Rest',results(row).animal,results(row).day)])

    subplot(1,2,2)
    g2 = plot(g_run,'NodeColor',[.3185,.3185,.3185],'EdgeColor',[0.4196 0.4196 0.4196],'XData',cen_run(:,1),'YData',cen_run(:,2),'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_run,'EdgeAlpha',1);%'EdgeAlpha',0.5,,'NodeLabel',cen_run(:,3),'EdgeLabel',round(g_rest.Edges.Weight,1));
    g2.NodeCData = results(row).heat_run;
    colormap jet
    axis square
    colorbar
    caxis([0.2 0.35])
    title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day)])
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.png',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.svg',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.fig',results(row).animal,results(row).day,tt{:})));


    % plot schematic network
    figure('Units','Normalized','Position',[0 0 1 1])
    subplot(1,2,1)
    NodeCData = [.3185,.3185,.3185]; %repmat([0,0,0],size(g_rest.Nodes,1),1);
    g3 = plot(g_rest,'Layout','force','WeightEffect','direct','NodeColor', NodeCData,'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_rest,'EdgeAlpha',1); %,'EdgeAlpha',0.5,'NodeLabel',cen_rest(:,3),'EdgeLabel',round(g_run.Edges.Weight,1));
    g3.NodeCData = results(row).heat_rest;
    colormap jet
    axis square
    colorbar
    caxis([0.2 0.35])
    title([tt, sprintf('%s day %i Rest',results(row).animal,results(row).day)])

    subplot(1,2,2)
    NodeCData = [.3185,.3185,.3185]; % repmat([0,0,0],size(g_run.Nodes,1),1);
    g4 = plot(g_run,'Layout','force','WeightEffect','direct','NodeColor',NodeCData,'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_run,'EdgeAlpha',1); %,'EdgeAlpha',0.5,'NodeLabel',cen_run(:,3),'EdgeLabel',round(g_run.Edges.Weight,1));
    g4.NodeCData = results(row).heat_run;
    colormap jet
    axis square
    colorbar
    caxis([0.2 0.35])
    title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day)])
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i schematic %s.png',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i schematic %s.svg',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i schematic %s.fig',results(row).animal,results(row).day,tt{:})));

    close all

end


% % %% add colorbar and set colorbar for all 8 plots included in paper
% % figs = struct;
% % 
% % % get files
% % [file,path] = uigetfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color\*.fig','Select figures','MultiSelect','on');
% % 
% % % open each figure and get its graph plot properties
% % for i = 1:size(file,2)
% %     filename=fullfile(path,file(i));
% %     open(string(filename));  
% %     g = get(gcf,'Children');
% %     figs(k) = g(1).Children;
% %     figs(k+1) = g(2).Children;
% %     k=k+2;
% % end

 