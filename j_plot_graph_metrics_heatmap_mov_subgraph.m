%% Setup
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init


%% Load graphs

load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','grpahStatsTable_10_17'),'graphStatsTable')

results = struct();
corrType = 3;  % 1 - Similarity, 2 - All corr, 3 - Sig corr

for row =  1:size(graphStatsTable,1)

    switch corrType
        case 1
            c_rest = graphStatsTable{row,4}{:};
            c_run = graphStatsTable{row,5}{:};
            tt_clust = {'Clustering coefficient','Cosine similarity'};
            tt_close = {'Closeness centrality coefficient','cosine similarity'};
        case 2
            c_rest = graphStatsTable{row,8}{:};
            c_run = graphStatsTable{row,9}{:};
            tt_clust = {'Clustering coefficient','Pearson correlation'};
            tt_close = {'Closeness centrality coefficient','Pearson correlation'};
        case 3
            c_rest = graphStatsTable{row,12}{:};
            c_run = graphStatsTable{row,13}{:};
            tt_clust = {'Clustering coefficient','Pearson correlation sig'};
            tt_close = {'Closeness centrality coefficient','Pearson correlation sig'};
            mov_rest =  graphStatsTable{row,16}{:};
            mov_run =  graphStatsTable{row,17}{:};
    end

    d_rest = sqrt(-log(c_rest+eye(size(c_rest))));
    d_rest = d_rest-diag(diag(d_rest));
    g_rest = graph(d_rest,'omitselfloops');
    g_rest = rmedge(g_rest,find(isinf(g_rest.Edges{:,2})));
    cl_coeff_rest = getClusteringCoeff(g_rest);
    wcc_rest = size(c_rest,1)*centrality(g_rest,'closeness','Cost',g_rest.Edges.Weight);
    clust_rest = mean((cl_coeff_rest)*100);
    close_rest_res = mean((wcc_rest(mov_rest))*100);
    close_rest_nonres = mean((wcc_rest(~mov_rest))*100);


    d_run = sqrt(-log(c_run+eye(size(c_run))));
    d_run = d_run-diag(diag(d_run));
    g_run = graph(d_run,'omitselfloops');
    % Running network 
    g_run = rmedge(g_run,find(isinf(g_run.Edges{:,2})));
    wcc_run = size(c_run,1)*centrality(g_run,'closeness','Cost',g_run.Edges.Weight);
    cl_coeff_run = getClusteringCoeff(g_run);
    clust_run = mean((cl_coeff_run)*100);
    close_run = mean((wcc_run)*100);

    % Running responsive network 
    g_run_res = subgraph(g_run,mov_run);
    wcc_run_res = size(g_run_res.Nodes,1)*centrality(g_run_res,'closeness','Cost',g_run_res.Edges.Weight);
    close_run_res = mean((wcc_run_res)*100);

    g_run_nonres = subgraph(g_run,~mov_run);
    wcc_run_nonres = size(g_run_nonres.Nodes,1)*centrality(g_run_nonres,'closeness','Cost',g_run_nonres.Edges.Weight);
    close_run_nonres = mean((wcc_run_nonres)*100);
   

    % Assign session details
    results(row).animal = graphStatsTable{row,1}{:};
    results(row).day = graphStatsTable{row,2};
    results(row).condition = graphStatsTable{row,3}{:};

    results(row).g_rest = g_rest;
    results(row).g_run = g_run;
    results(row).g_run_res = g_run_res;
    results(row).g_run_nonres = g_run_nonres;
    results(row).clust_rest = clust_rest;
    results(row).clust_run = clust_run;
    results(row).close_rest_res = close_rest_res;
    results(row).close_run_res = close_run_res;
    results(row).close_rest_nonres = close_rest_nonres;
    results(row).close_run_nonres = close_run_nonres;
    results(row).heat_rest = wcc_rest;
    results(row).heat_run = wcc_run;
end


%% Plot boxplots of closeness of mov responsive vs non 
close all

% responsive    
close_WT_rest_res = [results(21:end).close_rest_res];
close_WT_rest_nonres = [results(21:end).close_rest_nonres];
close_WT_rest_delta = close_WT_rest_res-close_WT_rest_nonres;

close_KO_rest_res = [results(1:20).close_rest_res];
close_KO_rest_nonres = [results(1:20).close_rest_nonres];
close_KO_rest_delta = close_KO_rest_res-close_KO_rest_nonres;

close_WT_run_res = [results(21:end).close_run_res];
close_WT_run_nonres = [results(21:end).close_run_nonres];
close_WT_run_delta = close_WT_run_res-close_WT_run_nonres;

close_KO_run_res = [results(1:20).close_run_res];
close_KO_run_nonres = [results(1:20).close_run_nonres];
close_KO_run_delta = close_KO_run_res-close_KO_run_nonres;

%% Boxplot - Closeness centrality 

%% T test for Run mov modulated - non modulated 

figure('units','normalized','position',[0.2,0.4,0.5,0.4])
b1 = boxchart(1*ones(1,numel(close_WT_run_delta)),...
    close_WT_run_delta,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'WhiskerLineColor', colors(2,:),'LineWidth',1);
hold on
s1 = scatter(1,close_WT_run_delta,48,colors(2,:),'filled');

b3 = boxchart(2*ones(1,numel(close_KO_run_delta)),...
    close_KO_run_delta,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'WhiskerLineColor', colors(4,:),'LineWidth',1);
hold on
s3 = scatter(2,close_KO_run_delta,48,colors(4,:),'filled');
xticks([1,2])
xticklabels({'WT Run','KO Run'})
yy = ylim;
[hrun, prun] = ttest2(close_WT_run_delta,close_KO_run_delta);
text(1.5, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prun),'HorizontalAlignment','center')
title({'Mean closeness centrality of','Mov modulated cells - Mov non modulated cells'})
axis square
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
box off
ylim([-inf,25])

%% T test for Run mov modulated vs non modulated 
figure('WindowStyle','docked')
b1 = boxchart(1*ones(1,numel(close_WT_run_res)),...
    close_WT_run_res,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'WhiskerLineColor', colors(2,:),'LineWidth',1);
hold on
s1 = scatter(1,close_WT_run_res,48,colors(2,:),'filled');

b2 = boxchart(2*ones(1,numel(close_WT_run_nonres)),...
    close_WT_run_nonres,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',0.3*[1,1,1],...
    'WhiskerLineColor',0.3*[1,1,1],'LineWidth',1);
hold on
s2 = scatter(2,close_WT_run_nonres,48,0.3*[1,1,1],'filled');

b3 = boxchart(3*ones(1,numel(close_KO_run_res)),...
    close_KO_run_res,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'WhiskerLineColor', colors(4,:),'LineWidth',1);
hold on
s3 = scatter(3,close_KO_run_res,48,colors(4,:),'filled');

b4 = boxchart(4*ones(1,numel(close_KO_run_nonres)),...
    close_KO_run_nonres,...
    'BoxFaceAlpha',0,...
    'MarkerStyle','none',...
    'BoxFaceColor',0.3*[1,1,1],...
    'WhiskerLineColor',0.3*[1,1,1],'LineWidth',1);
hold on
s4 = scatter(4,close_KO_run_nonres,48,0.3*[1,1,1],'filled');

xticks([1,2,3,4])
xticklabels({'WT-Run-MovMod','WT-Run-MovNonMod','KO-Run-MovMod','KO-Run-MovNonMod'})
[hWT, pWT] = ttest2(close_WT_run_res,close_WT_run_nonres);
yy = ylim;
text(1.5, yy(2)*1.05, sprintf('P-ttest: %.3f' ,pWT),'HorizontalAlignment','center')

[hKO, pKO] = ttest2(close_KO_run_res,close_KO_run_nonres);
yy = ylim;
text(3.5, yy(2)*1.05, sprintf('P-ttest: %.3f' ,pKO),'HorizontalAlignment','center')
title({'Mean closeness centrality of','movement modulated vs movement non modulated','Subnetworks'})
axis square
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
box off
ylim([-inf,45])

% means = [mean(close_WT_rest_delta),mean(close_WT_run),mean(close_KO_rest_delta),mean(close_KO_run_res)]';
% stds = [std(close_WT_rest_delta),std(close_WT_run),std(close_KO_rest_delta),std(close_KO_run_res)]';

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('closeness_centrality_%s.png',strjoin(strsplit(tt_clust{2}),'_'))));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('closeness_centrality_%s.epsc',strjoin(strsplit(tt_clust{2}),'_'))));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\boxplots',...
%     sprintf('closeness_centrality_%s.fig',strjoin(strsplit(tt_clust{2}),'_'))));
% 

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
    g_run_res = results(row).g_run_res;
    g_run_nonres =  results(row).g_run_nonres;
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
            mov_rest =  double(graphStatsTable{row,16}{:});
            mov_run =  double(graphStatsTable{row,17}{:});
    end
    
    % Get edge weight 
    e_rest = graph(c_rest,'omitselfloops');
    LWidths_rest = 0.75*e_rest.Edges.Weight/max(e_rest.Edges.Weight);

    e_run = graph(c_run,'omitselfloops');
    LWidths_run = 0.75*e_run.Edges.Weight/max(e_run.Edges.Weight);


    % plot network with physical roi location
%     figure('Units','Normalized','Position',[0 0 1 1])
%     subplot(1,2,1)
%     g1 = plot(g_rest,'NodeColor',[.3185,.3185,.3185],'EdgeColor',[0.4196 0.4196 0.4196],'XData',cen_rest(:,1),'YData',cen_rest(:,2),'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_rest,'EdgeAlpha',1); %,'NodeLabel',cen_rest(:,3),'EdgeLabel',round(g_rest.Edges.Weight,1));
%     g1.NodeCData = results(row).heat_rest;
%     colormap jet
%     axis square
%     colorbar
%     caxis([0.2 0.35])
%     title([tt, sprintf('%s day %i Rest',results(row).animal,results(row).day)])
% 
%     subplot(1,2,2)
%     g2 = plot(g_run,'NodeColor',[.3185,.3185,.3185],'EdgeColor',[0.4196 0.4196 0.4196],'XData',cen_run(:,1),'YData',cen_run(:,2),'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_run,'EdgeAlpha',1);%'EdgeAlpha',0.5,,'NodeLabel',cen_run(:,3),'EdgeLabel',round(g_rest.Edges.Weight,1));
%     g2.NodeCData = results(row).heat_run;
%     colormap jet
%     axis square
%     colorbar
%     caxis([0.2 0.35])
%     title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day)])
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.png',results(row).animal,results(row).day,tt{:})));
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.svg',results(row).animal,results(row).day,tt{:})));
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\weighted edges and node color",sprintf('%s day %i physical %s.fig',results(row).animal,results(row).day,tt{:})));


    %% plot schematic network
%   figure('Units','Normalized','Position',[0 0 1 1])
    figure('WindowStyle','docked')
    
    subplot(1,3,1)
    NodeCData = [.3185,.3185,.3185]; %repmat([0,0,0],size(g_rest.Nodes,1),1);
    g3 = plot(g_run,'Layout','force','WeightEffect','direct','NodeColor', NodeCData,'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',5,'NodeLabel', {},'LineWidth',LWidths_run,'EdgeAlpha',1); %,'EdgeAlpha',0.5,'NodeLabel',cen_rest(:,3),'EdgeLabel',round(g_run.Edges.Weight,1));
    g3.NodeCData = mov_run;
    colormap jet
    axis square
    colorbar
    caxis([0 1])
    title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day)])

    subplot(1,3,2)
    NodeCData = [.3185,.3185,.3185]; % repmat([0,0,0],size(g_run.Nodes,1),1);
    g4 = plot(g_run_res,'Layout','force','WeightEffect','direct','NodeColor',NodeCData,'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',5,'NodeLabel', {},'EdgeAlpha',1); %,'EdgeAlpha',0.5,'NodeLabel',cen_run(:,3),'EdgeLabel',round(g_run.Edges.Weight,1));
    g4.NodeCData = ones(size(g_run_res.Nodes,1),1);
    colormap jet
    axis square
    colorbar
    caxis([0 1])
    title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day),'Modulated Subgraph'])

    subplot(1,3,3)
    NodeCData = [.3185,.3185,.3185]; % repmat([0,0,0],size(g_run.Nodes,1),1);
    g5 = plot(g_run_nonres,'Layout','force','WeightEffect','direct','NodeColor',NodeCData,'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',5,'NodeLabel', {},'EdgeAlpha',1); %,'EdgeAlpha',0.5,'NodeLabel',cen_run(:,3),'EdgeLabel',round(g_run.Edges.Weight,1));
    g5.NodeCData = zeros(size(g_run_nonres.Nodes,1),1);
    colormap jet
    axis square
    colorbar
    caxis([0 1])
    title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day),'Non Modulated Subgraph'])

 
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\mov responsive cells subgraph",sprintf('%s day %i schematic %s.png',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\mov responsive cells subgraph",sprintf('%s day %i schematic %s.svg',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\j_graphs\network_sigcorr\mov responsive cells subgraph",sprintf('%s day %i schematic %s.fig',results(row).animal,results(row).day,tt{:})));

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

 