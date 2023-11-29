%% Setup
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init

%% Load graphs

load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','grpahStatsTable_02_23'),'graphStatsTable')


corr = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTablePearson_09_12'),'corrStatsTable');
corrStatsTable = corr.corrStatsTable;
% fieldNames = corrStatsTableSpeed.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

results = struct();
corrType = 3;  % 1 - Similarity, 2 - All corr, 3 - Sig corr

[divergent_rainbow_map] = diverging_map(linspace(0,1,256),rainbowmap(77,:),rainbowmap(end-10,:));

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
    inactive_cells_rest = find(sum(c_rest)==1); %graphStatsTable{row,16}{:};
    g_rest = rmnode(g_rest,inactive_cells_rest);
    g_rest = rmedge(g_rest,find(isinf(g_rest.Edges{:,2})));

    cl_coeff_rest = getClusteringCoeff(g_rest);
    wcc_rest = size(g_rest.Nodes,1)*centrality(g_rest,'closeness','Cost',g_rest.Edges.Weight);

    clust_rest = mean((cl_coeff_rest)*100);
    close_rest = mean((wcc_rest)*100);


    d_run = sqrt(-log(c_run+eye(size(c_run))));
    d_run = d_run-diag(diag(d_run));
    g_run = graph(d_run,'omitselfloops');
    % Add centroid here
    inactive_cells_run = find(sum(c_run)==1); %graphStatsTable{row,17}{:};
    g_run = rmnode(g_run,inactive_cells_run);
    g_run = rmedge(g_run,find(isinf(g_run.Edges{:,2})));
    wcc_run = size(g_run.Nodes,1)*centrality(g_run,'closeness','Cost',g_run.Edges.Weight);
    cl_coeff_run = getClusteringCoeff(g_run);
    %     [~ ,i_run]= min((cl_coeff_run)*100);

    clust_run = mean((cl_coeff_run)*100);
    close_run = mean((wcc_run)*100);

    %     heat_desynch =
    heat_run_all = zeros(max(size(c_run,1),size(c_rest,1)),1);
    heat_run_all(inactive_cells_run) = nan;
    heat_run_all(~isnan(heat_run_all)) = wcc_run;

    heat_rest_all = zeros(max(size(c_run,1),size(c_rest,1)),1);
    heat_rest_all(inactive_cells_rest) = nan;
    heat_rest_all(~isnan(heat_rest_all)) = wcc_rest;

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
    results(row).heat_desynch = heat_rest_all-heat_run_all;

end

%% Plot networks
% close all

for row = 1%1:32
    % decide color
    if row>=21 % WT
        edge_col = colors(2,:);
        %         edge_col = colors(1,:);
    else
        edge_col = colors(4,:);
        %         edge_col = colors(3,:);
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

            % network for whole trace is used to plot edges 
            corr_all  = corrStatsTable{row,4}{:};
            thresh_all  = corrStatsTable{row,5}{:};
            cell_close = double(corrStatsTable{row,7}{:}<20);
            n = size(corr_all);
            corr_all(logical(triu(ones(n)))) = nan;
            thresh_all(logical(triu(ones(n)))) = nan;

            corr_all_pos = corr_all.*(corr_all>0);
            sig_edges_mask_all = (corr_all-thresh_all >1e-5) & (cell_close == 0);
            corr_all_pos_sig = corr_all_pos.*sig_edges_mask_all;

            corr_all_pos(isnan(corr_all_pos)) = 0;
            corr_all_pos = 0.5*(corr_all_pos+corr_all_pos')+eye(size(corr_all_pos));
            corr_all_pos_sig(isnan(corr_all_pos_sig)) = 0;
            corr_all_pos_sig = 0.5*(corr_all_pos_sig+corr_all_pos_sig')+eye(size(corr_all_pos_sig));


    end

    % Get edge weight
    g_desynch = graph(corr_all_pos_sig,'omitselfloops');
    inactive_cells = unique([find(sum(c_run)==1),find(sum(c_rest)==1)]);
    g_desynch = rmnode(g_desynch,inactive_cells);
    xCord = cen_rest(:,1);
    yCord = cen_rest(:,2);
    xCord(inactive_cells) = [];
    yCord(inactive_cells) = [];
    %     LWidths = 0.75*e.Edges.Weight/max(e.Edges.Weight);



    % plot network with physical roi location
    figure('Units','Normalized','Position',[0.25 0.25 0.5 0.5])
    %     subplot(1,2,1)
    g1 = plot(g_desynch,'Layout','force','WeightEffect','direct','NodeColor',[.3185,.3185,.3185],'EdgeColor',[0.4196 0.4196 0.4196],'MarkerSize',4,'NodeLabel', {},'EdgeAlpha',1); %,'LineWidth',LWidths,,'NodeLabel',cen_rest(:,3),'EdgeLabel',round(g_rest.Edges.Weight,1));
    nodeColors = results(row).heat_desynch;
    nodeColors(inactive_cells) = [];
    g1.NodeCData = nodeColors;
    colormap(rainbowmap)
    axis square
    colorbar
    caxis([-0.1 0.1])
    title([tt, sprintf('%s day %i Desynchrony',results(row).animal,results(row).day)])
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs\maps\divergent color_052223",sprintf('%s day %i physical %s.png',results(row).animal,results(row).day,tt{:})));
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs\maps\divergent color_052223",sprintf('%s day %i physical %s.svg',results(row).animal,results(row).day,tt{:})));
%     saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs\maps\divergent color_052223",sprintf('%s day %i physical %s.fig',results(row).animal,results(row).day,tt{:})));
    saveas(gcf,fullfile("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\j_graphs\maps\divergent color_052223",sprintf('%s day %i force %s.pdf',results(row).animal,results(row).day,tt{:})));



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


%