%% Description

% Add mov responsive cells to graph stats table 

%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
addpath(genpath('J:\nexmif_paper\Utils'));
init % Initialize data directories and genotypes
close all

%% Begin code


% load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedAsym'),'corrStatsTable')
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')

fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

%% Load mov response 
load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22')
responsive_cells_table = struct2table(responsive_cells);
%% gather 
load(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','grpahStatsTable_08_26'))
graphWeights = table2struct(graphStatsTable);
for row = 1:size(corrStatsTable,1)
     % Distance between cells
    cell_dis = corrStatsTable{row,10}{:};
    cell_close = double(cell_dis<20);
   
    % resting correlation
    corr_rest  = corrStatsTable{row,4}{:};
    thresh_rest  = corrStatsTable{row,5}{:};
    n = size(corr_rest);
    cell_close(logical(triu(ones(n)))) = nan;
    corr_rest(logical(triu(ones(n)))) = nan;
    thresh_rest(logical(triu(ones(n)))) = nan;
    
    corr_rest_pos = corr_rest.*(corr_rest>0);
    sig_edges_mask_rest = (corr_rest-thresh_rest >1e-5) & (cell_close == 0);
    corr_rest_pos_sig = corr_rest_pos.*sig_edges_mask_rest;
    
    % Make Symmetric
    corr_rest_pos(isnan(corr_rest_pos)) = 0;
    corr_rest_pos = 0.5*(corr_rest_pos+corr_rest_pos')+eye(size(corr_rest_pos));
    corr_rest_pos_sig(isnan(corr_rest_pos_sig)) = 0;
    corr_rest_pos_sig = 0.5*(corr_rest_pos_sig+corr_rest_pos_sig')+eye(size(corr_rest_pos_sig));
    
    % running
    corr_run  = corrStatsTable{row,7}{:};
    thresh_run  = corrStatsTable{row,8}{:};
    n = size(corr_run);
    corr_run(logical(triu(ones(n)))) = nan;
    thresh_run(logical(triu(ones(n)))) = nan;
    
    corr_run_pos = corr_run.*(corr_run>0);
    sig_edges_mask_run = (corr_run-thresh_run >1e-5) & (cell_close == 0);
    corr_run_pos_sig = corr_run_pos.*sig_edges_mask_run;
    
    corr_run_pos(isnan(corr_run_pos)) = 0;
    corr_run_pos = 0.5*(corr_run_pos+corr_run_pos')+eye(size(corr_run_pos));
    corr_run_pos_sig(isnan(corr_run_pos_sig)) = 0;
    corr_run_pos_sig = 0.5*(corr_run_pos_sig+corr_run_pos_sig')+eye(size(corr_run_pos_sig));
    %% Add mov res
    mouse = corrStatsTable{row,1}{:};
    day = corrStatsTable{row,2}{:};
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_ball',mouse,day));
    
    % handle missing files
    try
        load(mPath)
    catch
        miceBad = [miceBad,[miceStudy(m);dIdx(d);{lower(conditionList{c})}]];
    end
    
    % display
    {mouse, day, 'ball'}
    
    % collect traces
    goodTraces = fullData.goodIdx;
    goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);
    binTraces = binarizeTrace(fullData.roi_list_minusBG_new);
    binTraces = binTraces(goodTraces,:);
    binTracesRun = binTraces(:,fullData.movBoutIdx);
    binTracesRest = binTraces(:,fullData.restBoutIdx);
    
    % Add mov responsive 
   
    isM = strcmp(responsive_cells_table.animal,{mouse});
    isD = [responsive_cells_table.day] == day;
    cell_list = [responsive_cells_table.cell];

    isRun = [responsive_cells_table.responsive_bout];
    isRun(isnan(isRun)) = 0;
    [mov_run_list,~] = ismember(goodTraces,cell_list(isM & isD & isRun));
        inactiveRun_sig = unique([find(sum(binTracesRun,2) == 0);find(sum(corr_run_pos_sig,2)==1)]);
    
    mov_run_list(inactiveRun_sig) = [];

    isRest = [responsive_cells_table.responsive_rest];
    isRest(isnan(isRest)) = 0;
    [mov_rest_list,~] = ismember(goodTraces,cell_list(isM & isD & isRest));
    inactiveRest_sig = unique([find(sum(binTracesRest,2) == 0);find(sum(corr_rest_pos_sig,2)==1)]);
    mov_rest_list(inactiveRest_sig) = [];


    graphWeights(row).mov_rest_list_sig = mov_rest_list;
    graphWeights(row).mov_run_list_sig = mov_run_list;
    
    
end

graphStatsTable = struct2table(graphWeights);
save(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','grpahStatsTable_10_17'),'graphStatsTable')

% %% Compare metrics;
% 
% load(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','grpahStatsTable_08_26'),'graphStatsTable')
% 
% results = struct();
% for corrType = 1:3  % 1 - Similarity, 2 - All corr, 3 - Sig corr
%     
%     for row =  1:size(graphStatsTable,1)
%         switch corrType
%             case 1
%                 c_rest = graphStatsTable{row,4}{:};
%                 c_run = graphStatsTable{row,5}{:};
%             case 2
%                 c_rest = graphStatsTable{row,8}{:};
%                 c_run = graphStatsTable{row,9}{:};
%             case 3
%                 c_rest = graphStatsTable{row,12}{:};
%                 c_run = graphStatsTable{row,13}{:};
%         end
%         
%         d_rest = sqrt(-log(c_rest+eye(size(c_rest))));
%         d_rest = d_rest-diag(diag(d_rest));
%         g_rest = graph(d_rest,'omitselfloops');
%         g_rest = rmedge(g_rest,find(isinf(g_rest.Edges{:,2})));
%         wcc_rest = size(c_rest,1)*centrality(g_rest,'closeness','Cost',g_rest.Edges.Weight);
%         [~,i_rest] = min((wcc_rest)*100);
%         [m_rest] = mean((wcc_rest)*100);
%         
%         
%         d_run = sqrt(-log(c_run+eye(size(c_run))));
%         d_run = d_run-diag(diag(d_run));
%         g_run = graph(d_run,'omitselfloops');
%         % Add centroid here
%         
%         g_run = rmedge(g_run,find(isinf(g_run.Edges{:,2})));
%         wcc_run = size(c_run,1)*centrality(g_run,'closeness','Cost',g_run.Edges.Weight);
%         [~ ,i_run]= min((wcc_run)*100);
%         [m_run]= mean((wcc_run)*100);
%           
%           
%         % Assign session details
%         results(row).animal = graphStatsTable{row,1}{:};
%         results(row).day = graphStatsTable{row,2};
%         results(row).condition = graphStatsTable{row,3}{:};
%         
%         switch corrType
%             case 1
%                 results(row).sim_m_rest = m_rest;
%                 results(row).sim_m_run = m_run;
%             case 2
%                 results(row).corr_m_rest = m_rest;
%                 results(row).corr_m_run = m_run;
%             case 3
%                 results(row).corr_sig_m_rest = m_rest;
%                 results(row).corr_sig_m_run = m_run;
%         end
%         
%         switch corrType
%             case 1
%                 results(row).sim_g_rest = g_rest;
%                 results(row).sim_g_run = g_run;
%             case 2
%                 results(row).corr_g_rest = g_rest;
%                 results(row).corr_g_run = g_run;
%             case 3
%                 results(row).corr_sig_g_rest = g_rest;
%                 results(row).corr_sig_g_run = g_run;
%         end
%         
%        switch corrType
%             case 1
%                 results(row).sim_i_rest = i_rest;
%                 results(row).sim_i_run = i_run;
%             case 2
%                 results(row).corr_i_rest = i_rest;
%                 results(row).corr_i_run = i_run;
%             case 3
%                 results(row).corr_sig_i_rest = i_rest;
%                 results(row).corr_sig_i_run = i_run;
%         end
%     end
% end
% 
% %% Plot
% close all
% resultsTable = struct2table(results);
% for corrType = 1:3
%     
%     switch corrType
%         case 1
%             cc_WT_rest = resultsTable{21:end,4};
%             cc_KO_rest = resultsTable{1:20,4};
%             cc_WT_run = resultsTable{21:end,5};
%             cc_KO_run = resultsTable{1:20,5};
%             tt = {'Closeness - Cosine similarity'};
%         case 2
%             cc_WT_rest = resultsTable{21:end,10};
%             cc_KO_rest = resultsTable{1:20,10};
%             cc_WT_run = resultsTable{21:end,11};
%             cc_KO_run = resultsTable{1:20,11};
%             tt = {'Closeness - Pearson correlation'};
%         case 3
%             cc_WT_rest = resultsTable{21:end,16};
%             cc_KO_rest = resultsTable{1:20,16};
%             cc_WT_run = resultsTable{21:end,17};
%             cc_KO_run = resultsTable{1:20,17};
%             tt = {'Closeness - Pearson correlation sig'};
%     end
%     
%     
%     % Scatter
%     
%     figure('units','normalized','position',[0.2,0.4,0.5,0.4])
%     subplot(1,2,1)
%     scatter(cc_WT_rest,cc_WT_run,'b','filled')
%     hold on
%     scatter(cc_KO_rest,cc_KO_run,'r','filled')
%     legend('WT','KO')
%     xlabel('Rest')
%     ylabel('Run')
%     title(tt)
%     hold on
%     mm = max([cc_WT_rest;cc_WT_run;cc_KO_rest;cc_KO_run]);
%     nn = min([cc_WT_rest;cc_WT_run;cc_KO_rest;cc_KO_run]);
%     plot([nn,mm],[nn,mm],'k');
%     legend('WT','KO','Run  = Rest','Location','SouthWest')
%     axis square
%     % Boxplot
%     
%     subplot(1,2,2)
%     b1 = boxchart(1*ones(1,numel(cc_WT_rest)),...
%         cc_WT_rest,...
%         'BoxFaceAlpha',0,...
%         'MarkerStyle','none',...
%         'BoxFaceColor',colors(1,:),...
%         'WhiskerLineColor', colors(1,:),'LineWidth',2);
%     hold on
%     s1 = scatter(1,cc_WT_rest,48,colors(1,:),'filled');
%     
%     b2 = boxchart(2*ones(1,numel(cc_WT_run)),...
%         cc_WT_run,...
%         'BoxFaceAlpha',0,...
%         'MarkerStyle','none',...
%         'BoxFaceColor',colors(2,:),...
%         'WhiskerLineColor', colors(2,:),'LineWidth',2);
%     hold on
%     s2 = scatter(2,cc_WT_run,48,colors(2,:),'filled');
%     
%     b3 = boxchart(3*ones(1,numel(cc_KO_rest)),...
%         cc_KO_rest,...
%         'BoxFaceAlpha',0,...
%         'MarkerStyle','none',...
%         'BoxFaceColor',colors(3,:),...
%         'WhiskerLineColor', colors(3,:),'LineWidth',2);
%     hold on
%     s3 = scatter(3,cc_KO_rest,48,colors(3,:),'filled');
%     
%     b4 = boxchart(4*ones(1,numel(cc_KO_run)),...
%         cc_KO_run,...
%         'BoxFaceAlpha',0,...
%         'MarkerStyle','none',...
%         'BoxFaceColor',colors(4,:),...
%         'WhiskerLineColor', colors(4,:),'LineWidth',2);
%     hold on
%     s4 = scatter(4,cc_KO_run,48,colors(4,:),'filled');
%     
%     xticks([1,2,3,4])
%     xticklabels({'WT Rest','WT Run','KO Rest','KO Run'})
%     [hrest, prest] = ttest2(cc_WT_rest,cc_KO_rest);
%     [hrun, prun] = ttest2(cc_WT_run,cc_KO_run);
%     [hWT, pWT] = ttest(cc_WT_rest,cc_WT_run);
%     [hKO, pKO] = ttest(cc_KO_rest,cc_KO_run);
%     yy = ylim;
%     text(1.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pWT),'HorizontalAlignment','center')
%     % text(2, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prest),'HorizontalAlignment','center')
%     % text(3, yy(2)*1.05, sprintf('P-ttest: %.3f' ,prun),'HorizontalAlignment','center')
%     text(3.5, yy(2)*0.95, sprintf('P-ttest: %.3f' ,pKO),'HorizontalAlignment','center')
%     title(tt')
%     axis square
% end
% 
% %% Plot some graphs
% close all
% for row = 1:32
%     for corrType = 1:3
%         switch corrType
%             case 1
%                 g_rest = results(row).sim_g_rest;
%                 g_run =  results(row).sim_g_run;
%                 tt = {'Graph - Cosine similarity'};
%                 cen_rest = graphStatsTable{row,6}{:};
%                 cen_run = graphStatsTable{row,7}{:};
%                 cc_rest = results(row).sim_m_rest;
%                 cc_run = results(row).sim_m_run;
%                 ii_rest = results(row).sim_i_rest;
%                 ii_run = results(row).sim_i_run;
%             case 2
%                 g_rest = results(row).corr_g_rest;
%                 g_run = results(row).corr_g_run;
%                 tt = {'Graph - Pearson correlation'};
%                 cen_rest = graphStatsTable{row,10}{:};
%                 cen_run = graphStatsTable{row,11}{:};
%                 cc_rest = results(row).corr_m_rest;
%                 cc_run = results(row).corr_m_run;
%                 ii_rest = results(row).corr_i_rest;
%                 ii_run = results(row).corr_i_run;
%             case 3
%                 g_rest = results(row).corr_sig_g_rest;
%                 g_run = results(row).corr_sig_g_run;
%                 tt = {'Graph - Pearson correlation sig'};
%                 cen_rest = graphStatsTable{row,14}{:};
%                 cen_run = graphStatsTable{row,15}{:};
%                 cc_rest = results(row).corr_sig_m_rest;
%                 cc_run = results(row).corr_sig_m_run;
%                 ii_rest = results(row).corr_sig_i_rest;
%                 ii_run = results(row).corr_sig_i_run;
%         end
%         
%         figure('Units','Normalized','Position',[0 0 1 1])
%         subplot(1,2,1)
%         g1 = plot(g_rest,'EdgeAlpha',0.75,'NodeColor','k','EdgeColor','#7E2F8E','XData',cen_rest(:,1),'YData',cen_rest(:,2),'NodeLabel',cen_rest(:,3),'MarkerSize',6);%,'EdgeLabel',round(g_rest.Edges.Weight,1));
%         axis square
%         title([tt, sprintf('%s day %i Rest',results(row).animal,results(row).day),sprintf('Mean Closeness Centrality - %.2f',cc_rest)])
%         subplot(1,2,2)
%         g2 = plot(g_run,'EdgeAlpha',0.75,'NodeColor','k','EdgeColor','#77AC30','XData',cen_run(:,1),'YData',cen_run(:,2),'NodeLabel',cen_run(:,3),'MarkerSize',6);%,'EdgeLabel',round(g_rest.Edges.Weight,1));
%         axis square
%         title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day),sprintf('Mean Closeness Centrality - %.2f',cc_run)])
%         saveas(gcf,fullfile("J:\nexmif_paper\code_ball\plots\Graph\Networks",sprintf('%s day %i Physical %s.png',results(row).animal,results(row).day,tt{:})));
%         
%         figure('Units','Normalized','Position',[0 0 1 1])
%         subplot(1,2,1)
%         NodeCData = repmat([0,0,0],size(g_rest.Nodes,1),1);
%         NodeCData(ii_rest,:) = [1,0,0];
%         g1 = plot(g_rest,'Layout','force','WeightEffect','direct','EdgeAlpha',0.5,'NodeColor', NodeCData,'NodeLabel',cen_rest(:,3),'MarkerSize',6); %,'EdgeLabel',round(g_run.Edges.Weight,1));
%         axis square
%         
%         title([tt, sprintf('%s day %i Rest',results(row).animal,results(row).day),sprintf('Mean Closeness Centrality - %.2f',cc_rest)])
%         subplot(1,2,2)
%         NodeCData = repmat([0,0,0],size(g_run.Nodes,1),1);
%         NodeCData(ii_run,:) = [1,0,0];
%         g2 = plot(g_run,'Layout','force','WeightEffect','direct','EdgeAlpha',0.5,'NodeColor',NodeCData,'NodeLabel',cen_run(:,3),'MarkerSize',6); %,'EdgeLabel',round(g_run.Edges.Weight,1));
%         axis square
%         
%         title([tt, sprintf('%s day %i Run',results(row).animal,results(row).day),sprintf('Mean Closeness Centrality - %.2f',cc_run)])
%         saveas(gcf,fullfile('J:\nexmif_paper\code_ball\plots\Graph\Networks',sprintf('%s day %i Schematic %s.png',results(row).animal,results(row).day,tt{:})));
%                 
%         close all
%     end
% end
