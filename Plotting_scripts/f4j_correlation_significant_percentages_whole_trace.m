%% Description

% Plots of cells correlated to each other marked by their physical location
% only positive correlation plotted
%% Initialize

addpath('J:\nexmif_paper\Utils')  % Add utilities
init % Initialize data directories and genotypes
% close all

%% Begin code
t_start = tic;
% load correlation stats
% load(fullfile('D:\nexmif_paper\code_ball\stats\correlation','corrStatsTableSpeedAsym'),'corrStatsTable')
corr = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTablePearson_09_12'),'corrStatsTable');
corrStatsTable = corr.corrStatsTable;
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];


load('J:\nexmif_paper\code_ball\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells_table = struct2table(responsive_cells);

miceBad = {};
list_sig = [];
list_non_sig = [];
listisKO = [];
listEs = [];
list_sig_mov = [];
list_mov = [];

for k =  1:numel(corrStatsTable.animal)
    mouseName = corrStatsTable{k,1}{:};
    day = corrStatsTable{k,2}{:};
    
    mPath = fullfile(dataFolder,sprintf('fullData_%s_D%i_%s',mouseName,day,'ball'));
    
    % handle missing files
    try
        data_var = load(mPath);
        fullData = data_var.fullData;
    catch
        miceBad = [miceBad,[mouseName,day,'ball']];
    end
    
    % display
    {mouseName, day, 'ball'}
    
    % get cell centroids
    
    goodTraces = fullData.goodIdx;
    goodTraces = setdiff(goodTraces,fullData.empty_idx_minusBG_new);
    centroids = fullData.centroids;
    centroids =  centroids(goodTraces,:);
    
    % get and threshold significant correlations
    
    n = numel(goodTraces);
    y_close = double(corrStatsTable{k,7}{:}<20);
    y_close(logical(triu(ones(n)))) = nan;
    
    w_corr = corrStatsTable{k,4}{:};
    w_corr(logical(triu(ones(n)))) = nan;
    w_thresh_pos = corrStatsTable{k,5}{:};
    w_thresh_neg = corrStatsTable{k,6}{:};
    adj_sig_corr_pos = (w_corr - w_thresh_pos > 1e-5) & (y_close == 0);
    idx_pos = find(adj_sig_corr_pos);
    adj_sig_corr_neg = w_corr - w_thresh_neg < 1e-5 & (y_close == 0);
    idx_neg = find(adj_sig_corr_neg);
    adj_sig_corr_non_sig = (w_corr < w_thresh_pos) & (y_close == 0);
    
    %% Mov responsive
    isM = strcmp({responsive_cells.animal},mouseName);
    isD = [responsive_cells.day] == day;
    isR = [responsive_cells.responsive_bout];
    cell_list = [responsive_cells.cell];
    isR(isnan(isR)) = 0;
    [~,res_cell_list_idx] = ismember(cell_list(isM & isD & isR),goodTraces);
    
    % Both cells have to be mov responsive
    res_cell_mat_both = zeros(n);
    p = nchoosek(res_cell_list_idx,2);
    res_cell_mat_both(sub2ind([n,n],p(:,1),p(:,2))) = 1;
    res_cell_mat_both(sub2ind([n,n],p(:,2),p(:,1))) = 1;
    
    %% Stats

    list_sig = [list_sig,sum(adj_sig_corr_pos(:))];  
    list_non_sig = [list_non_sig,sum(adj_sig_corr_non_sig(:))];
  
    listEs = [listEs, n*(n-1)/2 - sum(y_close(:) == 1)];
    listisKO = [listisKO,ismember(mouseName,miceKO)];
    

    list_sig_mov = [list_sig_mov,sum(sum(adj_sig_corr_pos & res_cell_mat_both))];
    list_mov = [list_mov,sum(sum(res_cell_mat_both))/2];

end
%%
listisKO = logical(listisKO);
save(fullfile('J:\nexmif_paper\code_ball\stats\correlation\','identity_stats_pearson_whole_trace_09_23'),...
   'list_sig','list_sig_mov','listEs','list_mov','list_non_sig','listisKO')

%% Fisher - cellwise

list_sig
ff = 1;
fisher_struct = struct();

fisher_struct(ff).type = 'WT vs KO';
fisher_struct(ff).fisherMat = zeros(2);
fisher_struct(ff).fisherMat(1,1) = sum(list_sig(~listisKO));
fisher_struct(ff).fisherMat(1,2) = sum(listEs(~listisKO))-sum(list_sig(~listisKO));
fisher_struct(ff).fisherMat(2,1) = sum(list_sig(listisKO));
fisher_struct(ff).fisherMat(2,2) = sum(listEs(listisKO))-sum(list_sig(listisKO));


[fisher_struct(ff).h_f,...
    fisher_struct(ff).p_f,...
    fisher_struct(ff).stats_f] = fishertest(fisher_struct(ff).fisherMat,'Alpha',0.05);%,'Tail','right'

fisher_struct(ff).pr_1 = fisher_struct(ff).fisherMat(1,1)/(fisher_struct(ff).fisherMat(1,1)+ fisher_struct(ff).fisherMat(1,2))*100;
fisher_struct(ff).SEP_1 = 1.96*sqrt(fisher_struct(ff).pr_1*(100-fisher_struct(ff).pr_1)/(fisher_struct(ff).fisherMat(1,1)+fisher_struct(ff).fisherMat(1,2)));
fisher_struct(ff).pr_2 = fisher_struct(ff).fisherMat(2,1)/(fisher_struct(ff).fisherMat(2,1)+ fisher_struct(ff).fisherMat(2,2))*100;
fisher_struct(ff).SEP_2 = 1.96*sqrt(fisher_struct(ff).pr_2*(100-fisher_struct(ff).pr_2)/(fisher_struct(ff).fisherMat(2,1)+ fisher_struct(ff).fisherMat(2,2)));




figure
b = bar([1,1.5],[fisher_struct(1).pr_1,fisher_struct(1).pr_2]);
b.FaceColor = 'flat';
b.EdgeColor = 'k';
b.LineWidth = 0.75;
b.BarWidth = 0.5;
b.CData = colors([2,4],:);
hold on
xticks(gca,[1,1.5])
xticklabels(gca, {'WT','KO'})
er = errorbar([1,1.5],[fisher_struct(1).pr_1,fisher_struct(1).pr_2],...
    [fisher_struct(1).SEP_1,fisher_struct(1).SEP_2]);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 0.2;
xlim([0.75,1.75])


ylabel('Significantly correlated cells (%)')
set(gca,'Units','inches','InnerPosition',[.8 .5 0.75 1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',0.5)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_08_23\Codes\plots\pdfs',sprintf('Fig4_whole_trace_fractions.pdf')));

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e1_correlation_significant_percentages',sprintf('Pearson Perc of sig corr out of all edges cellwise.fig')));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e1_correlation_significant_percentages',sprintf('Pearson Perc of sig corr out of all edges cellwise.png')));
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e1_correlation_significant_percentages',sprintf('Pearson Perc of sig corr out of all edges cellwise.epsc')));


