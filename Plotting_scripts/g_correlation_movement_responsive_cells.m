%% Description

% What happens to significantly correlated pairs when states change % speed low -> speed  high -- corr
% decreases and vice versa
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils'))  % Add utilities
init % Initialize data directories and genotypes
% close all

%% Begin code
t_start = tic;
% load correlation stats
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells_table = struct2table(responsive_cells);

miceBad = {};
listisKO = [];

list_run_sig = [];
list_rest_sig = [];
list_both_sig = [];
list_run_non_sig = [];
list_rest_non_sig = [];
list_both_non_sig = [];
listEs = [];

list_run_sig_mov = [];
list_run_non_sig_mov = [];
list_mov = [];
list_run_sig_mov_val_WT = [];
list_run_non_sig_mov_val_WT = [];
list_run_sig_mov_val_KO = [];
list_run_non_sig_mov_val_KO = [];

list_run_sig_non_mov = [];
list_run_non_sig_non_mov = [];
list_non_mov = [];
list_run_sig_non_mov_val_WT = [];
list_run_non_sig_non_mov_val_WT = [];
list_run_sig_non_mov_val_KO = [];
list_run_non_sig_non_mov_val_KO = [];

list_run_sig_mov_one =[];
list_run_non_sig_mov_one =[];
list_mov_one = [];
list_run_sig_mov_val_one_WT = [];
list_run_non_sig_mov_val_one_WT = [];
list_run_sig_mov_val_one_KO = [];
list_run_non_sig_mov_val_one_KO = [];
list_n_mov_cells = [];

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
    y_close = double(corrStatsTable{k,10}{:}<20);
    y_close(logical(triu(ones(n)))) = nan;
    
    w_corr_low = corrStatsTable{k,4}{:};
    %     w_corr_low(isnan(w_corr_low)) = 0;
    w_corr_low(logical(triu(ones(n)))) = nan;
    w_thresh_low_pos = corrStatsTable{k,5}{:};
    w_thresh_low_neg = corrStatsTable{k,6}{:};
    adj_sig_corr_low_pos = (w_corr_low - w_thresh_low_pos > 1e-5) & (y_close == 0);
    idx_low_pos = find(adj_sig_corr_low_pos);
    adj_sig_corr_low_neg = w_corr_low - w_thresh_low_neg < 1e-5 & (y_close == 0);
    idx_low_neg = find(adj_sig_corr_low_neg);
    adj_sig_corr_low_non_sig = (w_corr_low - w_thresh_low_pos < 1e-5) & (w_corr_low - w_thresh_low_neg > 1e-5) & (y_close == 0);
    %     adj_sig_corr_low_non_sig = (w_corr_low < w_thresh_low_pos) & (y_close == 0);
    
    w_corr_high = corrStatsTable{k,7}{:};
    %     w_corr_high(isnan(w_corr_high)) = 0;
    w_corr_high(logical(triu(ones(n)))) = nan;
    w_thresh_high_pos = corrStatsTable{k,8}{:};
    w_thresh_high_neg = corrStatsTable{k,9}{:};
    adj_sig_corr_high_pos = (w_corr_high - w_thresh_high_pos > 1e-5) & (y_close == 0);
    idx_high_pos = find(adj_sig_corr_high_pos);
    adj_sig_corr_high_neg = (w_corr_high - w_thresh_high_neg < 1e-5) & (y_close == 0);
    idx_high_neg = find(adj_sig_corr_high_neg);
    adj_sig_corr_high_non_sig = (w_corr_high - w_thresh_high_pos < 1e-5) & (w_corr_high - w_thresh_high_neg > 1e-5) & (y_close == 0);
    %      adj_sig_corr_high_non_sig = (w_corr_high < w_thresh_high_pos)  & (y_close == 0);
    idx_both_pos = intersect(idx_low_pos,idx_high_pos);
    
    
    % Mov responsive
    
    isM = strcmp({responsive_cells.animal},mouseName);
    isD = [responsive_cells.day] == day;
    isR = [responsive_cells.responsive_bout];
    cell_list = [responsive_cells.cell];
    isR(isnan(isR)) = 0;
    [~,res_cell_list_idx] = ismember(cell_list(isM & isD & isR),goodTraces);
    
    p = nchoosek(res_cell_list_idx,2);  % Both cells have to be mov responsive
    res_cell_mat_both = zeros(n);
    res_cell_mat_both(sub2ind([n,n],p(:,1),p(:,2))) = 1;
    res_cell_mat_both(sub2ind([n,n],p(:,2),p(:,1))) = 1;
    
    % Only one cell is mov responsive
    res_cell_mat_one = zeros(n);
    res_cell_mat_one(res_cell_list_idx,:) = 1;
    res_cell_mat_one(:,res_cell_list_idx) = 1;
    res_cell_mat_one(logical(res_cell_mat_both)) = 0;
    
    % count movement modulated cells
    
    list_n_mov_cells = [list_n_mov_cells ;[numel(res_cell_list_idx),n]];
    
    % is WT
    listisKO = [listisKO,ismember(mouseName,miceKO)];
    %%% All cells
    list_run_sig = [list_run_sig,sum(adj_sig_corr_high_pos(:))];
    list_rest_sig = [list_rest_sig,sum(adj_sig_corr_low_pos(:))];
    list_both_sig = [list_both_sig,sum(adj_sig_corr_low_pos(:) & adj_sig_corr_high_pos(:))];
    
    list_run_non_sig = [list_run_non_sig,sum(adj_sig_corr_high_non_sig(:))];
    list_rest_non_sig = [list_rest_non_sig,sum(adj_sig_corr_low_non_sig(:))];
    list_both_non_sig = [list_both_non_sig,sum(adj_sig_corr_low_non_sig(:) & adj_sig_corr_high_non_sig(:))];
    
    listEs = [listEs, n*(n-1)/2 - sum(y_close(:) == 1)];
    
    %%% Mov Modulated cells only - both cells
    
    % count of sig pairs
    list_run_sig_mov = [list_run_sig_mov,sum(sum(adj_sig_corr_high_pos & res_cell_mat_both))];
    list_run_non_sig_mov = [list_run_non_sig_mov,sum(sum(adj_sig_corr_high_non_sig & res_cell_mat_both))];
    list_mov = [list_mov,sum(sum(res_cell_mat_both))/2];
    
    % collect all values
    if ismember(mouseName,miceWT)
        list_run_sig_mov_val_WT = [list_run_sig_mov_val_WT;w_corr_high(adj_sig_corr_high_pos & res_cell_mat_both)];
        list_run_non_sig_mov_val_WT = [list_run_non_sig_mov_val_WT;w_corr_high(adj_sig_corr_high_non_sig & res_cell_mat_both)];
    else
        list_run_sig_mov_val_KO = [list_run_sig_mov_val_KO;w_corr_high(adj_sig_corr_high_pos & res_cell_mat_both)];
        list_run_non_sig_mov_val_KO = [list_run_non_sig_mov_val_KO;w_corr_high(adj_sig_corr_high_non_sig & res_cell_mat_both)];
    end
    %%% Non mov modulated cells only - both cells
    % count of sig pairs
    list_run_sig_non_mov = [list_run_sig_non_mov,sum(sum(adj_sig_corr_high_pos & ~res_cell_mat_both))];
    list_run_non_sig_non_mov = [list_run_non_sig_non_mov,sum(sum(adj_sig_corr_high_non_sig & ~res_cell_mat_both))];
    list_non_mov = [list_non_mov,sum(sum(~res_cell_mat_both))/2];
    
    % collect all values
    if ismember(mouseName,miceWT)
        list_run_sig_non_mov_val_WT = [list_run_sig_non_mov_val_WT;w_corr_high(adj_sig_corr_high_pos & ~res_cell_mat_both)];
        list_run_non_sig_non_mov_val_WT = [list_run_non_sig_non_mov_val_WT;w_corr_high(adj_sig_corr_high_non_sig & ~res_cell_mat_both)];
    else
        list_run_sig_non_mov_val_KO = [list_run_sig_non_mov_val_KO;w_corr_high(adj_sig_corr_high_pos & ~res_cell_mat_both)];
        list_run_non_sig_non_mov_val_KO = [list_run_non_sig_non_mov_val_KO;w_corr_high(adj_sig_corr_high_non_sig & ~res_cell_mat_both)];
        
    end
    %%% Mov Modulated cells only - only one cell
    % count of sig pairs
    list_run_sig_mov_one = [list_run_sig_mov_one,sum(sum(adj_sig_corr_high_pos & res_cell_mat_one))];
    list_run_non_sig_mov_one = [list_run_non_sig_mov_one,sum(sum(adj_sig_corr_high_non_sig & res_cell_mat_one))];
    list_mov_one = [list_mov_one,sum(sum(res_cell_mat_one))/2];
    
    % collect all values
    if ismember(mouseName,miceWT)
        list_run_sig_mov_val_one_WT = [list_run_sig_mov_val_one_WT;w_corr_high(adj_sig_corr_high_pos & res_cell_mat_one)];
        list_run_non_sig_mov_val_one_WT = [list_run_non_sig_mov_val_one_WT;w_corr_high(adj_sig_corr_high_non_sig & res_cell_mat_one)];
    else
        list_run_sig_mov_val_one_KO = [list_run_sig_mov_val_one_KO;w_corr_high(adj_sig_corr_high_pos & res_cell_mat_one)];
        list_run_non_sig_mov_val_one_KO = [list_run_non_sig_mov_val_one_KO;w_corr_high(adj_sig_corr_high_non_sig & res_cell_mat_one)];
    end
end

listisKO = logical(listisKO);

%% Plot

% compare Modulated to non modulated
frac_run_sig_mov_WT = list_run_sig_mov(~listisKO)./list_mov(~listisKO)*100;
frac_run_sig_mov_KO = list_run_sig_mov(listisKO)./list_mov(listisKO)*100;
frac_run_sig_non_mov_WT = list_run_sig_non_mov(~listisKO)./list_non_mov(~listisKO)*100;
frac_run_sig_non_mov_KO = list_run_sig_non_mov(listisKO)./list_non_mov(listisKO)*100;

meanMovWT = mean(frac_run_sig_mov_WT);
meanMovKO = mean(frac_run_sig_mov_KO);
meanNonMovWT = mean(frac_run_sig_non_mov_WT);
meanNonMovKO = mean(frac_run_sig_non_mov_KO);

stdMovWT = std(frac_run_sig_mov_WT);
stdMovKO = std(frac_run_sig_mov_KO);
stdNonMovWT = std(frac_run_sig_non_mov_WT);
stdNonMovKO = std(frac_run_sig_non_mov_KO);

figure
boxplot([frac_run_sig_mov_WT';...
    frac_run_sig_non_mov_WT';...
    frac_run_sig_mov_KO';...
    frac_run_sig_non_mov_KO'],...
    [ones(size(frac_run_sig_mov_WT'));...
    2*ones(size(frac_run_sig_non_mov_WT'));...
    3*ones(size(frac_run_sig_mov_KO'));...
    4*ones(size(frac_run_sig_non_mov_KO'))],'Colors','k')
hold on
s1 = scatter(1,frac_run_sig_mov_WT',[],'filled','MarkerFaceColor', colors(2,:));
hold on
s2 = scatter(2,frac_run_sig_non_mov_WT',[],'filled','MarkerFaceColor', colors(1,:));
hold on
s3 = scatter(3,frac_run_sig_mov_KO',[],'filled','MarkerFaceColor', colors(4,:));
hold on
s4 = scatter(4,frac_run_sig_non_mov_KO',[],'filled','MarkerFaceColor', colors(3,:));

xlim([0.5,4.5])
xticks([1,2,3,4])
xticklabels({'MovMod','NonMovMod','MovMod','NonMovMod'})
[h_WT,p_WT] = ttest(frac_run_sig_mov_WT',frac_run_sig_non_mov_WT');
[h_KO,p_KO] = ttest(frac_run_sig_mov_KO',frac_run_sig_non_mov_KO');
[h_Mod,p_Mod] = ttest2(frac_run_sig_mov_WT',frac_run_sig_mov_KO');
[h_NonMod,p_NonMod] = ttest2(frac_run_sig_non_mov_WT',frac_run_sig_non_mov_KO');
ttests = [p_WT,p_KO,p_Mod,p_NonMod]'*4;
means = [mean(frac_run_sig_mov_WT'),mean(frac_run_sig_mov_KO'),mean(frac_run_sig_non_mov_WT'),mean(frac_run_sig_non_mov_KO')]';
stds = [std(frac_run_sig_mov_WT'),std(frac_run_sig_mov_KO'),std(frac_run_sig_non_mov_WT'),std(frac_run_sig_non_mov_KO')]';

yy  = ylim;
text(1.5,1.05*yy(2),sprintf('p = %.4f',p_WT*4),'HorizontalAlignment','Center')
text(3.5,1.05*yy(2),sprintf('p = %.4f',p_KO*4),'HorizontalAlignment','Center')
text(2,1.3*yy(2),sprintf('p = %.4f',p_Mod*4),'HorizontalAlignment','Center')
text(3,1.3*yy(2),sprintf('p = %.4f',p_NonMod*4),'HorizontalAlignment','Center')
ylim([yy(1),yy(2)*1.5])

ylabel({'Correlated cells (% of modulated population)'})
title({'Fraction of signifcantly correlated', 'movement modulated or Non modulated cell pairs'})

% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 40])
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\g_correlation_modulated_cells\8.3.22',sprintf('%s.fig','Mod vs Non Mod')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\g_correlation_modulated_cells\8.3.22',sprintf('%s.png','Mod vs Non Mod')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\g_correlation_modulated_cells\8.3.22',sprintf('%s.epsc','Mod vs Non Mod')))

