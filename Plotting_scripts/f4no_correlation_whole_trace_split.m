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

corrSpeed = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable');
corrStatsTableSpeed = corrSpeed.corrStatsTable;
% fieldNames = corrStatsTableSpeed.Properties.VariableNames;
corrStatsTableSpeed([6,13,18:19,29,31:35,40:43,45],:) = [];

corr = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTablePearson_09_12'),'corrStatsTable');
corrStatsTable = corr.corrStatsTable;
% fieldNames = corrStatsTableSpeed.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells_table = struct2table(responsive_cells);

miceBad = {};
list_isKO = [];
list_sig = [];
list_sig_res_res = [];
list_sig_res_non = [];
list_sig_non_non = [];

% Values of these groups in rest
list_sig_res_res_val_rest = [];
list_sig_res_non_val_rest = [];
list_sig_non_non_val_rest = [];
list_sig_all_val_rest =  [];

% Values of these groups in run

list_sig_res_res_val_run = [];
list_sig_res_non_val_run = [];
list_sig_non_non_val_run  = [];
list_sig_all_val_run  = [];

list_n_mov_cells = [];

for row =  1:numel(corrStatsTable.animal)
    mouseName = corrStatsTable{row,1}{:};
    day = corrStatsTable{row,2}{:};

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
    % Correlation of cells significant in full duration
    corrMat = corrStatsTable{row,4}{:};
    threshMat_pos = corrStatsTable{row,5}{:};
    disMat_close = double(corrStatsTable{row,7}{:}<20);

    n = size(corrMat);
    corrMat(logical(triu(ones(n)))) = nan;
    threshMat_pos(logical(triu(ones(n)))) = nan;
    disMat_close(logical(triu(ones(n)))) = nan;

    m_sig = (corrMat-threshMat_pos >1e-5) & (disMat_close == 0); % binary mask of cell pairs correlated over the whole session
    m_non_sig = (corrMat-threshMat_pos <= 1e-5) & (disMat_close == 0);


    % Create masks for mvement responsive groups

    % Mov responsive
    isM = strcmp({responsive_cells.animal},mouseName);
    isD = [responsive_cells.day] == day;
    isR = [responsive_cells.responsive_bout];
    cell_list = [responsive_cells.cell];
    isR(isnan(isR)) = 0;
    [~,res_cell_list_idx] = ismember(cell_list(isM & isD & isR),goodTraces);

    % Res - Res
    p = nchoosek(res_cell_list_idx,2);
    m_res_res = zeros(n);
    m_res_res(sub2ind([n,n],p(:,1),p(:,2))) = 1;
    m_res_res(sub2ind([n,n],p(:,2),p(:,1))) = 1;
    close_sym = (disMat_close == 0)+(disMat_close == 0)';
    m_res_res(close_sym == 0) = 0;

    % Res - Non Res
    m_res_non = zeros(n);
    m_res_non(res_cell_list_idx,:) = 1;
    m_res_non(:,res_cell_list_idx) = 1;
    m_res_non(logical(m_res_res)) = 0;
    m_res_non(close_sym == 0) = 0;

    % Non Res - Non Res

    m_non_non = ones(n);
    m_non_non(res_cell_list_idx,:) = 0;
    m_non_non(:,res_cell_list_idx) = 0;
    m_non_non(close_sym == 0) = 0;


    m_res_res(logical(triu(ones(n)))) = nan;
    m_res_non(logical(triu(ones(n)))) = nan;
    m_non_non(logical(triu(ones(n)))) = nan;

    % %% Sanity check plots
    %    figure('WindowStyle','docked')
    %    subplot(2,4,1)
    %    imagesc(m_res_res)
    %    axis square
    %    title('Res-Res')
    %
    %    subplot(2,4,2)
    %    imagesc(m_res_non)
    %    axis square
    %    title('Res-non')
    %
    %    subplot(2,4,3)
    %    imagesc(m_non_non)
    %    axis square
    %    title('non-non')
    %
    %    subplot(2,4,4)
    %    imagesc(m_non_non + m_res_res + m_res_non)
    %    axis square
    %    title('Sum')
    %
    %    subplot(2,4,5)
    %    imagesc((m_sig + m_res_res) == 2)
    %    axis square
    %    title('Res-res sig')
    %
    %    subplot(2,4,6)
    %    imagesc((m_sig + m_res_non) == 2)
    %    axis square
    %    title('Res-non sig')
    %
    %    subplot(2,4,7)
    %    imagesc((m_sig + m_non_non) == 2)
    %    axis square
    %    title('non-non sig')
    %
    %    subplot(2,4,8)
    %    imagesc(m_sig)
    %    axis square
    %    title('All sig')
    %
    %% calculate fractions

    % What percentage of sig edges are in each group
    n_sig = sum(m_sig,'all','omitnan');
    list_sig = [list_sig, n_sig];
    list_sig_res_res = [list_sig_res_res, sum(m_sig+m_res_res == 2,'all','omitnan')/n_sig*100];
    list_sig_res_non = [list_sig_res_non, sum(m_sig+m_res_non == 2,'all','omitnan')/n_sig*100];
    list_sig_non_non = [list_sig_non_non, sum(m_sig+m_non_non == 2,'all','omitnan')/n_sig*100];

    % Values of these groups in rest
    corrMat_rest = corrStatsTableSpeed{row,4}{:};
    corrMat_rest(logical(triu(ones(n)))) = nan;
    list_sig_res_res_val_rest = [list_sig_res_res_val_rest, mean(corrMat_rest(m_sig+m_res_res == 2),'all','omitnan')];
    list_sig_res_non_val_rest = [list_sig_res_non_val_rest, mean(corrMat_rest(m_sig+m_res_non == 2),'all','omitnan')];
    list_sig_non_non_val_rest = [list_sig_non_non_val_rest, mean(corrMat_rest(m_sig+m_non_non == 2),'all','omitnan')];
    list_sig_all_val_rest = [list_sig_non_non_val_rest, mean(corrMat_rest(m_sig == 1),'all','omitnan')];

    % Values of these groups in run
    corrMat_run = corrStatsTableSpeed{row,7}{:};
    corrMat_run(logical(triu(ones(n)))) = nan;
    list_sig_res_res_val_run = [list_sig_res_res_val_run, mean(corrMat_run(m_sig+m_res_res == 2),'all','omitnan')];
    list_sig_res_non_val_run = [list_sig_res_non_val_run, mean(corrMat_run(m_sig+m_res_non == 2),'all','omitnan')];
    list_sig_non_non_val_run = [list_sig_non_non_val_run, mean(corrMat_run(m_sig+m_non_non == 2),'all','omitnan')];
    list_sig_all_val_run = [list_sig_non_non_val_run, mean(corrMat_run(m_sig == 1),'all','omitnan')];

    % count movement modulated cells
    list_n_mov_cells = [list_n_mov_cells ;[numel(res_cell_list_idx),n]];
    % is WT
    list_isKO = [list_isKO,ismember(mouseName,miceKO)];


end

list_isKO = logical(list_isKO);

%% Plot

%% Sessionwise plos of values

% plot 12 groups and do GLM with y - x1*x2*x3

G = ismember(corrStatsTableSpeed{:,1},miceWT);
plot_data =  {list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(G),...
    list_sig_res_res_val_rest(~G),...
    list_sig_res_res_val_run(~G),...
    list_sig_res_non_val_rest(G),...
    list_sig_res_non_val_run(G),...
    list_sig_res_non_val_rest(~G),...
    list_sig_res_non_val_run(~G),...
    list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(G),...
    list_sig_non_non_val_rest(~G),...
    list_sig_non_non_val_run(~G)};

xx_axis = [1,2,3,4,6,7,8,9,11,12,13,14,15];
cc_axis = [1,2,3,4,1,2,3,4,1,2,3,4];
figure
b = struct();
for cc = 1:12
    y_data = plot_data{cc}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    b(cc).bar = boxchart(x_data,y_data,...
        'MarkerStyle','none',...
        'BoxFaceColor',c_data,...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', c_data,...
        'LineWidth',1);
    hold on
    b(cc).scatter = scatter(x_data,y_data,[],'filled','MarkerFaceColor', c_data); %'MarkerEdgeColor', c_data);%
    hold on
    b(cc).mean = scatter(x_data(1),nanmean(y_data),30,'k','*','LineWidth',2);
    hold on
end

% Stats

stat_data = [list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(G),...
    list_sig_res_res_val_rest(~G),...
    list_sig_res_res_val_run(~G),...
    list_sig_res_non_val_rest(G),...
    list_sig_res_non_val_run(G),...
    list_sig_res_non_val_rest(~G),...
    list_sig_res_non_val_run(~G),...
    list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(G),...
    list_sig_non_non_val_rest(~G),...
    list_sig_non_non_val_run(~G)]';
stat_gen = [ones(size(list_sig_res_res_val_rest(G))),...
    ones(size(list_sig_res_res_val_run(G))),...
    zeros(size(list_sig_res_res_val_rest(~G))),...
    zeros(size(list_sig_res_res_val_run(~G))),...
    ones(size(list_sig_res_non_val_rest(G))),...
    ones(size(list_sig_res_non_val_run(G))),...
    zeros(size(list_sig_res_non_val_rest(~G))),...
    zeros(size(list_sig_res_non_val_run(~G))),...
    ones(size(list_sig_non_non_val_rest(G))),...
    ones(size(list_sig_non_non_val_run(G))),...
    zeros(size(list_sig_non_non_val_rest(~G))),...
    zeros(size(list_sig_non_non_val_run(~G)))]';

stat_mov = [zeros(size(list_sig_res_res_val_rest(G))),...
    ones(size(list_sig_res_res_val_run(G))),...
    zeros(size(list_sig_res_res_val_rest(~G))),...
    ones(size(list_sig_res_res_val_run(~G))),...
    zeros(size(list_sig_res_non_val_rest(G))),...
    ones(size(list_sig_res_non_val_run(G))),...
    zeros(size(list_sig_res_non_val_rest(~G))),...
    ones(size(list_sig_res_non_val_run(~G))),...
    zeros(size(list_sig_non_non_val_rest(G))),...
    ones(size(list_sig_non_non_val_run(G))),...
    zeros(size(list_sig_non_non_val_rest(~G))),...
    ones(size(list_sig_non_non_val_run(~G)))]';

stat_res = [2*ones(size(list_sig_res_res_val_rest(G))),...
    2*ones(size(list_sig_res_res_val_run(G))),...
    2*ones(size(list_sig_res_res_val_rest(~G))),...
    2*ones(size(list_sig_res_res_val_run(~G))),...
    ones(size(list_sig_res_non_val_rest(G))),...
    ones(size(list_sig_res_non_val_run(G))),...
    ones(size(list_sig_res_non_val_rest(~G))),...
    ones(size(list_sig_res_non_val_run(~G))),...
    zeros(size(list_sig_non_non_val_rest(G))),...
    zeros(size(list_sig_non_non_val_run(G))),...
    zeros(size(list_sig_non_non_val_rest(~G))),...
    zeros(size(list_sig_non_non_val_run(~G)))]';


% FIT GLM
sessionwise_stats(1).name = '3 factor analysis';
sessionwise_stats(1).GLM.mdl = fitglm([stat_gen,stat_mov,stat_res],stat_data, 'Mean_Corr ~ Genotype*Behavior*Network','Distribution','normal','CategoricalVars',[1,2,3],'VarNames',{'Genotype','Behavior','Network','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl);
sessionwise_stats(1).GLM.fit = deviance_test{2,4};

% sessionwise_stats(1).GLM.mdl = fitlm([stat_gen,stat_mov,stat_res],stat_data, 'y ~ x1*x3','CategoricalVars',[1,2,3]);
figure; interactionplot(stat_data,[stat_gen,stat_mov,stat_res],'varnames',{'Genotype','Behavior','Network'})

% figure; plotInteraction(sessionwise_stats(1).GLM.mdl,'x1','x3','predictions')
% figure; plotInteraction(sessionwise_stats(1).GLM.mdl,'x2','x3','predictions')


% Same information plotted in a different order

xx_axis = [1,2,3,5,6,7,9,10,11,13,14,15];
dd_axis = [1,5,9,2,6,10,3,7,11,4,8,12];
cc_axis = [1,1,1,2,2,2,3,3,3,4,4,4];
figure
b = struct();
for cc = 1:12
    y_data = plot_data{dd_axis(cc)}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    b(cc).bar = boxchart(x_data,y_data,...
        'MarkerStyle','none',...
        'BoxFaceColor',c_data,...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', c_data,...
        'LineWidth',1);
    hold on
    b(cc).scatter = scatter(x_data,y_data,[],'filled','MarkerFaceColor', c_data); %'MarkerEdgeColor', c_data);%
    hold on
    b(cc).mean = scatter(x_data(1),nanmean(y_data),30,'k','*','LineWidth',1);
    hold on
end
%%
% Same information plotted in a different order
figWT = figure;
figKO = figure;
xx_axis = [1,2,4,5,7,8,1,2,4,5,7,8];
dd_axis = [1,2,5,6,9,10,3,4,7,8,11,12];
cc_axis = [1,2,1,2,1,2,3,4,3,4,3,4];
b = struct();
for cc = 1:12
    y_data = plot_data{dd_axis(cc)}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    if cc<=6
        figWT
        tt = 'WT';
        b(cc).bar = boxchart(x_data,y_data,...
            'MarkerStyle','none',...
            'BoxFaceColor',c_data,...
            'BoxFaceAlpha',0,...
            'WhiskerLineColor', c_data,...
            'LineWidth',1);
        hold on
        b(cc).scatter = scatter(x_data,y_data,18,'filled','MarkerFaceColor', c_data); %'MarkerEdgeColor', c_data);%
        hold on
        % b(cc).mean = scatter(x_data(1),nanmean(y_data),30,'k','*','LineWidth',1);
        hold on
        ylim([-0.01,0.25])
        xticks([1.5,4.5,7.5])
        xticklabels({'Both cells \newlinemodulated','One cell\newlinemodulated','Both cells\newlinenonmodulated'})
        ylabel('Correlation coefficient')
        xlabel('Movement Modulation of the cell pair','Position',[4.5,-0.05,1])
        title(tt)
    else
        figKO
        tt = 'KO';
        b(cc).bar = boxchart(x_data,y_data,...
            'MarkerStyle','none',...
            'BoxFaceColor',c_data,...
            'BoxFaceAlpha',0,...
            'WhiskerLineColor', c_data,...
            'LineWidth',1);
        hold on
        b(cc).scatter = scatter(x_data,y_data,18,'filled','MarkerFaceColor', c_data); %'MarkerEdgeColor', c_data);%
        hold on
        % b(cc).mean = scatter(x_data(1),nanmean(y_data),30,'k','*','LineWidth',1);
        hold on
        ylim([-0.01,0.25])
        xticks([1.5,4.5,7.5])
        xticklabels({'Both cells \newlinemodulated','One cell\newlinemodulated','Both cells\newlinenonmodulated'})
        ylabel('Correlation coefficient')
        xlabel('Movement Modulation of the cell pair','Position',[4.5,-0.05,1])
        title(tt)
    end


    % str1 = {'Mod','Mod','Nonmod'};
    % str2 = {'Mod','Nonmod','Nonmod'};
    % labelArray = [str1;str2]
    % tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    % xticklabels({'Mod\newlineMod','Mod\newlineNonmod','Nonmod\newlineNonmod'})
    % xticklabels({'Both cells \newlinemodulated','One cell\newlinemodulated','Both cells\newlinenonmodulated'})
    % ylabel('Correlation coefficient')
    % xlabel('Movement Modulation of the cell pair','Position',[4.5,-0.05,1])
    % title(tt)
end
% subplot(1,2,1)
% set(gca,'OuterPosition',[0,0.1,0.43,0.9])
% legend([b(5).scatter,b(6).scatter],'Rest','Run')
% xlim([0,9])
% set(gca,'InnerPosition',[0.085,0.2,0.35,0.73])
% subplot(1,2,2)
% set(gca,'OuterPosition',[0.45,0.1,0.43,0.9])
% legend([b(11).scatter,b(12).scatter],'Rest','Run')
% set(gca,'InnerPosition',[0.53,0.2,0.35,0.73])
% xlim([0,9])


%% plot for figure
% Same information plotted in a different order
figure;
xx_axis = [1,2,4,5,7,8,1,2,4,5,7,8];
dd_axis = [1,2,5,6,9,10,3,4,7,8,11,12];
cc_axis = [1,2,1,2,1,2,3,4,3,4,3,4];
b = struct();
for cc = 1:12
    y_data = plot_data{dd_axis(cc)}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    if cc<=6
        subplot(1,2,1)
        tt = 'WT';
    else
        subplot(1,2,2)
        tt = 'KO';
    end
    b(cc).bar = boxchart(x_data,y_data,...
        'MarkerStyle','none',...
        'BoxFaceColor',c_data,...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', c_data,...
        'LineWidth',0.5, ...
        'BoxWidth',0.7);
    hold on
    b(cc).scatter = scatter(x_data,y_data,10,'filled','MarkerFaceColor', c_data); %'MarkerEdgeColor', c_data);%
    hold on
    % b(cc).mean = scatter(x_data(1),nanmean(y_data),30,'k','*','LineWidth',1);
    hold on
    ylim([-0.01,0.25])
    xticks([1.5,4.5,7.5])

    % str1 = {'Mod','Mod','Nonmod'};
    % str2 = {'Mod','Nonmod','Nonmod'};
    % labelArray = [str1;str2]
    % tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    % xticklabels({'Mod\newlineMod','Mod\newlineNonmod','Nonmod\newlineNonmod'})
%     xticklabels({'Both cells \newlinemodulated','One cell\newlinemodulated','Both cells\newlinenonmodulated'})
%     ylabel('Correlation coefficient')
%     title(tt)
end
subplot(1,2,1)
% set(gca,'OuterPosition',[0,0.1,0.43,0.9])
xlim([0,9])
% ylim([0,inf])
ylim([0,0.25])
set(gca,'Units','Inches','InnerPosition',[0.5,0.8,1.875,1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
subplot(1,2,2)
% set(gca,'OuterPosition',[0.45,0.1,0.43,0.9])
% set(gca,'InnerPosition',[0.53,0.2,0.35,0.73])
% ylim([0,inf])
ylim([0,0.25])
xlim([0,9])
set(gca,'Units','Inches','InnerPosition',[3,0.8,1.875,1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)


%% Now analyse the networks separately split into 4 groups matching the split 1


g1 = stat_gen== 0 & stat_mov == 0;
g2 = stat_gen== 1 & stat_mov == 0;
g3 = stat_gen== 0 & stat_mov == 1;
g4 = stat_gen== 1 & stat_mov == 1;
g5 = stat_res == 2;

sessionwise_stats(2).name = 'KO_Rest';
sessionwise_stats(2).GLM.mdl = fitglm(stat_res(g1),stat_data(g1), 'Mean_Corr ~ Responsive_network','Distribution','normal','CategoricalVars',[1],'VarNames',{'Responsive_network','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(2).GLM.mdl);
sessionwise_stats(2).GLM.fit = deviance_test{2,4};

sessionwise_stats(3).name = 'WT_Rest';
sessionwise_stats(3).GLM.mdl = fitglm(stat_res(g2),stat_data(g2), 'Mean_Corr ~ Responsive_network','Distribution','normal','CategoricalVars',[1],'VarNames',{'Responsive_network','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(3).GLM.mdl);
sessionwise_stats(3).GLM.fit = deviance_test{2,4};

sessionwise_stats(4).name = 'KO_Run';
sessionwise_stats(4).GLM.mdl = fitglm(stat_res(g3),stat_data(g3), 'Mean_Corr ~ Responsive_network','Distribution','normal','CategoricalVars',[1],'VarNames',{'Responsive_network','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(4).GLM.mdl);
sessionwise_stats(4).GLM.fit = deviance_test{2,4};

sessionwise_stats(5).name = 'WT_Run';
sessionwise_stats(5).GLM.mdl = fitglm(stat_res(g4),stat_data(g4), 'Mean_Corr ~ Responsive_network','Distribution','normal','CategoricalVars',[1],'VarNames',{'Responsive_network','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(5).GLM.mdl);
sessionwise_stats(5).GLM.fit = deviance_test{2,4};

sessionwise_stats(8).name = 'Res-Res';
sessionwise_stats(8).GLM.mdl = fitglm([stat_gen(g5),stat_mov(g5)],stat_data(g5), 'Mean_Corr ~ Genotype*Behavior','Distribution','normal','CategoricalVars',[1,2],'VarNames',{'Genotype','Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(8).GLM.mdl);
sessionwise_stats(8).GLM.fit = deviance_test{2,4};


%% Now analyse the networks separately split into 2 groups matching the split 1 -  WT and KO Separately


gKO = stat_gen == 0;
gWT = stat_gen == 1;

sessionwise_stats(9).name = 'KO';
sessionwise_stats(9).GLM.mdl = fitglm([stat_res(gKO),stat_mov(gKO)],stat_data(gKO), ...
    'Mean_Corr ~ Responsive_network*Behavior','Distribution','normal',...
    'CategoricalVars',[1,2],...
    'VarNames',{'Responsive_network','Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(9).GLM.mdl);
sessionwise_stats(9).GLM.fit = deviance_test{2,4};
% means 

sessionwise_stats(9).means = {mean(stat_data((stat_gen==0) & (stat_res ==2) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==0) & (stat_res ==2) & (stat_mov == 1)),'omitnan'),...
    mean(stat_data((stat_gen==0) & (stat_res ==1) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==0) & (stat_res ==1) & (stat_mov == 1)),'omitnan'),...
    mean(stat_data((stat_gen==0) & (stat_res ==0) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==0) & (stat_res ==0) & (stat_mov == 1)),'omitnan')}';

sessionwise_stats(9).stds = {std(stat_data((stat_gen==0) & (stat_res ==2) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==0) & (stat_res ==2) & (stat_mov == 1)),'omitnan'),...
    std(stat_data((stat_gen==0) & (stat_res ==1) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==0) & (stat_res ==1) & (stat_mov == 1)),'omitnan'),...
    std(stat_data((stat_gen==0) & (stat_res ==0) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==0) & (stat_res ==0) & (stat_mov == 1)),'omitnan')}';

sessionwise_stats(10).name = 'WT';
sessionwise_stats(10).GLM.mdl = fitglm([stat_res(gWT),stat_mov(gWT)],stat_data(gWT), 'Mean_Corr ~ Responsive_network*Behavior','Distribution','normal','CategoricalVars',[1,2],'VarNames',{'Responsive_network','Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(10).GLM.mdl);
sessionwise_stats(10).GLM.fit = deviance_test{2,4};

sessionwise_stats(10).means = {mean(stat_data((stat_gen==1) & (stat_res ==2) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==1) & (stat_res ==2) & (stat_mov == 1)),'omitnan'),...
    mean(stat_data((stat_gen==1) & (stat_res ==1) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==1) & (stat_res ==1) & (stat_mov == 1)),'omitnan'),...
    mean(stat_data((stat_gen==1) & (stat_res ==0) & (stat_mov == 0)),'omitnan'),...
    mean(stat_data((stat_gen==1) & (stat_res ==0) & (stat_mov == 1)),'omitnan')}';

sessionwise_stats(10).stds = {std(stat_data((stat_gen==0) & (stat_res ==2) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==1) & (stat_res ==2) & (stat_mov == 1)),'omitnan'),...
    std(stat_data((stat_gen==1) & (stat_res ==1) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==1) & (stat_res ==1) & (stat_mov == 1)),'omitnan'),...
    std(stat_data((stat_gen==1) & (stat_res ==0) & (stat_mov == 0)),'omitnan'),...
    std(stat_data((stat_gen==1) & (stat_res ==0) & (stat_mov == 1)),'omitnan')}';

% Post for WT
gKO_net2 = stat_gen == 0 & stat_res == 2;
gKO_net1 = stat_gen == 0 & stat_res == 1;
gWT_net2 = stat_gen == 1 & stat_res == 2;

sessionwise_stats(9).post(1).name = 'KOnet2';
sessionwise_stats(9).post(1).GLM.mdl = fitglm([stat_mov(gKO_net2)],stat_data(gKO_net2), ...
    'Mean_Corr ~ Behavior','Distribution','normal',...
    'CategoricalVars',[1],...
    'VarNames',{'Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(9).post(1).GLM.mdl);
sessionwise_stats(9).post(1).GLM.fit = deviance_test{2,4};

sessionwise_stats(9).post(2).name = 'KOnet1';
sessionwise_stats(9).post(2).GLM.mdl = fitglm([stat_mov(gKO_net1)],stat_data(gKO_net1), ...
    'Mean_Corr ~ Behavior','Distribution','normal',...
    'CategoricalVars',[1],...
    'VarNames',{'Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(9).post(2).GLM.mdl);
sessionwise_stats(9).post(2).GLM.fit = deviance_test{2,4};

sessionwise_stats(10).post(1).name = 'WTnet2';
sessionwise_stats(10).post(1).GLM.mdl = fitglm([stat_mov(gWT_net2)],stat_data(gWT_net2), ...
    'Mean_Corr ~ Behavior','Distribution','normal',...
    'CategoricalVars',[1],...
    'VarNames',{'Behavior','Mean_Corr'});
deviance_test = devianceTest(sessionwise_stats(10).post(1).GLM.mdl);
sessionwise_stats(10).post(1).GLM.fit = deviance_test{2,4};

sessionwise_stats(10).post(2).name = 'WTnet1';
%% Sessionwise plos of values - difference between rest and run


G = ismember(corrStatsTableSpeed{:,1},miceWT);
plot_data =  {list_sig_res_res_val_run(G)-list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(~G)-list_sig_res_res_val_rest(~G),...
    list_sig_res_non_val_run(G)-list_sig_res_non_val_rest(G),...
    list_sig_res_non_val_run(~G)-list_sig_res_non_val_rest(~G),...
    list_sig_non_non_val_run(G)-list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(~G)-list_sig_non_non_val_rest(~G),...
    };

xx_axis = [1,2,4,5,7,8];
cc_axis = [2,4,2,4,2,4];
figure
b = struct();
for cc = 1:6
    y_data = plot_data{cc}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    b(cc).bar = boxchart(x_data,y_data,...
        'MarkerStyle','none',...
        'BoxFaceColor',c_data,...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', c_data,...
        'LineWidth',1);
    hold on
    b(cc).scatter = scatter(x_data,y_data,[],'filled','MarkerFaceColor', c_data);
    hold on
end

% Stats

stat_data = [list_sig_res_res_val_run(G)-list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(~G)-list_sig_res_res_val_rest(~G),...
    list_sig_res_non_val_run(G)-list_sig_res_non_val_rest(G),...
    list_sig_res_non_val_run(~G)-list_sig_res_non_val_rest(~G),...
    list_sig_non_non_val_run(G)-list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(~G)-list_sig_non_non_val_rest(~G),...
    ]';
stat_gen = [ones(size(list_sig_res_res_val_rest(G))),...
    zeros(size(list_sig_res_res_val_rest(~G))),...
    ones(size(list_sig_res_non_val_rest(G))),...
    zeros(size(list_sig_res_non_val_rest(~G))),...
    ones(size(list_sig_non_non_val_rest(G))),...
    zeros(size(list_sig_non_non_val_rest(~G)))]';


stat_res = [2*ones(size(list_sig_res_res_val_rest(G))),...
    2*ones(size(list_sig_res_res_val_run(~G))),...
    ones(size(list_sig_res_non_val_rest(G))),...
    ones(size(list_sig_res_non_val_run(~G))),...
    zeros(size(list_sig_non_non_val_rest(G))),...
    zeros(size(list_sig_non_non_val_run(~G)))]';

% FIT GLM
sessionwise_stats(6).name = 'Difference Run - Rest';
sessionwise_stats(6).GLM.mdl = fitglm([stat_gen,stat_res],stat_data,'y ~ x1 + x2 + x1 * x2','Distribution','normal','CategoricalVars',[1,2]);
deviance_test = devianceTest(sessionwise_stats(6).GLM.mdl);
sessionwise_stats(6).GLM.fit = deviance_test{2,4};

%% Sessionwise plos of values - Difference like above but do not include one cell responsive

% plot 4 groups

G = ismember(corrStatsTableSpeed{:,1},miceWT);
plot_data =  {list_sig_res_res_val_run(G)-list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(~G)-list_sig_res_res_val_rest(~G),...
    list_sig_non_non_val_run(G)-list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(~G)-list_sig_non_non_val_rest(~G),...
    };

xx_axis = [1,2,4,5];
cc_axis = [2,4,2,4];
figure
b = struct();
for cc = 1:4
    y_data = plot_data{cc}';
    x_data = xx_axis(cc)*ones(size(y_data));
    c_data = colors(cc_axis(cc),:);
    b(cc).bar = boxchart(x_data,y_data,...
        'MarkerStyle','none',...
        'BoxFaceColor',c_data,...
        'BoxFaceAlpha',0,...
        'WhiskerLineColor', c_data,...
        'LineWidth',1);
    hold on
    b(cc).scatter = scatter(x_data,y_data,[],'filled','MarkerFaceColor', c_data);
    hold on
end

% Stats

stat_data = [list_sig_res_res_val_run(G)-list_sig_res_res_val_rest(G),...
    list_sig_res_res_val_run(~G)-list_sig_res_res_val_rest(~G),...
    list_sig_non_non_val_run(G)-list_sig_non_non_val_rest(G),...
    list_sig_non_non_val_run(~G)-list_sig_non_non_val_rest(~G),...
    ]';
stat_gen = [ones(size(list_sig_res_res_val_rest(G))),...
    zeros(size(list_sig_res_res_val_rest(~G))),...
    ones(size(list_sig_non_non_val_rest(G))),...
    zeros(size(list_sig_non_non_val_rest(~G)))]';

stat_res = [ones(size(list_sig_res_res_val_rest(G))),...
    ones(size(list_sig_res_res_val_run(~G))),...
    zeros(size(list_sig_non_non_val_rest(G))),...
    zeros(size(list_sig_non_non_val_run(~G)))]';


% FIT GLM
sessionwise_stats(7).name = 'Difference Run - Rest - Both Res';
sessionwise_stats(7).GLM.mdl = fitglm([stat_gen,stat_res],stat_data,'y ~ x1 * x2','Distribution','normal','CategoricalVars',[1,2]);
deviance_test = devianceTest(sessionwise_stats(7).GLM.mdl);
sessionwise_stats(7).GLM.fit = deviance_test{2,4};
