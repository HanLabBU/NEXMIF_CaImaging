%% Description

% This code runs through the list of all correlation values for all
% sessions and plots distribution

% Assymmetric correlation - intersection(A,B)/n(A); % custom (in utils)
% Jaccards Index - intersection(A,B)/union(A,B); % Built in

% Spet 15th 2021 - by Athif Mohamed

%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code


% load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedAsym'),'corrStatsTable')
corrSpeed = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable');
corrStatsTableSpeed = corrSpeed.corrStatsTable;
% fieldNames = corrStatsTableSpeed.Properties.VariableNames;
corrStatsTableSpeed([6,13,18:19,29,31:35,40:43,45],:) = [];


corr = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTablePearson_09_12'),'corrStatsTable');
corrStatsTable = corr.corrStatsTable;
% fieldNames = corrStatsTableSpeed.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];
%% plots

results = struct(); 

for row = 1:size(corrStatsTable,1)
    row;
    
    results(row).animal = corrStatsTable{row,1}{:};
    results(row).day = corrStatsTable{row,2}{:};
    % Correlation of cells significant in full duration
    corrMat = corrStatsTable{row,4}{:};
    threshMat_pos = corrStatsTable{row,5}{:};
    disMat_close = double(corrStatsTable{row,7}{:}<20);
    
    n = size(corrMat);
    corrMat(logical(triu(ones(n)))) = nan;
    threshMat_pos(logical(triu(ones(n)))) = nan;
    disMat_close(logical(triu(ones(n)))) = nan;
    
    m_sig = (corrMat-threshMat_pos >1e-5) & (disMat_close == 0);
    m_non_sig = (corrMat-threshMat_pos <= 1e-5) & (disMat_close == 0);
    
    results(row).c_sig = mean(corrMat(m_sig)); % Correlation of cells significant in full duration
    
    % Correlation of cells significant in resting duration
    corrMat_rest = corrStatsTableSpeed{row,4}{:};
    corrMat_rest(logical(triu(ones(n)))) = nan;
    results(row).c_sig_rest = mean(corrMat_rest(m_sig),'omitnan');
    
    % Correlation of cells significant in running duration
    corrMat_run = corrStatsTableSpeed{row,7}{:};
    corrMat_run(logical(triu(ones(n)))) = nan;
    results(row).c_sig_run = mean(corrMat_run(m_sig),'omitnan');
    
    
    %% Fraction that remained significant 
    
    
     % Correlation of cells significant in resting duration
    threshMat_rest = corrStatsTableSpeed{row,5}{:};
    threshMat_rest(logical(triu(ones(n)))) = nan;
    m_sig_rest = (corrMat_rest-threshMat_rest >1e-5) & (disMat_close == 0);
    results(row).frac_sig_rest = sum(sum(m_sig_rest.*m_sig))./sum(sum(m_sig))*100;
    
    % Correlation of cells significant in running duration
    threshMat_run = corrStatsTableSpeed{row,8}{:};
    threshMat_run(logical(triu(ones(n)))) = nan;
    m_sig_run = (corrMat_run-threshMat_run >1e-5) & (disMat_close == 0);
    results(row).frac_sig_run = sum(sum(m_sig_run.*m_sig))./sum(sum(m_sig))*100;
    
    
end



%% Sessionwise

G = ismember(corrStatsTableSpeed{:,1},miceWT);
resultsTable = struct2table(results);

WTRun = resultsTable{G,5};
KORun = resultsTable{~G,5};
WTRest = resultsTable{G,4};
KORest = resultsTable{~G,4};

figure

b2 = boxchart(2*ones(size(WTRun')),...
    WTRun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',1);
hold on
b1 = boxchart(ones(size(WTRest')),...
    WTRest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',1);
hold on

s2 = scatter(2,WTRun,[],colors(2,:),'filled');
hold on
s1 = scatter(1,WTRest,[],colors(1,:),'filled');
hold on
line([1,2],[WTRest,WTRun],'color','k')
b4 = boxchart(4*ones(size(KORun')),...
    KORun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',1);
hold on
b3 = boxchart(3*ones(size(KORest')),...
    KORest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',1);
hold on
s4 = scatter(4,KORun,[],colors(4,:),'filled');
hold on
s3 = scatter(3,KORest,[],colors(3,:),'filled');
hold on
line([3,4],[KORest,KORun],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest','Run','Rest','Run'})
ylabel('Mean correlation')

% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 inf])

[p,tbl,stats] = anova1([WTRest',WTRun',KORest',KORun'],...
    [ones(size(WTRest')),2*ones(size(WTRun')),3*ones(size(KORest')),4*ones(size(KORun'))],'off');
[hWT,pWT] = ttest(WTRest,WTRun)
[hKO,pKO] = ttest(KORest,KORun)
[hRest,pRest] = ttest2(WTRest,KORest)
[hRun,pRun] = ttest2(WTRun,KORun)
t_test_stats = [pWT;pKO;pRest;pRun]
mean_stats = [mean(WTRest);mean(WTRun);mean(KORest);mean(KORun)];
std_stats = [std(WTRest);std(WTRun);std(KORest);std(KORun)];
% [cs_a,ms_a] = multcompare(stats,'CType','bonferroni');

set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1, 'YLim',[0 inf])
box off
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.epsc'))


%% Fractions 

WTRun = resultsTable{G,7};
KORun = resultsTable{~G,7};
WTRest = resultsTable{G,6};
KORest = resultsTable{~G,6};

figure

b2 = boxchart(2*ones(size(WTRun')),...
    WTRun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',1);
hold on
b1 = boxchart(ones(size(WTRest')),...
    WTRest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',1);
hold on

s2 = scatter(2,WTRun,[],colors(2,:),'filled');
hold on
s1 = scatter(1,WTRest,[],colors(1,:),'filled');
hold on
line([1,2],[WTRest,WTRun],'color','k')
b4 = boxchart(4*ones(size(KORun')),...
    KORun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',1);
hold on
b3 = boxchart(3*ones(size(KORest')),...
    KORest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',1);
hold on
s4 = scatter(4,KORun,[],colors(4,:),'filled');
hold on
s3 = scatter(3,KORest,[],colors(3,:),'filled');
hold on
line([3,4],[KORest,KORun],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest','Run','Rest','Run'})
ylabel('Fraction (%)')

% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 inf])
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1, 'YLim',[0 inf])
box off

[p,tbl,stats] = anova1([WTRest',WTRun',KORest',KORun'],...
    [ones(size(WTRest')),2*ones(size(WTRun')),3*ones(size(KORest')),4*ones(size(KORun'))],'off');
[hWT,pWT] = ttest(WTRest,WTRun)
[hKO,pKO] = ttest(KORest,KORun)
[hRest,pRest] = ttest2(WTRest,KORest)
[hRun,pRun] = ttest2(WTRun,KORun)
t_test_stats = [pWT;pKO;pRest;pRun]
mean_stats = [mean(WTRest);mean(WTRun);mean(KORest);mean(KORun)];
std_stats = [std(WTRest);std(WTRun);std(KORest);std(KORun)];
% [cs_a,ms_a] = multcompare(stats,'CType','bonferroni');


% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Fraction of sig correlated cells rest run.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Fraction of sig correlated cells rest run.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\k_corr_speed_all\Fraction of sig correlated cells rest run.epsc'))
% 

%% Sessionwise all 

G = ismember(corrStatsTableSpeed{:,1},miceWT);
resultsTable = struct2table(results);

WTAll = resultsTable{G,3};
KOAll = resultsTable{~G,3};
WTRun = resultsTable{G,5};
KORun = resultsTable{~G,5};
WTRest = resultsTable{G,4};
KORest = resultsTable{~G,4};

figure

b3 = boxchart(3*ones(size(WTAll')),...
    WTAll',...
    'MarkerStyle','none',...
    'BoxFaceColor',0.5*(colors(2,:)+colors(1,:)),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', 0.5*(colors(2,:)+colors(1,:)),...
    'LineWidth',2);
hold on

b2 = boxchart(2*ones(size(WTRun')),...
    WTRun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
b1 = boxchart(ones(size(WTRest')),...
    WTRest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on

s2 = scatter(2,WTRun,[],colors(2,:),'filled');
hold on
s1 = scatter(1,WTRest,[],colors(1,:),'filled');
hold on
line([1,2],[WTRest,WTRun],'color','k')
hold on
s3 = scatter(3,WTAll,[],0.5*(colors(2,:)+colors(1,:)),'filled');
hold on
line([2,3],[WTRun,WTAll],'color','k')


b6 = boxchart(6*ones(size(KOAll')),...
    KOAll',...
    'MarkerStyle','none',...
    'BoxFaceColor',0.5*(colors(3,:)+colors(4,:)),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', 0.5*(colors(3,:)+colors(4,:)),...
    'LineWidth',2);
hold on

b5 = boxchart(5*ones(size(KORun')),...
    KORun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b4 = boxchart(4*ones(size(KORest')),...
    KORest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
s5 = scatter(5,KORun,[],colors(4,:),'filled');
hold on
s4 = scatter(4,KORest,[],colors(3,:),'filled');
hold on
line([4,5],[KORest,KORun],'color','k')
hold on
s6 = scatter(6,KOAll,[],0.5*(colors(3,:)+colors(4,:)),'filled');
hold on
line([5,6],[KORun,KOAll],'color','k')

xlim([0.5,6.5])
% ylim([0,27])
xticks([1,2,3,4,5,6])
xticklabels({'Rest','Run','All','Rest','Run','All'})
ylabel('Mean correlation')

