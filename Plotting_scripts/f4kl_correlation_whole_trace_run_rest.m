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
    results(row).c_non_sig_rest = mean(corrMat_rest(m_non_sig),'omitnan');


    % Correlation of cells significant in running duration
    corrMat_run = corrStatsTableSpeed{row,7}{:};
    corrMat_run(logical(triu(ones(n)))) = nan;
    results(row).c_sig_run = mean(corrMat_run(m_sig),'omitnan');
    results(row).c_non_sig_run = mean(corrMat_run(m_non_sig),'omitnan');

   
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



%% Sessionwise plot sig

% plot 4 groups 

G = ismember(corrStatsTableSpeed{:,1},miceWT);
resultsTable = struct2table(results);

WTrun = resultsTable.c_sig_run(G);
KOrun = resultsTable.c_sig_run(~G);
WTrest = resultsTable.c_sig_rest(G);
KOrest = resultsTable.c_sig_rest(~G);

figure
b2 = boxchart(2*ones(size(WTrun')),...
    WTrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',1);
hold on
b1 = boxchart(ones(size(WTrest')),...
    WTrest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',1);
hold on

s2 = scatter(2,WTrun,18,colors(2,:),'filled');
hold on
s1 = scatter(1,WTrest,18,colors(1,:),'filled');
hold on
% line([1,2],[WTRest,WTRun],'color','k')
b4 = boxchart(4*ones(size(KOrun')),...
    KOrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',1);
hold on
b3 = boxchart(3*ones(size(KOrest')),...
    KOrest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',1);
hold on
s4 = scatter(4,KOrun,18,colors(4,:),'filled');
hold on
s3 = scatter(3,KOrest,18,colors(3,:),'filled');
hold on
% line([3,4],[KORest,KORun],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
% xticklabels({'Rest','Run','Rest','Run'})
% ylabel('Mean correlation')

% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .75 1.25 1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 inf])
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1, 'YLim',[0 inf])
box off

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.epsc'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Correlation of sig correlated cells rest run.pdf'))



%% Stats - 4 Way comparison sig

sessionwise_stats = struct();

% Do  ANOVA sessionwise
yData = [WTrest;WTrun;KOrest;KOrun];
gType = [ones(size(WTrest));ones(size(WTrun));zeros(size(KOrest));zeros(size(KOrun))]; % WT - 1, KO - 0;
movSp = [zeros(size(WTrest));ones(size(WTrun));zeros(size(KOrest));ones(size(KOrun))]; % Running - 1, Resting - 0;

% Test for normality  - Shapiro - Wilks
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
sessionwise_stats(1).isNormal = isNormal';
sessionwise_stats(1).means = [mean(WTrest);mean(WTrun);mean(KOrest);mean(KOrun)];
sessionwise_stats(1).stds= [std(WTrest);std(WTrun);std(KOrest);std(KOrun)];

% Do two way anova
[sessionwise_stats(1).ANOVA.p...
    sessionwise_stats(1).ANOVA.tbl...
    sessionwise_stats(1).ANOVA.stats...
    sessionwise_stats(1).ANOVA.terms] = anovan(yData,{gType,movSp},"model","interaction",'display','off');


% FIT GLM
sessionwise_stats(1).GLM.mdl = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl);
sessionwise_stats(1).GLM.fit = deviance_test{2,4};

% GLM separately for WT and KO 

sessionwise_stats(1).GLM.mdl_WT = fitglm([movSp(gType==1)],yData(gType==1),' y ~ x1','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl_WT);
sessionwise_stats(1).GLM.fit_WT = deviance_test{2,4};

sessionwise_stats(1).GLM.mdl_KO = fitglm([movSp(gType==0)],yData(gType==0),' y ~ x1','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl_KO);
sessionwise_stats(1).GLM.fit_KO = deviance_test{2,4};

% GLM separately for Rest and Run 

sessionwise_stats(1).GLM.mdl_Rest = fitglm([gType(movSp==0)],yData(movSp==0),' y ~ x1','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl_Rest);
sessionwise_stats(1).GLM.fit_Rest = deviance_test{2,4};

sessionwise_stats(1).GLM.mdl_Run = fitglm([gType(movSp==1)],yData(movSp==1),' y ~ x1','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(1).GLM.mdl_Run);
sessionwise_stats(1).GLM.fit_Run = deviance_test{2,4};

% Ranksum test separately for WT and KO 
sessionwise_stats(1).GLM.ramsum_KO_p = ranksum(yData(gType==0 & movSp==0),yData(gType==0 & movSp==1));
sessionwise_stats(1).GLM.ramsum_WT_p = ranksum(yData(gType==1 & movSp==0),yData(gType==1 & movSp==1));


% Separate T tests 


[sessionwise_stats(1).Ttest.hWT,sessionwise_stats(1).Ttest.pWT] = ttest(WTrest,WTrun);
[sessionwise_stats(1).Ttest.hKO,sessionwise_stats(1).Ttest.pKO] = ttest(KOrest,KOrun);
[sessionwise_stats(1).Ttest.hRest,sessionwise_stats(1).Ttest.pRest] = ttest2(WTrest,KOrest);
[sessionwise_stats(1).Ttest.hRun,sessionwise_stats(1).Ttest.pRun] = ttest2(WTrun,KOrun);


%% Sessionwise plot non sig

% plot 4 groups 


WTrun = resultsTable.c_non_sig_run(G);
KOrun = resultsTable.c_non_sig_run(~G);
WTrest = resultsTable.c_non_sig_rest(G);
KOrest = resultsTable.c_non_sig_rest(~G);

figure
b2 = boxchart(2*ones(size(WTrun')),...
    WTrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',1);
hold on
b1 = boxchart(ones(size(WTrest')),...
    WTrest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',1);
hold on

s2 = scatter(2,WTrun,18,colors(2,:),'filled');
hold on
s1 = scatter(1,WTrest,18,colors(1,:),'filled');
hold on
% line([1,2],[WTRest,WTRun],'color','k')
b4 = boxchart(4*ones(size(KOrun')),...
    KOrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',1);
hold on
b3 = boxchart(3*ones(size(KOrest')),...
    KOrest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',1);
hold on
s4 = scatter(4,KOrun,18,colors(4,:),'filled');
hold on
s3 = scatter(3,KOrest,18,colors(3,:),'filled');
hold on
% line([3,4],[KORest,KORun],'color','k')

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
% xticklabels({'Rest','Run','Rest','Run'})
% ylabel('Mean correlation')

% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .75 1.25 1],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
set(gca,'TickDir', 'out','TickLength',[0.03, 0.025], 'Color','none','LineWidth',1)
box off


%% Do stats for non-sig 


yData = [WTrest;WTrun;KOrest;KOrun];
gType = [ones(size(WTrest));ones(size(WTrun));zeros(size(KOrest));zeros(size(KOrun))]; % WT - 1, KO - 0;
movSp = [zeros(size(WTrest));ones(size(WTrun));zeros(size(KOrest));ones(size(KOrun))]; % Running - 1, Resting - 0;


% FIT GLM
sessionwise_stats(2).GLM.mdl = fitglm([gType,movSp],yData,' y ~ x1 + x2 + x1 * x2','Distribution','normal');
deviance_test = devianceTest(sessionwise_stats(2).GLM.mdl);
sessionwise_stats(2).GLM.fit = deviance_test{2,4};

sessionwise_stats(2).means = [mean(WTrest);mean(WTrun);mean(KOrest);mean(KOrun)];
sessionwise_stats(2).stds= [std(WTrest);std(WTrun);std(KOrest);std(KOrun)];


%% % Plot differences
figure
b1 = boxchart(ones(size(WTrest'-WTrun')),...
    WTrest'-WTrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',1);
hold on
b2 = boxchart(2*ones(size(KOrest'-KOrun')),...
     KOrest'-KOrun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',1);
hold on

s1 = scatter(1,WTrest-WTrun,[],colors(2,:),'filled');
hold on
s2 = scatter(2,KOrest-KOrun,[],colors(4,:),'filled');
hold on

xlim([0.5,2.5])
% ylim([0,27])
xticks([1,2])
xticklabels({'WT','KO'})
ylabel('Mean correlation diff')

% set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
set(gca,'Units','inches','InnerPosition',[.8 .75 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)

box off

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Diff Correlation of sig correlated cells rest run.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Diff Correlation of sig correlated cells rest run.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\k_corr_speed_all\Diff Correlation of sig correlated cells rest run.epsc'))

%% Stats - diff 2  Way comparison


% Do  ANOVA sessionwise
yData = [WTrest-WTrun;KOrest-KOrun];
gType = [ones(size(WTrest));zeros(size(KOrest))]; % WT - 1, KO - 0;

% Test for normality  - Shapiro - Wilks
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
sessionwise_stats(2).isNormal = isNormal';
sessionwise_stats(2).means = [mean(WTrest-WTrun);mean(KOrest-KOrun)];
sessionwise_stats(2).stds= [std(WTrest-WTrun);std(KOrest-KOrun)];
[sessionwise_stats(2).Ttest.h,sessionwise_stats(2).Ttest.p] = ttest2((WTrest-WTrun)',(KOrest-KOrun)');

save("U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\stats\k_corr_speed_all.mat",'sessionwise_stats')






