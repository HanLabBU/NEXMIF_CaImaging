
%%
addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_ball','eventStats_ext'),'eventStatsAll')
eventStatsTable = struct2table(eventStatsAll);


maskFewM = false(size(eventStatsTable,1),1); % A flag for sessions with less than 10% mov
for ii = 1:numel(fewMovSessionsM)
    maskFewM(ismember([eventStatsTable.animal],fewMovSessionsM(ii))& ismember([eventStatsTable.day],fewMovSessionsD(ii)))= 1;
end

%% For each feature plots PDF and CDF accross all the neurons

c = 1;%:3 % ball

maskC = ismember([eventStatsTable.condition],'Ball'); %Ball
maskG1 = ismember([eventStatsTable.animal],miceStudy(maskWT)); %WT
maskG2 = ismember([eventStatsTable.animal],miceStudy(maskKO)); %KO

t3 = struct();
t3.WTRest = eventStatsTable{maskC&maskG1,21};
t3.KORest= eventStatsTable{maskC&maskG2,21};
t3.WTRun = eventStatsTable{maskC&maskG1,20};
t3.KORun = eventStatsTable{maskC&maskG2,20};
% t3.KORun(t3.KORun > 10 | t3.KORun == 0) = [];

figure
vp = violinplot(t3,[],'ShowMean',true);
vp(1,2).ViolinColor = colors(3,:);
vp(1,1).ViolinColor =  colors(1,:);
vp(1,2).BoxColor =  colors(3,:);
vp(1,1).BoxColor =  colors(1,:);
vp(1,2).EdgeColor =  colors(3,:);
vp(1,1).EdgeColor =  colors(1,:);
vp(1,4).ViolinColor = colors(4,:);
vp(1,3).ViolinColor =  colors(2,:);
vp(1,4).BoxColor =  colors(4,:);
vp(1,3).BoxColor =  colors(2,:);
vp(1,4).EdgeColor =  colors(4,:);
vp(1,3).EdgeColor =  colors(2,:);

vp(1,1).ViolinAlpha = 0.1;
vp(1,2).ViolinAlpha = 0.1;
title(sprintf('%s','Entropy - Cellwise - violin'),'FontSize',11,'FontName','Arial','FontWeight','Bold')
xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest','Rest','Run','Run'})
ylabel('Normalized entropy')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)

saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.png','Cellwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.fig','Cellwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.svg','Cellwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.epsc','Cellwise')))


% stats

% [h,p] = kstest(t3.WTRun);
% [h,p] = kstest(t3.KORun);
% [h,p] = kstest(t3.WTRest);
% [h,p] = kstest(t3.KORest);

% Kruskalwallis test

[p_k,tbl_k,stats_k] = kruskalwallis([t3.WTRun;...
                               t3.KORun;...
                               t3.WTRest;...
                               t3.KORest],...
                               [1*ones(size(t3.WTRun));...
                                2*ones(size(t3.KORun));...
                                3*ones(size(t3.WTRest));...
                                4*ones(size(t3.KORest))],'off');
[c_k,m_k] = multcompare(stats_k,'CType','bonferroni','display','off');

% Anova test 
[p_a,tbl_a,stats_a] = anova1([t3.WTRun;...
                               t3.KORun;...
                               t3.WTRest;...
                               t3.KORest],...
                               [1*ones(size(t3.WTRun));...
                                2*ones(size(t3.KORun));...
                                3*ones(size(t3.WTRest));...
                                4*ones(size(t3.KORest))],'off');
[c_a,m_a] = multcompare(stats_a,'CType','bonferroni','display','off');

%% sessionwise 

GWT = findgroups(eventStatsTable{maskC&maskG1,2});
GKO = findgroups(eventStatsTable{maskC&maskG2,2});

meanWTRun = splitapply(@nanmean,eventStatsTable{maskC&maskG1,20},GWT);
meanKORun = splitapply(@nanmean,eventStatsTable{maskC&maskG2,20},GKO);
meanWTRest = splitapply(@nanmean,eventStatsTable{maskC&maskG1,21},GWT);
meanKORest = splitapply(@nanmean,eventStatsTable{maskC&maskG2,21},GKO);

figure
b1 = boxchart(ones(size(meanWTRest)),...
    meanWTRest,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on
s1 = scatter(1,meanWTRest,[],colors(1,:),'filled');
hold on
b2 = boxchart(2*ones(size(meanKORest)),...
    meanKORest,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
s2 = scatter(2,meanKORest,[],colors(3,:),'filled');

hold on
b3 = boxchart(3*ones(size(meanWTRun)),...
    meanWTRun,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
s3 = scatter(3,meanWTRun,[],colors(2,:),'filled');


b4 = boxchart(4*ones(size(meanKORun)),...
    meanKORun,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
s4 = scatter(4,meanKORun,[],colors(4,:),'filled');

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest','Rest','Run','Run'})
ylabel('Mean normalized entropy of all cells')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
title(sprintf('%s','Entropy - Sessionwise'),'FontSize',11,'FontName','Arial','FontWeight','Bold')

saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.png','Sessionwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.fig','Sessionwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.svg','Sessionwise')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.epsc','Sessionwise')))


% stats
[p,tbl,stats] = anova1([meanWTRest;meanKORest;meanWTRun;meanKORun],...
                       [ones(size(meanWTRest));2*ones(size(meanKORest));3*ones(size(meanWTRun));4*ones(size(meanKORun))],'off');
% [h,p] = ttest2(meanWTRest,meanKORest);
% [h,p] = ttest2(meanWTRun,meanKORun);

if p<0.05
        [c,m] = multcompare(stats,'Display','off');   
end

%% WT vs KO 

meanWT = splitapply(@nanmean,eventStatsTable{maskC&maskG1,13},GWT);
meanKO = splitapply(@nanmean,eventStatsTable{maskC&maskG2,13},GKO);

figure
b1 = boxchart(ones(size(meanWT)),...
    meanWT,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
s1 = scatter(1,meanWT,[],colors(2,:),'filled');
hold on
b2 = boxchart(2*ones(size(meanKO)),...
    meanKO,...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
s2 = scatter(2,meanKO,[],colors(4,:),'filled');

xticks([1,2])
xticklabels({'WT','KO'})
ylabel('Mean normalized entropy of all cells')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
title(sprintf('%s','Entropy - Sessionwise - Wt vs KO'),'FontSize',11,'FontName','Arial','FontWeight','Bold')

[h,p] = ttest2(meanWT,meanKO);

saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.png','Sessionwise WT - KO')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.fig','Sessionwise WT - KO')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.svg','Sessionwise WT - KO')))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\h_entropy',sprintf('Entropy %s.epsc','Sessionwise WT - KO')))
