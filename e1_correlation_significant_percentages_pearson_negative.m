%% Description

% Plots of cells correlated to each other marked by their physical location
% only positive correlation plotted
%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% Begin code
t_start = tic;
% load correlation stats
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];


% load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle','responsive_cells')
% responsive_cells_table = struct2table(responsive_cells);


load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\identity_stats_pearson_negative_08_22')

pLow_WT = listLow(~listisKO)./listEs(~listisKO)*100; % _2 - As a fraction of all possible edges
pHigh_WT = listHigh(~listisKO)./listEs(~listisKO)*100;
pLow_KO = listLow(listisKO)./listEs(listisKO)*100;
pHigh_KO = listHigh(listisKO)./listEs(listisKO)*100;

figure

b2 = boxchart(2*ones(size(pHigh_WT')),...
    pHigh_WT',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(2,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(2,:),...
    'LineWidth',2);
hold on
b1 = boxchart(ones(size(pLow_WT')),...
    pLow_WT',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(1,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(1,:),...
    'LineWidth',2);
hold on

s2 = scatter(2,pHigh_WT,[],colors(2,:),'filled');
hold on
s1 = scatter(1,pLow_WT,[],colors(1,:),'filled');
hold on
line([1,2],[pLow_WT',pHigh_WT'],'color','k')
b4 = boxchart(4*ones(size(pHigh_KO')),...
    pHigh_KO',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b3 = boxchart(3*ones(size(pLow_KO')),...
    pLow_KO',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
hold on
line([3,4],[pLow_KO',pHigh_KO'],'color','k')

s4 = scatter(4,pHigh_KO,[],colors(4,:),'filled');
hold on
s3 = scatter(3,pLow_KO,[],colors(3,:),'filled');

xlim([0.5,4.5])
% ylim([0,27])
xticks([1,2,3,4])
xticklabels({'Rest','Run','Rest','Run'})
ylabel('Correlated pairs (% of all possible pairs)')
% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.5 .3 4 4],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
[p,tbl,stats] = anova1([pLow_WT,pHigh_WT,pLow_KO,pHigh_KO],...
                       [ones(size(pLow_WT)),2*ones(size(pHigh_WT)),3*ones(size(pLow_KO)),4*ones(size(pHigh_KO))],'off');
[hWT,pWT] = ttest(pLow_WT,pHigh_WT)
[hKO,pKO] = ttest(pLow_KO,pHigh_KO)

if p<0.05
        [c,m] = multcompare(stats,'Display','off','cType','bonferroni');
end
                   %     if p<0.05
%         [c,m] = multcompare(stats,'Display','off');
%         title(sprintf('Perc of sig corr out of all edges - %s \n P-anova: %.3f \n P-mov-rest: %.3f, P-mov-shared: %.3f,  P-rest-shared: %.3f,' ,ttNames{tt},p, c(1,6),c(2,6),c(3,6)))
%     else
%         title(sprintf('Perc of sig corr out of all edges - %s \n P-anova: %.3f' ,ttNames{tt},p))
%     end
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges.fig')));
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges.png')));
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges.epsc')));


%% Fisher - cellwise

% cLow_WT = sum(listLow(~listisKO))./sum(listEs(~listisKO))*100; % _2 - As a fraction of all possible edges
% cHigh_WT = sum(listHigh(~listisKO))./sum(listEs(~listisKO))*100;
% cLow_KO = sum(listLow(listisKO))./sum(listEs(listisKO))*100;
% cHigh_KO = sum(listHigh(listisKO))./sum(listEs(listisKO))*100;

fisher_struct = struct();
for ff = 1:4
    switch ff
        case 1 % All low 
            fisher_struct(ff).type = 'Resting'
            fisher_struct(ff).fisherMat = zeros(2);
            fisher_struct(ff).fisherMat(1,1) = sum(listLow(~listisKO));
            fisher_struct(ff).fisherMat(1,2) = sum(listEs(~listisKO))-sum(listLow(~listisKO));
            fisher_struct(ff).fisherMat(2,1) = sum(listLow(listisKO));
            fisher_struct(ff).fisherMat(2,2) = sum(listEs(listisKO))-sum(listLow(listisKO));
        case 2 % All high 
            fisher_struct(ff).type = 'Running'
            fisher_struct(ff).fisherMat = zeros(2);
            fisher_struct(ff).fisherMat(1,1) = sum(listHigh(~listisKO));
            fisher_struct(ff).fisherMat(1,2) = sum(listEs(~listisKO))-sum(listHigh(~listisKO));
            fisher_struct(ff).fisherMat(2,1) = sum(listHigh(listisKO));
            fisher_struct(ff).fisherMat(2,2) = sum(listEs(listisKO))-sum(listHigh(listisKO));
        case 3 % All WT
            fisher_struct(ff).type = 'WT'
            fisher_struct(ff).fisherMat = zeros(2);
            fisher_struct(ff).fisherMat(1,1) = sum(listLow(~listisKO));
            fisher_struct(ff).fisherMat(1,2) = sum(listEs(~listisKO))-sum(listLow(~listisKO));
            fisher_struct(ff).fisherMat(2,1) = sum(listHigh(~listisKO));
            fisher_struct(ff).fisherMat(2,2) = sum(listEs(~listisKO))-sum(listHigh(~listisKO));
        case 4 % All KO
            fisher_struct(ff).type = 'KO'
            fisher_struct(ff).fisherMat = zeros(2);
            fisher_struct(ff).fisherMat(1,1) = sum(listLow(listisKO));
            fisher_struct(ff).fisherMat(1,2) = sum(listEs(listisKO))-sum(listLow(listisKO));
            fisher_struct(ff).fisherMat(2,1) = sum(listHigh(listisKO));
            fisher_struct(ff).fisherMat(2,2) = sum(listEs(listisKO))-sum(listHigh(listisKO));
    end
    
   [fisher_struct(ff).h_f,...
    fisher_struct(ff).p_f,...
    fisher_struct(ff).stats_f] = fishertest(fisher_struct(ff).fisherMat,'Alpha',0.05);%,'Tail','right'
 
    fisher_struct(ff).pr_1 = fisher_struct(ff).fisherMat(1,1)/(fisher_struct(ff).fisherMat(1,1)+ fisher_struct(ff).fisherMat(1,2))*100;
    fisher_struct(ff).SEP_1 = sqrt(fisher_struct(ff).pr_1*(100-fisher_struct(ff).pr_1)/(fisher_struct(ff).fisherMat(1,1)+ fisher_struct(ff).fisherMat(1,2)));
    fisher_struct(ff).pr_2 = fisher_struct(ff).fisherMat(2,1)/(fisher_struct(ff).fisherMat(2,1)+ fisher_struct(ff).fisherMat(2,2))*100;
    fisher_struct(ff).SEP_2 = sqrt(fisher_struct(ff).pr_2*(100-fisher_struct(ff).pr_2)/(fisher_struct(ff).fisherMat(2,1)+ fisher_struct(ff).fisherMat(2,2)));
end



figure
b = bar([1,2,3,4],[fisher_struct(1).pr_1,...
                   fisher_struct(2).pr_1,...
                   fisher_struct(1).pr_2,...
                   fisher_struct(2).pr_2]);
b.FaceColor = 'none';
b.EdgeColor = 'flat';
b.LineWidth = 2;
b.CData = colors;
hold on
xticks(gca,[1,2,3,4])
xticklabels(gca, {'WTRest','WTRun','KORest','KORun'})
er = errorbar([1,2,3,4],[fisher_struct(1).pr_1,...
                   fisher_struct(2).pr_1,...
                   fisher_struct(1).pr_2,...
                   fisher_struct(2).pr_2],...
              1.96*[fisher_struct(1).SEP_1,...
                   fisher_struct(2).SEP_1,...
                   fisher_struct(1).SEP_2,...
                   fisher_struct(2).SEP_2]);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 2;
xlim([0.5,4.5])
yy = ylim;


text(1.5,0.95*yy(2),sprintf('p* = %.4f',fisher_struct(3).p_f*4),'HorizontalAlignment','Center')
text(3.5,0.95*yy(2),sprintf('p* = %.4f',fisher_struct(4).p_f*4),'HorizontalAlignment','Center')
text(2,1.15*yy(2),sprintf('p* = %.4f',fisher_struct(1).p_f*4),'HorizontalAlignment','Center')
text(3,1.15*yy(2),sprintf('p* = %.4f',fisher_struct(2).p_f*4),'HorizontalAlignment','Center')
ylim([0,yy(2)*1.2])
ylabel('Significantly correlated cells (%)')
% set(gcf,'Color','none')
% set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',2)

%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges cellwise.fig')));
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges cellwise.png')));
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_analysis_ball',sprintf('Pearson Perc of sig negative corr out of all edges cellwise.epsc')));
