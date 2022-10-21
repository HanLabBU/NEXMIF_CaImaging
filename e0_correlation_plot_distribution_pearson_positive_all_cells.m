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

% cd(pathData)
corrTypeId = 3;

% load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedAsym'),'corrStatsTable')
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable')

fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

%% plots
corrDir = {'neg','pos','both'};
for ccc = 2
    
    corrMat_low = 4; % Low speed correlation all (matrix)
    corrThresh_low_pos = 5; % Low speed threshold for positive corr
    corrThresh_low_neg = 6; % Low speed threshold for neg corr
    corrMat_high = 7;
    corrThresh_high_pos = 8;
    corrThresh_high_neg = 9;
    
    
    xlabels = corrNames;
    
    
    maskWTTable = ismember([corrStatsTable.animal],miceStudy(maskWT)); %WT
    maskKOTable = ismember([corrStatsTable.animal],miceStudy(maskKO)); %KO
    maskCTable = ismember([corrStatsTable.condition],'Ball'); %Ball
    masksAll = [maskWTTable,maskKOTable,maskCTable];
    
    
    count_pairs = struct();
    pdf_list = struct();
    kk = 0;
    for sp = 1:2 % low - high
        switch sp
            case 1
                corrMat = corrMat_low;
                corrThresh_neg = corrThresh_low_neg;
                corrThresh_pos = corrThresh_low_pos;
            case 2
                corrMat = corrMat_high;
                corrThresh_neg = corrThresh_high_neg;
                corrThresh_pos = corrThresh_high_pos;
        end
        
        ss = 1;
        for subG = 1:2
            switch subG
                case 1
                    nMice = sum(maskWT);
                    miceName = miceStudy(maskWT);
                case 2
                    nMice = sum(maskKO);
                    miceName = miceStudy(maskKO);
            end
            
            
            maskG = masksAll(:,subG);
            
            
            featureMax = 1;
            featureMin = -1;
            featureMaxPlot = 0.5;
            featureMinPlot = 0;
            
            
            hBins = 250;
            hEdges = linspace(featureMin,featureMax,hBins+1);
            hMat_pdf1 = zeros(nMice,hBins);
            hMat_cdf1 = zeros(nMice,hBins);
            hMat_pdf2 = zeros(nMice,hBins);
            hMat_cdf2 = zeros(nMice,hBins);
            
            c_non_sig_gen = [];
            c_sig_gen = [];
            for m = 1:nMice
                maskM = ismember([corrStatsTable.animal], miceName(m));
                maskY =  maskM & maskG;
                y_cellAll = corrStatsTable{maskY,corrMat}; % All correlations for the animal
                y_cellThresh_neg = corrStatsTable{maskY,corrThresh_neg};
                y_cellThresh_pos = corrStatsTable{maskY,corrThresh_pos};
                y_cell_close = corrStatsTable{maskY,10};
                % account for upper triangular matrix
                c_non_sig = [];
                c_sig = [];
                rows = find(maskY);
                for cll = 1:numel(y_cellAll)
                    row = rows(cll);
                    y_AllMat  = y_cellAll{cll};
%                   y_AllMat(isnan(y_AllMat)) = 0;
                    
                    n = size(y_AllMat);
                    y_AllMat(logical(triu(ones(n)))) = nan;
                    
                    % Correct the matrices for nans in lower triangular matrix
                    y_ThreshMat_neg = y_cellThresh_neg{cll};
%                     y_ThreshMat_neg(isnan(y_ThreshMat_neg)) = 0;
%                     y_ThreshMat_neg(logical(triu(ones(n)))) = 0;
                    y_ThreshMat_neg(logical(triu(ones(n)))) = nan;
                    
                    y_ThreshMat_pos = y_cellThresh_pos{cll};
%                     y_ThreshMat_pos(isnan(y_ThreshMat_pos)) = 0;
%                     y_ThreshMat_pos(logical(triu(ones(n)))) = 0;
                    y_ThreshMat_pos(logical(triu(ones(n)))) = nan;
                    
                    y_close = double(y_cell_close{cll}<20);
                    y_close(logical(triu(ones(n)))) = nan;
                       
                    switch ccc
                        case  1
%                             c_non_sig = [c_non_sig;y_AllMat(y_AllMat>=y_ThreshMat_neg & (y_close == 0) & y_AllMat <0)];
%                             c_sig = [c_sig;y_AllMat(y_AllMat<y_ThreshMat_neg & (y_close == 0 & y_AllMat <0))];
                              c_non_sig = [c_non_sig;y_AllMat((y_AllMat -y_ThreshMat_neg >= -1e-5) & (y_close == 0))];
                              c_sig = [c_sig;y_AllMat((y_AllMat-y_ThreshMat_neg <-1e-5) & (y_close == 0))];
                        case  2
%                             c_non_sig = [c_non_sig;y_AllMat(y_AllMat<=y_ThreshMat_pos & (y_close == 0) & y_AllMat >0 )];
%                             c_sig = [c_sig;y_AllMat(y_AllMat>y_ThreshMat_pos & (y_close == 0) & y_AllMat >0)];

                              c_non_sig = [c_non_sig;y_AllMat((y_AllMat-y_ThreshMat_pos<=1e-5) & (y_close == 0) )];
                              c_sig = [c_sig;y_AllMat((y_AllMat-y_ThreshMat_pos >1e-5) & (y_close == 0))];   
                              
                              corrStatsTable{row,sp+10} = mean(y_AllMat(y_close == 0 & y_AllMat>=0));
                             

                        case  3
                            c_non_sig = [c_non_sig;y_AllMat(y_AllMat-y_ThreshMat_neg >=-1e-5 & y_AllMat - y_ThreshMat_pos <=1e-5 & (y_close == 0))];
                            c_sig = [c_sig;y_AllMat((y_AllMat-y_ThreshMat_neg <-1e-5 | y_AllMat - y_ThreshMat_pos > 1e-5) & (y_close == 0))];
                            corrStatsTable{row,sp+10} = mean(y_AllMat((y_AllMat-y_ThreshMat_neg <-1e-5 | y_AllMat - y_ThreshMat_pos > 1e-5) & (y_close == 0)));
                    end
                end
                c_non_sig_gen = [c_non_sig_gen;c_non_sig];
                c_sig_gen = [c_sig_gen;c_sig];
            end
            
            kk = kk+1;
            count_pairs(kk).genotype = mType{subG};
            count_pairs(kk).speed_label = splabel{sp};
            count_pairs(kk).sig_pairs = numel(c_sig_gen);
            count_pairs(kk).non_sig_pairs = numel(c_non_sig_gen);
            
            hEdgesPlot = 0.5*(hEdges(1:end-1)+ hEdges(2:end));
            
            hMat_gen_cdf1 = [histcounts(c_non_sig_gen,hEdges,'Normalization','cdf')];%/numel(y);
            hMat_gen_pdf1 = [histcounts(c_non_sig_gen,hEdges,'Normalization','probability')];%/numel(y);
            hMat_gen_cdf2 = [histcounts(c_sig_gen,hEdges,'Normalization','cdf')];%/numel(y);
            hMat_gen_pdf2 = [histcounts(c_sig_gen,hEdges,'Normalization','probability')];%/numel(y);
            
            
            pdf_list(kk).subG = mType{subG};
            pdf_list(kk).sp = splabel{sp};
            pdf_list(kk).non_sig = hMat_gen_pdf1;
            pdf_list(kk).sig = hMat_gen_pdf2;
            pdf_list(kk).non_sig_mean = mean(c_non_sig_gen);
            pdf_list(kk).sig_mean = mean(c_sig_gen);
            pdf_list(kk).sig_data = c_sig_gen;
            ss = ss+1;
        end
    end
end


%% Sessionwise 

G = ismember(corrStatsTable{:,1},miceWT);

WTRun = corrStatsTable{G,12};
KORun = corrStatsTable{~G,12};
WTRest = corrStatsTable{G,11};
KORest = corrStatsTable{~G,11};

figure

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
b4 = boxchart(4*ones(size(KORun')),...
    KORun',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(4,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(4,:),...
    'LineWidth',2);
hold on
b3 = boxchart(3*ones(size(KORest')),...
    KORest',...
    'MarkerStyle','none',...
    'BoxFaceColor',colors(3,:),...
    'BoxFaceAlpha',0,...
    'WhiskerLineColor', colors(3,:),...
    'LineWidth',2);
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

set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)

[p,tbl,stats] = anova1([WTRest',WTRun',KORest',KORun'],...
                       [ones(size(WTRest')),2*ones(size(WTRun')),3*ones(size(KORest')),4*ones(size(KORun'))],'off');
[hWT,pWT] = ttest(WTRest,WTRun)
[hKO,pKO] = ttest(KORest,KORun)
[hRest,pRest] = ttest2(WTRest,KORest)
[hRun,pRun] = ttest2(WTRun,KORun)
t_test_stats = [pWT;pKO;pRest;pRun]


mean_stats = [mean(WTRest);mean(WTRun);mean(KORest);mean(KORun)];

% [cs_a,ms_a] = multcompare(stats,'CType','bonferroni');

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Sessionwise %s - %s - speedwise - %s - all.fig',corrNames{corrTypeId},'Ball', corrDir{ccc})))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Sessionwise %s - %s - speedwise - %s - all.png',corrNames{corrTypeId},'Ball', corrDir{ccc})))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Sessionwise %s - %s - speedwise - %s - all.epsc',corrNames{corrTypeId},'Ball', corrDir{ccc})))




%% %% Fisher test within speed, between WT and KO

count_pairs_low = [count_pairs(1:2).sig_pairs; count_pairs(1:2).non_sig_pairs].';
[h_f_low,p_f_low,stats_f_low] = fishertest(count_pairs_low,'Alpha',0.05)
[h_f_low_l,p_f_low_l,stats_f_low_l] = fishertest(count_pairs_low,'Alpha',0.05,'Tail','left')

figure
% set(gcf,'units','normalized','outerposition',[0 0 0.5 0.7])
plotFisher_cor(count_pairs_low,sprintf('Proportion of sig. corr. pairs - low speed - %s', corrDir{ccc}) ,colors(1,:),colors(3,:),mType)
ylabel('Correlated cell pairs (%)')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
%     hold on
%     text(0,mean(ylim)+0.1*diff(ylim),sprintf('Fisher sig: p = %d',p_f_low),'HorizontalAlignment','center')
%     text(0,mean(ylim)-0.1*diff(ylim),sprintf('Fisher sig left: p = %d',p_f_low_l),'HorizontalAlignment','center')

%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - Low Speed - %s.fig',corrNames{corrTypeId},'Ball', corrDir{ccc} )))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - Low Speed - %s.png',corrNames{corrTypeId},'Ball', corrDir{ccc} )))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - Low Speed - %s.epsc',corrNames{corrTypeId},'Ball', corrDir{ccc} )))

%%

count_pairs_high = [count_pairs(3:4).sig_pairs; count_pairs(3:4).non_sig_pairs].';
[h_f_high,p_f_high,stats_f_high] = fishertest(count_pairs_high,'Alpha',0.05)
[h_f_high_l,p_f_high_l,stats_f_high_l] = fishertest(count_pairs_high,'Alpha',0.05,'Tail','left')
figure
% set(gcf,'units','normalized','outerposition',[0 0 0.5 0.7])
plotFisher_cor(count_pairs_high,sprintf('Proportion of sig. corr. pairs - high speed - %s', corrDir{ccc}) ,colors(2,:),colors(4,:),mType)
ylabel('Correlated cell pairs (%)')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1)
%     hold on
%     text(0,mean(ylim)+0.1*diff(ylim),sprintf('Fisher sig: p = %d',p_f_high),'HorizontalAlignment','center')
%     text(0,mean(ylim)-0.1*diff(ylim),sprintf('Fisher sig left: p = %d',p_f_high_l),'HorizontalAlignment','center')

% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - High Speed - %s.fig',corrNames{corrTypeId},'Ball', corrDir{ccc} )))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - High Speed - %s.png',corrNames{corrTypeId},'Ball', corrDir{ccc} )))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball',sprintf('Fisher test %s - %s - High Speed - %s.epsc',corrNames{corrTypeId},'Ball', corrDir{ccc} )))
