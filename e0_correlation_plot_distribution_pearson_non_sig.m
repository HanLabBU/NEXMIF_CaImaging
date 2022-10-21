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
    
    
    %     fig2 = figure;
    %     set(fig2,'units','normalized','outerposition',[0 0 1 1])
    %
    %     fig3 = figure;
    %     set(fig3,'units','normalized','outerposition',[0 0 1 1])
    %
    %     fig4 = figure;
    %     set(fig4,'units','normalized','outerposition',[0 0 1 1])
    
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
            featureMaxPlot = 0.3;
            featureMinPlot = -0.2;
            
            
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

                              c_non_sig = [c_non_sig;y_AllMat((y_AllMat-y_ThreshMat_pos<=1e-5)& (y_AllMat>=1e-5) & (y_close == 0) )];
                              c_sig = [c_sig;y_AllMat((y_AllMat-y_ThreshMat_pos >1e-5) & (y_close == 0))];  
                              corrStatsTable{row,sp+10} = mean(y_AllMat((y_AllMat-y_ThreshMat_pos<=1e-5)& (y_AllMat>=1e-5) & (y_close == 0)));
                              
                        case  3
                            c_non_sig = [c_non_sig;y_AllMat((y_AllMat-y_ThreshMat_neg >=-1e-5 & y_AllMat - y_ThreshMat_pos <=1e-5) & (y_close == 0))];
                            c_sig = [c_sig;y_AllMat((y_AllMat-y_ThreshMat_neg <-1e-5 | y_AllMat - y_ThreshMat_pos > 1e-5) & (y_close == 0))];
                            corrStatsTable{row,sp+10} = mean(y_AllMat((y_AllMat-y_ThreshMat_neg >=-1e-5 & y_AllMat - y_ThreshMat_pos <=1e-5) & (y_close == 0)));
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
            pdf_list(kk).non_sig_data = c_non_sig_gen;
            ss = ss+1;
        end
    end
end


%% overlayed
% fig1 = figure;
% % set(fig1,'units','normalized','outerposition',[0 0 1 1])
% 
% % WT
% 
% x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
% y = [0.1,pdf_list(1).sig+0.1,0.1];
% patch(x,y,colors(1,:),'FaceAlpha',0.5)
% hold on
% 
% x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
% y = [0.1,pdf_list(3).sig+0.1,0.1];
% patch(x,y,colors(2,:),'FaceAlpha',0.5)
% hold on
% % KO
% 
% x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
% y = [0 ,pdf_list(2).sig,0];
% patch(x,y,colors(3,:),'FaceAlpha',0.5)
% hold on
% 
% x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
% y = [0,pdf_list(4).sig,0];
% patch(x,y,colors(4,:),'FaceAlpha',0.5)
% 
% legend({'WT - rest','WT - running','KO - rest','KO - running'})
% set(gca,'ytick',[])
% xlabel('Asymmetric correlations')
%% stacked
fig2 = figure;
% set(fig2,'units','normalized','outerposition',[0 0 1 1])

% WT rest -- top, +0.3
x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0.3,pdf_list(1).non_sig+0.3,0.3];
patch(x,y,colors(1,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% WT run -- third, +0.1
x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0.1,pdf_list(3).non_sig+0.1,0.1];
patch(x,y,colors(2,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% KO rest -- second, +0.2
x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0+0.2,pdf_list(2).non_sig+0.2,0.2];
patch(x,y,colors(3,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% KO Run -- bottom, +0
x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0 ,pdf_list(4).non_sig,0];
patch(x,y,colors(4,:),'FaceAlpha',1,'EdgeColor','none')
hold on
% plot means of each distribution
xline(pdf_list(1).non_sig_mean,'color',colors(1,:))
xline(pdf_list(2).non_sig_mean,'color',colors(3,:))
xline(pdf_list(3).non_sig_mean,'color',colors(2,:))
xline(pdf_list(4).non_sig_mean,'color',colors(4,:))


xlabel(xlabels{corrTypeId})
ylabel('Normalized probability')
%title(sprintf('%s - %s normalized probability distribution - %s',corrNames{corrTypeId},'Ball',corrDir{ccc}))
xlim([featureMinPlot featureMaxPlot])

% legend({'WT rest','WT run','KO rest','KO run'})
% set(gcf,'Color','none')
set(gcf,'Units','inches','Position',[5,7,4,6])
set(gca,'Units','inches','InnerPosition',[1 1 2 4.5],'TickDir','out','TickLength',[0.015, 0.025],'Color','none','Box','off','LineWidth',1)

% Y ticks 
yt = get(gca,'ytick');
ytl =  get(gca,'yticklabel');
ytl(yt>0.1) = strsplit(num2str(round(yt(yt>0.1)-0.1,2)));
ytl(yt>0.2) = strsplit(num2str(round(yt(yt>0.2)-0.2,2)));
ytl(yt>0.3) = strsplit(num2str(round(yt(yt>0.3)-0.3,2)));
set(gca,'yticklabel',ytl)


% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig',sprintf('PDF %s - %s - speedwise - non sig.fig',corrNames{corrTypeId},'Ball')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig',sprintf('PDF %s - %s - speedwise - non sig.png',corrNames{corrTypeId},'Ball')))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig',sprintf('PDF %s - %s - speedwise - non sig.epsc',corrNames{corrTypeId},'Ball')))

% %% Stats label
% disp('percentage of significant cell pairs out of all pairs')
% labels = {};
% percentage = zeros(1,4);
% 
% for kk = 1:4
%     labels = [labels,count_pairs(kk).genotype]
%     percentage(kk) = count_pairs(kk).sig_pairs/(count_pairs(kk).non_sig_pairs +count_pairs(kk).sig_pairs)*100
% end

%% Compare distributions 
% Anova test 
[p_a,tbl_a,stats_a] = anova1([pdf_list(1).non_sig_data;...
                               pdf_list(2).non_sig_data;...
                               pdf_list(3).non_sig_data;...
                               pdf_list(4).non_sig_data],...
                               [1*ones(size(pdf_list(1).non_sig_data));...
                                2*ones(size(pdf_list(2).non_sig_data));...
                                3*ones(size(pdf_list(3).non_sig_data));...
                                4*ones(size(pdf_list(4).non_sig_data))]);
[c_a,m_a] = multcompare(stats_a,'CType','bonferroni');


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
set(gca,'Units','inches','InnerPosition',[.8 .5 3 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',1,'YLim',[0 inf])

[p,tbl,stats] = anova1([WTRest',WTRun',KORest',KORun'],...
                       [ones(size(WTRest')),2*ones(size(WTRun')),3*ones(size(KORest')),4*ones(size(KORun'))],'off');

[hWT,pWT] = ttest(WTRest,WTRun)
[hKO,pKO] = ttest(KORest,KORun)
[hRest,pRest] = ttest2(WTRest,KORest)
[hRun,pRun] = ttest2(WTRun,KORun)
t_test_stats = [pWT;pKO;pRest;pRun]
mean_stats = [mean(WTRest);mean(WTRun);mean(KORest);mean(KORun)];
std_stats = [std(WTRest);std(WTRun);std(KORest);std(KORun)];

saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig\8.3.22',sprintf('Sessionwise %s - %s - speedwise - %s.fig',corrNames{corrTypeId},'Ball', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig\8.3.22',sprintf('Sessionwise %s - %s - speedwise - %s.png',corrNames{corrTypeId},'Ball', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_analysis_ball_non_sig\8.3.22',sprintf('Sessionwise %s - %s - speedwise - %s.epsc',corrNames{corrTypeId},'Ball', corrDir{ccc})))
