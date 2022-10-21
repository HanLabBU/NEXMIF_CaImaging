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
corrTypeId = 2;

load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_platform','corrStatsTableAsym'),'corrStatsTable')
fieldNames = corrStatsTable.Properties.VariableNames;

%% plots
corrDir = {'neg','pos','both'};
for ccc = 2
    
    corrMat = 4; % Correlation all (matrix)
    corrThresh_pos = 5; %  threshold for positive corr
    corrThresh_neg = 6; %  threshold for neg corr
    
    
    
    xlabels = corrType;
    
    
    maskWTTable = ismember([corrStatsTable.animal],miceStudy(maskWT)); %WT
    maskKOTable = ismember([corrStatsTable.animal],miceStudy(maskKO)); %KO
    maskCTable = ismember([corrStatsTable.condition],'Ball'); %Ball
    masksAll = [maskWTTable,maskKOTable,maskCTable];
    
    
    count_pairs = struct();
    pdf_list = struct();
    kk = 0;
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
        featureMin = 0;
        featureMaxPlot = 0.5;
        featureMinPlot = 0;
        
        
        hBins = 100;
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
            y_cell_close = corrStatsTable{maskY,7};
            % account for upper triangular matrix
            c_non_sig = [];
            c_sig = [];
            for cll = 1:numel(y_cellAll)
                y_AllMat  = y_cellAll{cll};
                y_AllMat(isnan(y_AllMat)) = 0;
                
                n = size(y_AllMat);
                y_AllMat(logical(triu(ones(n)))) = nan;
                
                % Correct the matrices for nans in lower triangular matrix
                y_ThreshMat_neg = y_cellThresh_neg{cll};
                y_ThreshMat_neg(isnan(y_ThreshMat_neg)) = 0;
                y_ThreshMat_neg(logical(triu(ones(n)))) = 0;
                
                y_ThreshMat_pos = y_cellThresh_pos{cll};
                y_ThreshMat_pos(isnan(y_ThreshMat_pos)) = 0;
                y_ThreshMat_pos(logical(triu(ones(n)))) = 0;
                
                y_close = double(y_cell_close{cll}<20);
                y_close(logical(triu(ones(n)))) = nan;
                
                switch ccc
                    case  1
                        c_non_sig = [c_non_sig;y_AllMat(y_AllMat>=y_ThreshMat_neg & (y_close == 0))];
                        c_sig = [c_sig;y_AllMat(y_AllMat<y_ThreshMat_neg & (y_close == 0))];
                    case  2
                        c_non_sig = [c_non_sig;y_AllMat(y_AllMat<=y_ThreshMat_pos & (y_close == 0))];
                        c_sig = [c_sig;y_AllMat(y_AllMat>y_ThreshMat_pos & (y_close == 0))];
                    case  3
                        c_non_sig = [c_non_sig;y_AllMat(y_AllMat>=y_ThreshMat_neg & y_AllMat<= y_ThreshMat_pos & (y_close == 0))];
                        c_sig = [c_sig;y_AllMat(y_AllMat<y_ThreshMat_neg | y_AllMat > y_ThreshMat_pos & (y_close == 0))];
                end
            end
            c_non_sig_gen = [c_non_sig_gen;c_non_sig];
            c_sig_gen = [c_sig_gen;c_sig];
        end
        
        kk = kk+1;
        count_pairs(kk).genotype = mType{subG};
        count_pairs(kk).sig_pairs = numel(c_sig_gen);
        count_pairs(kk).non_sig_pairs = numel(c_non_sig_gen);
        
        hEdgesPlot = 0.5*(hEdges(1:end-1)+ hEdges(2:end));
        
        hMat_gen_cdf1 = [histcounts(c_non_sig_gen,hEdges,'Normalization','cdf')];%/numel(y);
        hMat_gen_pdf1 = [histcounts(c_non_sig_gen,hEdges,'Normalization','probability')];%/numel(y);
        hMat_gen_cdf2 = [histcounts(c_sig_gen,hEdges,'Normalization','cdf')];%/numel(y);
        hMat_gen_pdf2 = [histcounts(c_sig_gen,hEdges,'Normalization','probability')];%/numel(y);
        
        
        pdf_list(kk).subG = mType{subG};
        pdf_list(kk).non_sig = hMat_gen_pdf1;
        pdf_list(kk).sig = hMat_gen_pdf2;
        pdf_list(kk).non_sig_mean = mean(c_non_sig_gen);
        pdf_list(kk).sig_mean = mean(c_sig_gen);
        
        ss = ss+1;
    end
end


%% stacked -- sig
fig2 = figure;
% set(fig2,'units','normalized','outerposition',[0 0 1 1])

% WT

x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0.15,pdf_list(1).sig+0.15,0.15];
patch(x,y,colors(2,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% KO

x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0 ,pdf_list(2).sig,0];
patch(x,y,colors(4,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% plot means of each distribution
xline(pdf_list(1).sig_mean,'color',colors(2,:))
xline(pdf_list(2).sig_mean,'color',colors(4,:))

xlabel([xlabels{corrTypeId}, ' correlation'])
ylabel('Normalized probability')
% title(sprintf('%s - %s normalized probability distribution - %s',corrNames{corrTypeId},'Platform',corrDir{ccc}))
xlim([featureMinPlot featureMaxPlot])

yt = get(gca,'ytick');
ytl =  get(gca,'yticklabel');
ytl(yt>=0.15) = strsplit(num2str(round(yt(yt>=0.15)-0.15,2)));
set(gca,'yticklabel',ytl)

set(gcf,'Units','inches','Position',[2,1,4,6])
set(gca,'Units','inches','InnerPosition',[1 1 2 2.25],'TickDir','out','TickLength',[0.015, 0.025],'Color','none','Box','off','LineWidth',2)

% legend({'WT','KO'})
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - %s.fig',corrNames{corrTypeId},'Platform', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - %s.png',corrNames{corrTypeId},'Platform', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - %s.epsc',corrNames{corrTypeId},'Platform', corrDir{ccc})))

%% stacked -- non sig
fig3 = figure;
% set(fig2,'units','normalized','outerposition',[0 0 1 1])

% WT

x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0.15,pdf_list(1).non_sig+0.15,0.15];
patch(x,y,colors(2,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% KO

x = [hEdgesPlot(1),hEdgesPlot,hEdgesPlot(end)];
y = [0 ,pdf_list(2).non_sig,0];
patch(x,y,colors(4,:),'FaceAlpha',1,'EdgeColor','none')
hold on

% plot means of each distribution
xline(pdf_list(1).non_sig_mean,'color',colors(2,:))
xline(pdf_list(2).non_sig_mean,'color',colors(4,:))

xlabel([xlabels{corrTypeId}, ' correlation'])
ylabel('Normalized probability')
% title(sprintf('%s - %s normalized probability distribution - %s',corrNames{corrTypeId},'Platform',corrDir{ccc}))
xlim([featureMinPlot featureMaxPlot])

yt = get(gca,'ytick');
ytl =  get(gca,'yticklabel');
ytl(yt>=0.15) = strsplit(num2str(round(yt(yt>=0.15)-0.15,2)));
set(gca,'yticklabel',ytl)

set(gcf,'Units','inches','Position',[2,1,4,6])
set(gca,'Units','inches','InnerPosition',[1 1 2 2.25],'TickDir','out','TickLength',[0.015, 0.025],'Color','none','Box','off','LineWidth',2)

% legend({'WT','KO'})
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - non_sig - %s.fig',corrNames{corrTypeId},'Platform', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - non_sig - %s.png',corrNames{corrTypeId},'Platform', corrDir{ccc})))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('PDF %s - %s - non_sig - %s.epsc',corrNames{corrTypeId},'Platform', corrDir{ccc})))


%% Fisher test within speed, between WT and KO

count_pairs_plot = [count_pairs(1:2).sig_pairs; count_pairs(1:2).non_sig_pairs].';
[h_f_low,p_f_low,stats_f_low] = fishertest(count_pairs_plot,'Alpha',0.05)
[h_f_low_l,p_f_low_l,stats_f_low_l] = fishertest(count_pairs_plot,'Alpha',0.05,'Tail','left')

figure
% set(gcf,'units','normalized','outerposition',[0 0 0.5 0.7])
plotFisher_cor(count_pairs_plot,sprintf('Proportion of sig. corr. pairs - %s', corrDir{ccc}) ,colors(2,:),colors(4,:),mType)
ylabel('Correlated cell pairs (%)')
set(gcf,'Color','none')
set(gca,'Units','inches','InnerPosition',[.5 .3 2 2],'TickDir','out','TickLength',[0.03, 0.025],'Color','none','Box','off','LineWidth',2)
%     hold on
%     text(0,mean(ylim)+0.1*diff(ylim),sprintf('Fisher sig: p = %d',p_f_low),'HorizontalAlignment','center')
%     text(0,mean(ylim)-0.1*diff(ylim),sprintf('Fisher sig left: p = %d',p_f_low_l),'HorizontalAlignment','center')
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('Fisher test %s - %s - %s.fig',corrNames{corrTypeId},'PLatform', corrDir{ccc} )))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('Fisher test %s - %s - %s.png',corrNames{corrTypeId},'Platform', corrDir{ccc} )))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\s3',sprintf('Fisher test %s - %s - %s.png',corrNames{corrTypeId},'Platform', corrDir{ccc} )))

