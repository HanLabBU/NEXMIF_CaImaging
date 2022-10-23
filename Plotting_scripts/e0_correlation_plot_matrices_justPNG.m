% Description

% This code runs through the list of all correlation values for all
% sessions and plots distribution

% Assymmetric correlation - intersection(A,B)/n(A); % custom (in utils)
% Jaccards Index - intersection(A,B)/union(A,B); % Built in

% Spet 15th 2021 - by Athif Mohamed
% edited 9.29.21 RAM for non-paper edits

% Initialize
%% Create matrices for python

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

% Begin code


corrTypeId = 3;

load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson'),'corrStatsTable')

fieldNames = corrStatsTable.Properties.VariableNames;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];


corrDir = {'neg','pos','both'};
for ccc = 2
    corrMat_low = 4; % Low speed correlation all (matrix)
    corrMat_high = 7;
    xlabels = corrType;
    for kk = 1:30

        miceName = corrStatsTable.animal{kk};
        day =  corrStatsTable.day{kk};

        % Low speed

        y_AllMat_low  = corrStatsTable{kk,corrMat_low}{:}; % All correlations for the animal
        y_AllMat_low(isnan(y_AllMat_low)) = 0;
        n = size(y_AllMat_low);
        yyy = rot90(fliplr(y_AllMat_low));
        y_AllMat_low(logical(triu(ones(n)))) = yyy(logical(triu(ones(n))));
        y_AllMat_low(logical(eye(n))) = 1;


        save(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\correlation_matrices_pearson\unclustered',sprintf('%s-%s-Day%i-%s-corr-matrix',miceName,mType{ismember(miceName,miceKO)+1},day,'low-speed')),'y_AllMat_low')

        % High speed

        y_AllMat_high  = corrStatsTable{kk,corrMat_high}{:}; % All correlations for the animal
        y_AllMat_high(isnan(y_AllMat_high)) = 0;

        n = size(y_AllMat_high);
        yyy = rot90(fliplr(y_AllMat_high));
        y_AllMat_high(logical(triu(ones(n)))) = yyy(logical(triu(ones(n))));
        y_AllMat_high(logical(eye(n))) = 1;

        save(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\correlation_matrices_pearson\unclustered',sprintf('%s-%s-Day%i-%s-corr-matrix',miceName,mType{ismember(miceName,miceKO)+1},day,'high-speed')),'y_AllMat_high')
    end
end

%% Plot  clustered matrices

% Sorted by high movement 
filedir_high = 'U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\correlation_matrices_pearson\clustered\high';
matfiles_high = dir(fullfile(filedir_high, '*.mat'));
nfiles_high = length(matfiles_high);

filedir_low = 'U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball\correlation_matrices_pearson\clustered\low';
matfiles_low = dir(fullfile(filedir_low, '*.mat'));


for i = 1 : nfiles_high
    s_h = load(fullfile(filedir_high, matfiles_high(i).name));
    s_l = load(fullfile(filedir_low, matfiles_low(i).name));
    
% %     clustered_H_H = s_h.clusteredMat;   % cluster order in high speed  correlation
% %     clustered_H_L = s_l.corrMat(s_h.row_labels+1,s_h.row_labels+1);
% %     
% %     figure
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% %     colormap('turbo')
% %     imagesc(clustered_H_H)
% %     axis square
% %     colorbar('Box','off','TickDirection','out','TickLength',0.03,'LineWidth',2)
% % %     set(gca,'Units','inches','InnerPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.03, 0.025],'LineWidth',2)
% % 
% %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.png',matfiles_high(i).name(1:end-4))))
% % %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.fig',matfiles_high(i).name(1:end-4))))
% % %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.epsc',matfiles_high(i).name(1:end-4))))
% % 
% % 
% %     figure
% %     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% %     colormap('turbo')
% %     imagesc(clustered_H_L)
% %     axis square
% %     colorbar('Box','off','TickDirection','out','TickLength',0.03,'LineWidth',2)
% % %     set(gca,'Units','inches','InnerPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.03, 0.025],'LineWidth',2)
% %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.png',matfiles_low(i).name(1:end-4))))
% % %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.fig',matfiles_low(i).name(1:end-4))))
% % %     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\high', sprintf('%s.epsc',matfiles_low(i).name(1:end-4))))
% %     
    clustered_L_L = s_l.clusteredMat;      % cluster order in low speed correlation
    clustered_L_H = s_h.corrMat(s_l.row_labels+1,s_l.row_labels+1);
    
     
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    colormap('turbo')
    imagesc(clustered_L_L)
    axis square
    colorbar('Box','off','TickDirection','out','TickLength',0.03,'LineWidth',2)
%     set(gca,'Units','inches','OuterPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.03, 0.025],'LineWidth',2)
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_matrices\all_spectral_low', sprintf('%s.png',matfiles_low(i).name(1:end-4))))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\low', sprintf('%s.fig',matfiles_low(i).name(1:end-4))))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\low', sprintf('%s.epsc',matfiles_low(i).name(1:end-4))))


    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    colormap('turbo')
    imagesc(clustered_L_H)
    axis square
    colorbar('Box','off','TickDirection','out','TickLength',0.03,'LineWidth',2)
%     set(gca,'Units','inches','OuterPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.03, 0.025],'LineWidth',2)
    saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_matrices\all_spectral_low', sprintf('%s.png',matfiles_high(i).name(1:end-4))))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\low', sprintf('%s.fig',matfiles_high(i).name(1:end-4))))
%     saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\correlation_matrices\spectral\low', sprintf('%s.epsc',matfiles_high(i).name(1:end-4))))
  
    close all
end
