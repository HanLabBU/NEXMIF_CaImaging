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
    for k = 1:30
        
        miceName = corrStatsTable.animal{k};
        day =  corrStatsTable.day{k};
        
        % Take significant correlations and rank by number of cells a cell
        % is signicantly correlated
        
        % get and threshold significant correlations
        
        y_close = double(corrStatsTable{k,10}{:}<20);
        n = size( y_close,1);
        y_close(logical(triu(ones(n)))) = nan;
        
        %% resting
        w_corr_low = corrStatsTable{k,4}{:};
        %     w_corr_low(isnan(w_corr_low)) = 0;
        w_corr_low(logical(triu(ones(n)))) = nan;
        w_thresh_low_pos = corrStatsTable{k,5}{:};
        w_thresh_low_neg = corrStatsTable{k,6}{:};
        adj_sig_corr_low_pos = (w_corr_low - w_thresh_low_pos > 1e-5) & (y_close == 0);
        
        count_low_pos = sum(adj_sig_corr_low_pos,1)+sum(adj_sig_corr_low_pos,2)'; % count sig cells
        [~,sort_low] = sort(count_low_pos,'descend');
        
        w_corr_low_sorted = w_corr_low;
        yyy = rot90(fliplr(w_corr_low_sorted));
        w_corr_low_sorted(logical(triu(ones(n)))) = yyy(logical(triu(ones(n))));
        w_corr_low_sorted(logical(eye(n))) = 1;
        w_corr_low_sorted = w_corr_low_sorted(sort_low,:);
        w_corr_low_sorted = w_corr_low_sorted(:,sort_low);
        
        figure
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        colormap('turbo')
        imagesc(w_corr_low_sorted)
        axis square
        caxis([-0.1 1])
        colorbar('Box','off','TickDirection','out','TickLength',0.02,'LineWidth',1)
%         set(gca,'Units','inches','OuterPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.02, 0.025],'LineWidth',1)
        set(gca,'TickDir','out','LineWidth',1, 'Box','off')
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_matrices\sorted_by_sig', sprintf('Rest %s-%i.png',miceName,day)))
        
        %% Running
        
        w_corr_high = corrStatsTable{k,7}{:};
        %     w_corr_high(isnan(w_corr_high)) = 0;
        w_corr_high(logical(triu(ones(n)))) = nan;
        w_thresh_high_pos = corrStatsTable{k,8}{:};
        w_thresh_high_neg = corrStatsTable{k,9}{:};
        adj_sig_corr_high_pos = (w_corr_high - w_thresh_high_pos > 1e-5) & (y_close == 0);
        
        count_high_pos = sum(adj_sig_corr_high_pos,1)+sum(adj_sig_corr_high_pos,2)'; % count sig cells
        [~,sort_high] = sort(count_high_pos,'descend');
        
        w_corr_high_sorted = w_corr_high;
        yyy = rot90(fliplr(w_corr_high_sorted));
        w_corr_high_sorted(logical(triu(ones(n)))) = yyy(logical(triu(ones(n))));
        w_corr_high_sorted(logical(eye(n))) = 1;
        w_corr_high_sorted = w_corr_high_sorted(sort_high,:);
        w_corr_high_sorted = w_corr_high_sorted(:,sort_high);
        
        figure
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        colormap('turbo')
        imagesc(w_corr_high_sorted)
        axis square
        caxis([-0.1 1])
        colorbar('Box','off','TickDirection','out','TickLength',0.02,'LineWidth',1)
%         set(gca,'Units','inches','OuterPosition',[0.5 0.5 2.7 2.3],'TickDir','out','TickLength',[0.02, 0.025],'LineWidth',1)
        set(gca,'TickDir','out','LineWidth',1, 'Box','off')
        saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_matrices\sorted_by_sig', sprintf('Run %s-%i.png',miceName,day)))
%         saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\e0_correlation_matrices\sorted_by_sig', sprintf('Run %s-%i.fig',miceName,day)))

        close all

    end
end


