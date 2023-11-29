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
% load correlation speedwise 
corrSpeed = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTableSpeedPearson_08_22'),'corrStatsTable');
corrStatsTableSpeed = corrSpeed.corrStatsTable;
corrStatsTableSpeed([6,13,18:19,29,31:35,40:43,45],:) = [];

% load correlation whole trace
corr = load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\correlation_ball','corrStatsTablePearson_09_12'),'corrStatsTable');
corrStatsTable = corr.corrStatsTable;
corrStatsTable([6,13,18:19,29,31:35,40:43,45],:) = [];

% load events stats new
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\event_features_ball','eventStats_02_23'),'eventStatsAll')
eventStatsTable = struct2table(eventStatsAll);


%% Sessionwise stats
% Get Sessionwise means
corrMat_low = 4; % Low speed correlation all (matrix)
corrThresh_low_pos = 5; % Low speed threshold for positive corr
corrThresh_low_neg = 6; % Low speed threshold for neg corr
corrMat_high = 7;
corrThresh_high_pos = 8;
corrThresh_high_neg = 9;
corrNames = {'Positive significant correlation','Positive random correlation'};

sessionwise_stats = struct;
plot_struct = table2struct(corrStatsTableSpeed(:,1:3));

for mov = 1:3 % rest-run-whole_trace
    for typ = 1:2 % sig-rand
        for rr = 1:size(corrStatsTableSpeed,1)
            eventStatsIdx = find(strcmp(eventStatsTable.animal,corrStatsTableSpeed.animal{rr}) & ...
                            eventStatsTable.day == corrStatsTableSpeed.day{rr});
            switch mov
                case 1  % rest
                    corrMat = corrStatsTableSpeed{rr,4}{:};
                    thresh = corrStatsTableSpeed{rr,5}{:};
                    cell_close = double(corrStatsTableSpeed{rr,10}{:}<20);
                    eData = eventStatsTable.EventRateLow(eventStatsIdx);
%                     movSp = [movSp;0];
                case 2  % run
                    corrMat = corrStatsTableSpeed{rr,7}{:};
                    thresh = corrStatsTableSpeed{rr,8}{:};
                    cell_close = double(corrStatsTableSpeed{rr,10}{:}<20);
                    eData = eventStatsTable.EventRateHigh(eventStatsIdx);
%                     movSp = [movSp;1];
                case 3  % whole trace
                    corrMat = corrStatsTable{rr,4}{:};
                    thresh = corrStatsTable{rr,5}{:};
                    cell_close = double(corrStatsTable{rr,7}{:}<20);
                    eData = eventStatsTable.EventRate(eventStatsIdx);

%                     movSp = [movSp;1];
            end

            n = size(corrMat);
            corrMat(logical(triu(ones(n)))) = nan;
            thresh(logical(triu(ones(n)))) = nan;
            cell_close(logical(triu(ones(n)))) = nan;
            [eData1all,eData2all] = meshgrid(eData);
            goodEIdx = ~(eData1all == 0 & eData2all == 0);

            switch typ
                case 1  % Positive significant correlation
                    corrData = corrMat((corrMat-thresh >1e-5) & (cell_close == 0) & (corrMat >1e-5) & goodEIdx);
                    eData1 = eData1all((corrMat-thresh >1e-5) & (cell_close == 0) & (corrMat >1e-5) & goodEIdx);
                    eData2 = eData2all((corrMat-thresh >1e-5) & (cell_close == 0) & (corrMat >1e-5) & goodEIdx);
                case 2  % Positive random correlations
                    corrData =corrMat((corrMat-thresh<=1e-5) & (cell_close == 0) & (corrMat>1e-5) & goodEIdx);
                    eData1 = eData1all((corrMat-thresh<=1e-5) & (cell_close == 0) & (corrMat>1e-5) & goodEIdx);
                    eData2 = eData2all((corrMat-thresh<=1e-5) & (cell_close == 0) & (corrMat>1e-5) & goodEIdx);
            end
            threshData = thresh((cell_close == 0) & (corrMat>1e-5) & goodEIdx);
            eData1pos = eData1all((cell_close == 0) & (corrMat>1e-5) & goodEIdx);
            eData2pos = eData2all((cell_close == 0) & (corrMat>1e-5) & goodEIdx);

%             if sum(threshData<0) >=1
%                 break
%             end
            switch mov
                case 1
                     switch typ
                         case 1
                             plot_struct(rr).sig_rest = [eData1,eData2,corrData];
                         case 2
                             plot_struct(rr).non_sig_rest = [eData1,eData2,corrData];
                     end
                     plot_struct(rr).thresh_rest = [eData1pos,eData2pos,threshData];
                case 2
                     switch typ
                         case 1
                             plot_struct(rr).sig_run = [eData1,eData2,corrData];
                         case 2
                             plot_struct(rr).non_sig_run = [eData1,eData2,corrData];
                     end
                     plot_struct(rr).thresh_run = [eData1pos,eData2pos,threshData];
                case 3
                     switch typ
                         case 1
                             plot_struct(rr).sig_whole_trace = [eData1,eData2,corrData];
                         case 2
                             plot_struct(rr).non_sig_whole_trace = [eData1,eData2,corrData];
                     end
                     plot_struct(rr).thresh_whole = [eData1pos,eData2pos,threshData];
            end      
%             gType = [gType;ismember(corrStatsTableSpeed{rr,1},miceWT)];
        end
    end
end


%% Plot correlation vs mean event rate for all combined 

plot_table = struct2table(plot_struct);
% 
% sig_rest = cell2mat(plot_table.sig_rest);
% non_sig_rest = cell2mat(plot_table.non_sig_rest);
% thresh_rest = cell2mat(plot_table.thresh_rest);
% 
% sig_run = cell2mat(plot_table.sig_run);
% non_sig_run = cell2mat(plot_table.non_sig_run);
% thresh_run = cell2mat(plot_table.thresh_run);

sig_whole = cell2mat(plot_table.sig_whole_trace);
non_sig_whole = cell2mat(plot_table.non_sig_whole_trace);
thresh_whole = cell2mat(plot_table.thresh_whole);

N = 20;
% figure('units','inches','position',[1,1,10,10]); plotEvCorr(N,sig_whole,'signifcant correlations') % Scatter plot for significant correlation vs event rate
% 
% figure('units','inches','position',[1,1,10,4]); plotEvCorr_Av(sig_whole,non_sig_whole) % Scatter plot for Non significant correlation vs event rate
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig.epsc'))
% 
% figure('units','inches','position',[1,1,4,4]);  plotEvThresh_Av(thresh_whole) % Scatter plot for significance threshold  vs event rate
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh.epsc'))


% Descretize event rates and plot violins (swarm plots)
% figure('units','inches','position',[1,1,9,4]);  
[fsig,fnonsig] = plotEvCorr_AvDisc(sig_whole,non_sig_whole,0); % Swarm plot for Non significant correlation vs event rate
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig_swarm_mean.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig_swarm_mean.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig_non_sig_swarm_median.epsc'))
saveas(fsig,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_sig.pdf'))
saveas(fnonsig,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_nonsig.pdf'))

% figure('units','inches','position',[1,1,2.7,2.7]);  
plotEvThresh_AvDisc(thresh_whole,0)% Swarm plot for significance threshold  vs event rate
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh_swarm_mean.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh_swarm_mean.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh_swarm_median.epsc'))
saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\Reviews_02_23\plots\s5\ev_rate_thresh_swarm_mean.pdf'))


% Plot average threshold for each discretized average event rate 
% [D,E] = discretize(mean(thresh_whole(:,1:2),2),20);
% E = E(2:end);
% [G,ID] = findgroups(D);
% M = splitapply(@mean,thresh_whole(:,3),G);
% figure
% b = plot(E(ID(2:end)-1),M(2:end));
% xlabel('Average event rate');
% ylabel('Mean threshold')
%%  Test for 1 session

% figure; scatter(mean([eData1,eData2],2),corrData,[],'b','filled')
% 
% 
% % Plot average threshold for each discretized average event rate 
% [D,E] = discretize(mean([eData1pos,eData2pos],2),20);
% E = E(1:end-1)+mean(diff(E));
% [G,ID] = findgroups(D);
% M = splitapply(@mean,threshData,G);
% figure
% b = plot(E(ID),M);
% xlabel('Average event rate');
% ylabel('Mean threshold')


% figure; 
% subplot(2,2,1); scatter3(sig_rest(:,1),sig_rest(:,2),sig_rest(:,3),10,'*');
% 
% [e,x,y] = histcounts2(sig_rest(:,1),sig_rest(:,2),N);
% dx = (x(1,2)-x(1,1))/2;
% dy = (y(1,2)-y(1,1))/2;
% subplot(2,2,2); imagesc(y(1:end-1)+dy,x(1:end-1)+dx,e,'AlphaData',1); axis xy
% [e,x,y] = histcounts2(sig_rest(:,1),sig_rest(:,3),N);
% dx = (x(1,2)-x(1,1))/2;
% dy = (y(1,2)-y(1,1))/2;
% subplot(2,2,3); imagesc(x(1:end-1)+dx,y(1:end-1)+dy,e','AlphaData',1); axis xy 
% set(gca, 'XDir','reverse')
% [e,x,y] = histcounts2(sig_rest(:,2),sig_rest(:,3),N);
% dx = (x(1,2)-x(1,1))/2;
% dy = (y(1,2)-y(1,1))/2;
% subplot(2,2,4); imagesc(x(1:end-1)+dx,y(1:end-1)+dy,e','AlphaData',1); axis xy
% 
% %% Plotting 3 D 
% % Method 1 
% N = 50;
% figure; scatter3(corrData,eData1,eData2,[],'k','filled')
% [e,x,y] = histcounts2(eData1,corrData,N);
% dx = (x(1,2)-x(1,1))/2;
% dy = (y(1,2)-y(1,1))/2;
% hold on; imagesc(y(1:end-1)+dy,x(1:end-1)+dx,e,'AlphaData',0.5)
% 
% % Method 2
% N = 50;
% figure; scatter3(eData1,eData2,corrData,[],'k','filled')
% [e,x,y] = histcounts2(eData1,corrData,N);
% [X,Y,Z] = meshgrid([-0.1,0,0.1],y(1:end-1)+dy,linspace(min(corrData),max(corrData),N));
% V = reshape(cat(3,zeros(size(e)),e,zeros(size(e))),[N,3,N]);
% hold on
% s = slice(X,Y,Z,V,[],0,[],'nearest');
% set(s,'FaceAlpha',0.5,'EdgeColor','none')

%%

