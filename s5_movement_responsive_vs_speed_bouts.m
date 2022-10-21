%% Initialize

addpath(genpath('U:\eng_research_handata\Athif Mohamed\nexmif_paper\Utils\'))  % Add utilities
init % Initialize data directories and genotypes
close all

%% begin code 
load('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\responsive_cells_bout_shuffle_08_22','responsive_cells')
responsive_cells = struct2table(responsive_cells);
load(fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\stats\movement_responsive_cells\move_dur_speed_perc_08_22'));

figure
hold on
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
s1 = scatter(res_cell_perc(1:20),  move_perc(1:20),'filled','MarkerFaceColor',colors(4,:),'SizeData',75);
s2 = scatter(res_cell_perc(21:end),  move_perc(21:end),'filled','MarkerFaceColor',colors(2,:),'SizeData',75);
l = lsline;
l(1).Color = colors(2,:);
l(1).LineWidth = 2;
l(2).Color = colors(4,:);
l(2).LineWidth = 2;

[r1,p1]= corrcoef(res_cell_perc(1:20),  move_perc(1:20));
[r2,p2]= corrcoef(res_cell_perc(21:end),  move_perc(21:end));

hold on % add another text for second line fit
text(mean(res_cell_perc(1:20))*1.75,mean(move_perc(1:20)),sprintf('KO: r^2 = %.2f p = %.4f',r1(2).^2,p1(2)),'Fontsize', 14)
text(mean(res_cell_perc(21:end))*1.75,mean(move_perc(21:end)),sprintf('WT: r^2 = %.2f p = %.4f',r2(2).^2,p2(2)),'Fontsize', 14)

xlabel('Movement modulated cells (% of total cells)'); ylabel('Running duration (% of total time)')
xlim([0 inf])
ylim([0 inf])
axis square
% set(gca,'Units','inches','OuterPosition',[0.5 0.5 4 4],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_bout.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_bout.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_bout.epsc'))

stats_running_dur = [r2(2).^2, p2(2);...
         r1(2).^2, p1(2)];

hold off

figure
hold on
%  set(gcf,'units','normalized','outerposition',[0 0 1 1])
s1 = scatter(res_cell_perc(1:20),  move_speed(1:20),'filled','MarkerFaceColor',colors(4,:),'SizeData',75);
s2 = scatter(res_cell_perc(21:end),  move_speed(21:end),'filled','MarkerFaceColor',colors(2,:),'SizeData',75);
l = lsline;
l(1).Color = colors(2,:);
l(1).LineWidth = 2;
l(2).Color = colors(4,:);
l(2).LineWidth = 2;
% text(res_cell_perc*1.01, move_speed*1.01, label, 'Fontsize', 10)

[r1,p1]= corrcoef(res_cell_perc(1:20),  move_speed(1:20));
[r2,p2]= corrcoef(res_cell_perc(21:end),  move_speed(21:end));

hold on % add another text for second line fit
text(mean(res_cell_perc(1:20))*1.75,mean(move_speed(1:20)),sprintf('KO: r^2 = %.2f p = %.4f',r1(2).^2,p1(2)),'Fontsize', 14)
text(mean(res_cell_perc(21:end))*1.75,mean(move_speed(21:end)),sprintf('WT: r^2 = %.2f p = %.4f',r2(2).^2,p2(2)),'Fontsize', 14)

xlabel('Movement modulated cells (% of total cells)'); ylabel('Average speed (cm/s)')
xlim([0 inf])
ylim([0 inf])
axis square
% set(gca,'Units','inches','OuterPosition',[0.5 0.5 4 4],'TickDir','out','TickLength',[0.03, 0.025],'Box','off','LineWidth',2)
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_speed.png'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_speed.fig'))
% saveas(gcf,fullfile('U:\eng_research_handata\Athif Mohamed\nexmif_paper\code_final\plots\supplementary\8.3.22\','responsive_vs_speed.epsc'))

stats_av_speed = [r2(2).^2, p2(2);...
         r1(2).^2, p1(2)];
% mdl = fitlm(res_cell_perc(1:20),  move_speed(1:20));

%% Stats 

genotype = [repmat({'KO'},[20,1]);repmat({'WT'},[12,1])];
t = table(label',genotype,res_cell_perc', move_speed',move_perc',VariableNames={'Mouse','Genotype','MM','Av_speed','Run_frac'});

mdl_1 = fitlm(t,'MM ~ Genotype * Av_speed','CategoricalVar',{'Genotype'})

mdl_2 = fitlm(t,'MM ~ Genotype * Run_frac','CategoricalVar',{'Genotype'})
