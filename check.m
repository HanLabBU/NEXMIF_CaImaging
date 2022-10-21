

% figure
% boxplot([list_run_sig_mov_val_WT';...
%     list_run_sig_non_mov_val_WT';...
%     list_run_sig_mov_val_KO';...
%     list_run_sig_non_mov_val_KO'],...
%     [ones(size(list_run_sig_mov_val_WT'));...
%     2*ones(size(list_run_sig_non_mov_val_WT'));...
%     3*ones(size(list_run_sig_mov_val_KO'));...
%     4*ones(size(list_run_sig_non_mov_val_KO'))],'Colors','k')
% hold on
% s1 = scatter(1,list_run_sig_mov_val_WT',[],'filled','MarkerFaceColor', colors(2,:));
% hold on
% s2 = scatter(2,list_run_sig_non_mov_val_WT',[],'filled','MarkerFaceColor', colors(1,:));
% hold on
% s3 = scatter(3,list_run_sig_mov_val_KO',[],'filled','MarkerFaceColor', colors(4,:));
% hold on
% s4 = scatter(4,list_run_sig_non_mov_val_KO',[],'filled','MarkerFaceColor', colors(3,:));
% 
% xlim([0.5,4.5])
% xticks([1,2,3,4])
% xticklabels({'MovMod','NonMovMod','MovMod','NonMovMod'})
% [h_WT,p_WT] = ttest(list_run_sig_mov_val_WT',list_run_sig_non_mov_val_WT');
% [h_KO,p_KO] = ttest(list_run_sig_mov_val_KO',list_run_sig_non_mov_val_KO');
% [h_Mod,p_Mod] = ttest2(list_run_sig_mov_val_WT',list_run_sig_mov_val_KO');
% [h_NonMod,p_NonMod] = ttest2(list_run_sig_non_mov_val_WT',list_run_sig_non_mov_val_KO');
% 
% yy  = ylim;
% text(1.5,1.05*yy(2),sprintf('p = %.4f',p_WT),'HorizontalAlignment','Center')
% text(3.5,1.05*yy(2),sprintf('p = %.4f',p_KO),'HorizontalAlignment','Center')
% text(2,1.3*yy(2),sprintf('p = %.4f',p_Mod),'HorizontalAlignment','Center')
% text(3,1.3*yy(2),sprintf('p = %.4f',p_NonMod),'HorizontalAlignment','Center')
% ylim([yy(1),yy(2)*1.5])
% 
% ylabel({'Correlation strength'})


figure
boxplot([list_rest_sig_mov_val_WT';...
    list_rest_sig_non_mov_val_WT';...
    list_rest_sig_mov_val_KO';...
    list_rest_sig_non_mov_val_KO'],...
    [ones(size(list_rest_sig_mov_val_WT'));...
    2*ones(size(list_rest_sig_non_mov_val_WT'));...
    3*ones(size(list_rest_sig_mov_val_KO'));...
    4*ones(size(list_rest_sig_non_mov_val_KO'))],'Colors','k')
hold on
s1 = scatter(1,list_rest_sig_mov_val_WT',[],'filled','MarkerFaceColor', colors(2,:));
hold on
s2 = scatter(2,list_rest_sig_non_mov_val_WT',[],'filled','MarkerFaceColor', colors(1,:));
hold on
s3 = scatter(3,list_rest_sig_mov_val_KO',[],'filled','MarkerFaceColor', colors(4,:));
hold on
s4 = scatter(4,list_rest_sig_non_mov_val_KO',[],'filled','MarkerFaceColor', colors(3,:));

xlim([0.5,4.5])
xticks([1,2,3,4])
xticklabels({'MovMod','NonMovMod','MovMod','NonMovMod'})
[h_WT,p_WT] = ttest(list_rest_sig_mov_val_WT',list_rest_sig_non_mov_val_WT');
[h_KO,p_KO] = ttest(list_rest_sig_mov_val_KO',list_rest_sig_non_mov_val_KO');
[h_Mod,p_Mod] = ttest2(list_rest_sig_mov_val_WT',list_rest_sig_mov_val_KO');
[h_NonMod,p_NonMod] = ttest2(list_rest_sig_non_mov_val_WT',list_rest_sig_non_mov_val_KO');

yy  = ylim;
text(1.5,1.05*yy(2),sprintf('p = %.4f',p_WT),'HorizontalAlignment','Center')
text(3.5,1.05*yy(2),sprintf('p = %.4f',p_KO),'HorizontalAlignment','Center')
text(2,1.3*yy(2),sprintf('p = %.4f',p_Mod),'HorizontalAlignment','Center')
text(3,1.3*yy(2),sprintf('p = %.4f',p_NonMod),'HorizontalAlignment','Center')
ylim([yy(1),yy(2)*1.5])

ylabel({'Correlation strength'})
title({'Correlation of signifcantly correlated', 'Movement modulated or Non modulated cell pairs','During resting'})
