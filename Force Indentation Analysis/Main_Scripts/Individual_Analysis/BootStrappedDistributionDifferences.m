clear all

% Experimental Young's Moduli for each droplet. 
EXP_noRNA = [429.5446622
1185.726681
1248.769751
986.1089876];

EXP_FCRNA = [673.7842532
1917.686488
2150.513619
289.3375154
622.3275504
245.7371907
556.7518901
153.6452528
2127.689019
236.2173016
243.5051721];

EXP_F3RNA = [692.0447202
452.3665188
1276.538839
286.2404141
538.0961728
454.7643879
515.8766377
610.464363
804.2182637
1181.118309
522.9076897
302.5738321
1289.287856
1050.089198
981.37091
872.8331812
1007.425995
2243.698278
713.9423803];

EXP_vgRNA = [1331.109636
2444.278377
1022.51857
1528.202196
1009.29014
849.1857818
957.4421248
1279.890055
956.3697345
875.8567701
1096.461098
2326.664484
1429.426296
716.8574367
5609.046336
2997.3021
777.8419659
1074.665126
2187.627767
3043.075416
820.8677293
840.8750229
2434.603849
2254.253504
1137.329741
1350.161289
2387.898216
1248.080912
4297.324909
1570.542681
466.0239541];

%Rounding experimental values to the close power of 10
EXP_noRNA = round(EXP_noRNA,-1);
EXP_FCRNA = round(EXP_FCRNA,-1);
EXP_F3RNA = round(EXP_F3RNA,-1);
EXP_vgRNA = round(EXP_vgRNA,-1);

% Contrusctions of bootstrapped distributions.
nBootstrap = 10000;
[bootstrapMean_noRNA,ci_noRNA,bootstrapMeans_noRNA] = boostrap(EXP_noRNA,nBootstrap);
[bootstrapMean_FCRNA,ci_FCRNA,bootstrapMeans_FCRNA] = boostrap(EXP_FCRNA,nBootstrap);
[bootstrapMean_F3RNA,ci_F3RNA,bootstrapMeans_F3RNA] = boostrap(EXP_F3RNA,nBootstrap);
[bootstrapMean_vgRNA,ci_vgRNA,bootstrapMeans_vgRNA] = boostrap(EXP_vgRNA,nBootstrap);

% Compared to noRNA
[difference_FCRNA] = DifferenceOfDistributions(bootstrapMeans_FCRNA,bootstrapMeans_noRNA);
[difference_F3RNA] = DifferenceOfDistributions(bootstrapMeans_F3RNA,bootstrapMeans_noRNA);
[difference_vgRNA] = DifferenceOfDistributions(bootstrapMeans_vgRNA,bootstrapMeans_noRNA);
% Compared to FCRNA
[difference_F3RNA_FCRNA] = DifferenceOfDistributions(bootstrapMeans_F3RNA,bootstrapMeans_FCRNA);
[difference_vgRNA_FCRNA] = DifferenceOfDistributions(bootstrapMeans_vgRNA,bootstrapMeans_FCRNA);

xmin = 100;
xmax = 2500;

figure(1);
subplot(4,1,1);
[counts1, edges1] = histcounts(bootstrapMeans_noRNA, 'BinMethod', 'auto');
centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
plot(centers1, counts1/nBootstrap, '-', "Color",[0 0 0],"LineWidth",1.5);
hold on 
subplot(4,1,2);
[counts1, edges1] = histcounts(bootstrapMeans_FCRNA, 'BinMethod', 'auto');
centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
plot(centers1, counts1/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
hold on 
subplot(4,1,3);
[counts1, edges1] = histcounts(bootstrapMeans_F3RNA, 'BinMethod', 'auto');
centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
plot(centers1, counts1/nBootstrap, '-', "Color",[0.5 0.5 0.5],"LineWidth",1.5);
hold on 
subplot(4,1,4);
[counts1, edges1] = histcounts(bootstrapMeans_vgRNA, 'BinMethod', 'auto');
centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
hold on 

% Align the x-axis of all subplots
ax = findobj(gcf, 'Type', 'axes'); % Get handles to all axes
linkaxes(ax, 'x'); % Link the x-axis limits of all axes
set(ax, 'YAxisLocation', 'left'); % Move y-axis to left
xlim([xmin, xmax]); % Set the x-axis limits for all subplots







% figure(1)
% subplot(4,1,2);
% [counts1, edges1] = histcounts(difference_FCRNA, 'BinMethod', 'auto');
% centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
% plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
% hold on 
% [counts2, edges2] = histcounts(bootstrapMeans_noRNA, 'BinMethod', 'auto');
% centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% plot(centers2, counts2/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% hold on
% [counts3, edges3] = histcounts(bootstrapMeans_FCRNA, 'BinMethod', 'auto');
% centers3 = (edges3(1:end-1) + edges3(2:end)) / 2;
% plot(centers3, counts3/nBootstrap, '-', "Color",[0.75 0.75 0.75],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% xline(0,"--")
% hold on
% 
% subplot(4,1,3);
% [counts1, edges1] = histcounts(difference_F3RNA, 'BinMethod', 'auto');
% centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
% plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
% hold on 
% [counts2, edges2] = histcounts(bootstrapMeans_noRNA, 'BinMethod', 'auto');
% centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% plot(centers2, counts2/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% hold on
% [counts3, edges3] = histcounts(bootstrapMeans_F3RNA, 'BinMethod', 'auto');
% centers3 = (edges3(1:end-1) + edges3(2:end)) / 2;
% plot(centers3, counts3/nBootstrap, '-', "Color",[0.75 0.75 0.75],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% xline(0,"--")
% hold on
% 
% subplot(4,1,4);
% [counts1, edges1] = histcounts(difference_vgRNA, 'BinMethod', 'auto');
% centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
% plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
% hold on 
% [counts2, edges2] = histcounts(bootstrapMeans_noRNA, 'BinMethod', 'auto');
% centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% plot(centers2, counts2/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% hold on
% [counts3, edges3] = histcounts(bootstrapMeans_vgRNA, 'BinMethod', 'auto');
% centers3 = (edges3(1:end-1) + edges3(2:end)) / 2;
% plot(centers3, counts3/nBootstrap, '-', "Color",[0.75 0.75 0.75],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% xline(0,"--")
% hold off
% 
% 
% subplot(4,1,1);
% [counts1, edges1] = histcounts(difference_F3RNA_FCRNA, 'BinMethod', 'auto');
% centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
% plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
% hold on 
% [counts2, edges2] = histcounts(bootstrapMeans_FCRNA, 'BinMethod', 'auto');
% centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% plot(centers2, counts2/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% hold on
% [counts3, edges3] = histcounts(bootstrapMeans_F3RNA, 'BinMethod', 'auto');
% centers3 = (edges3(1:end-1) + edges3(2:end)) / 2;
% plot(centers3, counts3/nBootstrap, '-', "Color",[0.75 0.75 0.75],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% xline(0,"--")
% hold on
% 
% subplot(5,1,1);
% [counts1, edges1] = histcounts(difference_vgRNA_FCRNA, 'BinMethod', 'auto');
% centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
% plot(centers1, counts1/nBootstrap, '-', "Color",[1 0 0],"LineWidth",1.5);
% hold on 
% [counts2, edges2] = histcounts(bootstrapMeans_FCRNA, 'BinMethod', 'auto');
% centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
% plot(centers2, counts2/nBootstrap, '-', "Color",[0.25 0.25 0.25],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% hold on
% [counts3, edges3] = histcounts(bootstrapMeans_vgRNA, 'BinMethod', 'auto');
% centers3 = (edges3(1:end-1) + edges3(2:end)) / 2;
% plot(centers3, counts3/nBootstrap, '-', "Color",[0.75 0.75 0.75],"LineWidth",1.5);
% hold on
% xlim([xmin xmax]);
% xline(0,"--")
% hold on

% Set the same y-axis scale for all subplots
allAxes = findobj(gcf, 'type', 'axes');
set(allAxes, 'YLim', [0, 0.15]);

% % Shift y-axis ticks to the right for all subplots
% set(allAxes, 'YAxisLocation', 'right');

P_count_noRNA_FCRNA = sum(difference_FCRNA <= 0);
P_count_noRNA_F3RNA = sum(difference_F3RNA <= 0);
P_count_noRNA_vgRNA = sum(difference_vgRNA <= 0);
P_count_FCRNA_F3RNA = sum(difference_F3RNA_FCRNA <= 0);
P_count_FCRNA_vgRNA = sum(difference_vgRNA_FCRNA <= 0);