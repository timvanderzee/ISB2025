function s = PlotReflexes_Force(R)
% Plot Reflexes

% Info on figure
lw = 1.5; % linewidth
s  = figure('Name','Simulated reflexes');
color_ext = [148 93 94]./255;
color_flex= [221 167 123]./255;

subplot(221)
plot(R.C.a_refl,'Color',color_ext,'LineWidth',lw); hold on
title('A refl'); box off;

subplot(223)
plot(R.C.Fsrs,'Color',color_ext,'LineWidth',lw); hold on
%plot(R.Fsrs_del,'Color',[0.7 0.7 0.7],'LineWidth',lw); hold on
title('Fsrs'); box off; legend({'Fsrs','Delayed'})

subplot(224)
plot(R.C.Fce_ext,'Color',color_ext,'LineWidth',lw); hold on
plot(R.Fce_del,'Color',[0.7 0.7 0.7],'LineWidth',lw); hold on
title('Fce'); box off; legend({'Fce ext','Delayed'})