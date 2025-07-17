function[] = effect_of_rate_parms(parms)

% make up some parameters
parms.f = 500;
parms.k11 = 20;
parms.k21 = 20;
parms.k12 = 2;
parms.k22 = 2;

xs = linspace(.5, 2, 10);
color = parula(length(xs));

for i = 1:length(xs)
    subplot(321)
    plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, parms.f * xs(i), parms.w), 'color', color(i,:)); hold on

    subplot(322)
    plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, parms.f, parms.w * xs(i)), 'color', color(i,:)); hold on

    subplot(323)
    plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11* xs(i), -parms.k12) + parms.g_func(parms.xi, parms.k21, parms.k22),'-', 'color', color(i,:)); hold on

    subplot(324)
    plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11, -parms.k12 * xs(i)) + parms.g_func(parms.xi, parms.k21, parms.k22),'-', 'color', color(i,:)); hold on
    
    subplot(325)
    plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11, -parms.k12) + parms.g_func(parms.xi, parms.k21* xs(i), parms.k22),'-', 'color', color(i,:)); hold on

    subplot(326)
    plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11, -parms.k12) + parms.g_func(parms.xi, parms.k21, parms.k22 * xs(i)),'-', 'color', color(i,:)); hold on

end

fmax = max(parms.f_func(parms.xi, parms.f * xs(end), parms.w));

%%
titles = {'Attachment function', 'Attachment function', 'Detachment function', 'Detachment function', 'Detachment function', 'Detachment function'};
subtitles = {'Effect of scaling f','Effect of scaling w','Effect of scaling k_{11}', 'Effect of scaling k_{12}','Effect of scaling k_{21}','Effect of scaling k_{22}'};

for i = 1:6
    subplot(3,2,i)
    title(titles{i})
    subtitle(subtitles{i})
    ylim([0 fmax])
%     axis tight
    xlim([-20 20])
    box off
    xlabel('Cross-bridge strain (nm)')
    ylabel('Rate (s^{-1})')
end
end