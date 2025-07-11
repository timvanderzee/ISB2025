close all
figure(1)
color = get(gca,'colororder');

xs = linspace(.5, 2, 10);
color = parula(length(xs));


for i = 1:length(xs)
subplot(221)
plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, parms.f * xs(i), parms.w), 'color', color(i,:)); hold on

subplot(222)
plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, parms.f, parms.w * xs(i)), 'color', color(i,:)); hold on

subplot(223)
plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11* xs(i), -parms.k12) + parms.g_func(parms.xi, parms.k21* xs(i), parms.k22),'-', 'color', color(i,:)); hold on

subplot(224)
plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11, -parms.k12 * xs(i)) + parms.g_func(parms.xi, parms.k21, parms.k22 * xs(i)),'-', 'color', color(i,:)); hold on

end

%%
titles = {'Attachment', 'Attachment', 'Detachment', 'Detachment'};
subtitles = {'Effect of f','Effect of w','Effect of k_{11} and k_{21}', 'Effect of k_{12} and k_{22}'};

figure(1)
for i = 1:4
    subplot(2,2,i)
    title(titles{i})
    subtitle(subtitles{i})
    ylim([0 800])
    xlim([-20 20])
    box off
    xlabel('Cross-bridge strain (nm)')
    ylabel('Rate (s^{-1})')
end