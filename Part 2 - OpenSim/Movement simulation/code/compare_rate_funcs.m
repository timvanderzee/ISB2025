function[] = compare_rate_funcs(newparms, parms)

color = get(gca,'colororder');

plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, parms.f, parms.w), 'color', color(1,:)); hold on
plot(parms.xi*parms.h*1e9, parms.f_func(parms.xi, newparms.f, newparms.w), 'color', color(2,:));
plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, parms.k11, -parms.k12) + parms.g_func(parms.xi, parms.k21, parms.k22),'--', 'color', color(1,:)); hold on
plot(parms.xi*parms.h*1e9, parms.g_func(parms.xi, newparms.k11, -newparms.k12) + parms.g_func(parms.xi, newparms.k21, newparms.k22),'--', 'color', color(2,:)); hold on
ylim([0 5000])
box off
legend('Attachment (new)', 'Attachment (old)', 'Detachment (new)', 'Detachment (old)', 'location', 'best')
legend boxoff
xlabel('Cross-bridge strain (nm)')
ylabel('Rate (s^{-1})')

end