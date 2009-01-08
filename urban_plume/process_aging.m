load -ascii aging_a_wc
load -ascii aging_f_wc
load -ascii aging_ea_wc
load -ascii aging_ef_wc
load -ascii aging_h_wc

time = aging_a_wc(:,1);
height = aging_h_wc(:,2);
activ_001 = aging_a_wc(:,2);
fresh_001  = aging_f_wc(:,2);
total_001 = activ_001 + fresh_001;

activ_01 = aging_a_wc(:,11);
fresh_01 = aging_f_wc(:,11);
total_01 = activ_01 + fresh_01;

e_activ_001 = aging_ea_wc(:,2);
e_fresh_001 = aging_ef_wc(:,2);

e_activ_01 = aging_ea_wc(:,11);
e_fresh_01 = aging_ef_wc(:,11);
 
lambda = 1.5e-5; % s^{-1}

activ_001_dot = (activ_001(2:end) - activ_001(1:(end-1))) ./ (time(2:end) - time(1:(end-1)));
activ_01_dot  = (activ_01(2:end)  - activ_01(1:(end-1)))  ./ (time(2:end) - time(1:(end-1)));

fresh_001_dot = (fresh_001(2:end) - fresh_001(1:(end-1))) ./ (time(2:end) - time(1:(end-1)));
fresh_01_dot  = (fresh_01(2:end)  - fresh_01(1:(end-1)))  ./ (time(2:end) - time(1:(end-1)));

activ_001_plot = (activ_001(1:(end-1)) + activ_001(2:end)) / 2 ;
fresh_001_plot = (fresh_001(1:(end-1)) + fresh_001(2:end)) / 2 ;

activ_01_plot = (activ_01(1:(end-1)) + activ_01(2:(end))) / 2;
fresh_01_plot = (fresh_01(1:(end-1)) + fresh_01(2:(end))) / 2;

total_001_plot = (total_001(1:(end-1)) + total_001(2:(end))) / 2;
total_01_plot = (total_01(1:(end-1)) + total_01(2:(end))) / 2;

e_activ_001_plot = e_activ_001(2:end) ./ (time(2:end) - time(1:(end-1))) ;
e_fresh_001_plot = e_fresh_001(2:end) ./ (time(2:end) - time(1:(end-1))) ;

e_activ_01_plot = e_activ_01(2:end) ./ (time(2:end) - time(1:(end-1)));
e_fresh_01_plot = e_fresh_01(2:end) ./(time(2:end) - time(1:(end-1)));

time_plot = (time(1:(end-1)) + time(2:end)) / 2;

height_dot  = (height(2:end) - height(1:(end-1))) ./ (time(2:end) - time(1:(end-1)));
height_plot = (height(1:(end-1)) + height(2:end)) / 2;

lambda_eff = lambda + max(0, height_dot ./ height_plot);

k_activ_001 = (activ_001_dot + lambda_eff .* activ_001_plot - e_activ_001_plot) ./ fresh_001_plot;
k_activ_01 = (activ_01_dot + lambda_eff .* activ_01_plot - e_activ_01_plot) ./ fresh_01_plot;  

k_fresh_001 = -(fresh_001_dot + lambda_eff .* fresh_001_plot - e_fresh_001_plot) ./ fresh_001_plot;
k_fresh_01 = -(fresh_01_dot + lambda_eff .* fresh_01_plot - e_fresh_01_plot) ./ fresh_01_plot;

total_01_dot = (total_01(2:end)  - total_01(1:(end-1)))  ./ (time(2:end) - time(1:(end-1)));
e_total_01_plot = e_activ_01_plot + e_fresh_01_plot;
total_01_rhs = - lambda_eff .* total_01_plot + e_total_01_plot;

tau_activ_001 = 1 ./ k_activ_001;
tau_fresh_001 = 1 ./ k_fresh_001;

tau_activ_01 = 1 ./ k_activ_01;
tau_fresh_01 = 1 ./ k_fresh_01;

figure
plot(time_plot/3600+6, total_01_dot, time_plot/3600+6, total_01_rhs, time_plot/3600+6, e_total_01_plot, time_plot/3600+6, lambda_eff .* total_01_plot)
legend('total_dot', 'rhs', 'e', 'dilution')
grid on
saveas(gcf, 'rhs_wc.pdf')

figure
plot(time_plot/3600+6, activ_01_dot, time_plot/3600+6, e_activ_01_plot )
legend('activ_01_dot', 'e_activ_01' )
grid on
axis([6 30 -1e6, 5e6])
saveas(gcf, 'components_wc.pdf')

figure
plot(time_plot/3600+6, activ_01_plot, time_plot/3600+6, fresh_01_plot, time_plot/3600+6, total_01_plot)
legend('activating ss = 0.01', 'not-activating ss = 0.01', 'total')
grid on
saveas(gcf,'number_wc.pdf')

figure
plot(time_plot/3600+6, k_activ_01, time_plot/3600+6, k_fresh_01)
legend('activ ss = 0.01', 'fresh ss = 0.01')
title('k')
%axis([6 30 -100 100])
grid on
saveas(gcf,'k_wc.pdf')

figure
plot(time_plot/3600+6, tau_activ_01/3600, time_plot/3600+6, tau_fresh_01/3600)
legend('activ ss = 0.01', 'fresh ss = 0.01')
title('tau in hours')
axis([6 30 -10 10])
grid on
saveas(gcf,'tau_wc.pdf')

figure
plot(time_plot/3600+6, e_activ_01_plot, time_plot/3600+6, e_fresh_01_plot)
legend('activ ss = 0.01', 'fresh ss = 0.01')
title('e')
%axis([6 30 -100 100])
grid on
saveas(gcf,'e_wc.pdf')


