load -ascii aging_a_nc
load -ascii aging_f_nc
load -ascii aging_ea_nc
load -ascii aging_ef_nc

time = aging_a_nc(:,1);
activ_001 = aging_a_nc(:,2);
fresh_001  = aging_f_nc(:,2);
total_001 = activ_001 + fresh_001;

activ_01 = aging_a_nc(:,11);
fresh_01 = aging_f_nc(:,11);
total_01 = activ_01 + fresh_01;

e_activ_01 = aging_ea_nc(:,2);
e_fresh_01 = aging_ef_nc(:,2);

e_activ_001 = aging_ea_nc(:,11);
e_activ_001 = aging_ef_nc(:,2);
 
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

e_activ_001_plot = e_activ_001(2:end);
e_fresh_001_plot = e_fresh_001(2:end);

e_activ_01_plot = e_activ_01(2:end);
e_fresh_01_plot = e_fresh_01(2:end);

time_plot = (time(1:(end-1)) + time(2:end)) / 2;

k_activ_001 = (activ_001_dot + lambda * activ_001_plot - e_activ_001_plot) ./ fresh_001_plot;
k_activ_01 = (activ_01_dot + lambda * activ_01_plot - e_activ_01_plot) ./ fresh_01_plot;  

k_fresh_001 = -(fresh_001_dot + lambda * fresh_001_plot - e_fresh_001_plot) ./ fresh_001_plot;
k_activ_01 = -(fresh_01_dot + lambda * fresh_01_plot - e_fresh_01_plot) ./ fresh_01_plot;

figure
plot(time_plot/3600+6, activ_001_plot, time_plot/3600+6, fresh_001_plot, time_plot/3600+6, total_001_plot, 
time_plot/3600+6, activ_01_plot, time_plot/3600+6, fresh_plot_01, time_plot/3600+6, total_01_plot)
legend('activating ss = 0.001', 'not-activating ss = 0.001', 'total', 'activating ss = 0.01', 'not-activating ss = 0.01', 'total')
grid on
saveas(gcf,'number_nc.pdf')

figure
plot(time_plot/3600+6, k_activ_001/3600, time_plot/3600+6, k_activ_01/3600, time_plot/3600+6, k_fresh_001/3600, time_plot/3600+6, k_fresh_01/3600)
legend('activ ss = 0.001', 'activ ss = 0.01', 'fresh ss = 0.001', 'fresh ss = 0.01')
title('k')
%axis([6 30 -100 100])
grid on
saveas(gcf,'k_nc.pdf')


