load -ascii aging_a_wc
load -ascii aging_p_wc

time = aging_a_wc(:,1);
activ_001 = aging_a_wc(:,2);
pure_001  = aging_p_wc(:,2);
total_001 = activ_001 + pure_001;
activ_01 = aging_a_wc(:,11);
pure_01  = aging_p_wc(:,11);
total_01 = activ_01 + pure_01;

lambda = 1.5e-5; % s^{-1}

activ_001_dot = (activ_001(60:end) - activ_001(1:(end-59))) ./ (time(60:end) - time(1:(end-59)));
activ_01_dot = (activ_01(60:end) - activ_01(1:(end-59))) ./ (time(60:end) - time(1:(end-59)));

activ_001_plot = activ_001(1:(end-59));
pure_001_plot = pure_001(1:(end-59));

activ_01_plot = activ_01(1:(end-59));
pure_01_plot = pure_01(1:(end-59));

time_plot = time(1:(end-59));

%k_no_correct = activ_dot ./ pure_plot;
k_001 = activ_001_dot ./ pure_001_plot + lambda * activ_001_plot ./ pure_001_plot;
k_01 = activ_01_dot ./ pure_01_plot + lambda * activ_01_plot ./ pure_01_plot;

%tau_no_correct = 1.0 ./ k_no_correct;
tau_001 = 1.0 ./ k_001;
tau_01  = 1.0 ./ k_01;

transfer_001 = activ_001_dot + lambda * activ_001_plot;
transfer_01 = activ_01_dot + lambda * activ_01_plot;

plot(time/3600+6, activ_001, time/3600+6, pure_001, time/3600+6, total_001, time/3600+6, activ_01, time/3600+6, pure_01, time/3600+6, total_01)
legend('activating ss = 0.001', 'not-activating ss = 0.001', 'total', 'activating ss = 0.01', 'not-activating ss = 0.01', 'total')
grid on
saveas(gcf,'number_wc.pdf')
%figure
%plot(time_plot/3600+6, tau_no_correct/3600)
%title('tau no correction')
%axis([6 30 -100 100])
%grid on
%figure
%plot(time_plot/3600+6, k*3600)
%title('k')
%grid on
figure
plot(time_plot/3600+6, tau_001/3600, time_plot/3600+6, tau_01/3600)
legend('ss = 0.001', 'ss = 0.01')
title('tau')
axis([6 30 -100 100])
grid on
saveas(gcf,'tau_wc.pdf')

figure
plot(time_plot/3600+6, k_001/3600, time_plot/3600+6, k_01/3600)
legend('ss = 0.001', 'ss = 0.01')
title('k')
%axis([6 30 -100 100])
grid on
saveas(gcf,'k_wc.pdf')

figure
plot(time_plot/3600+6, transfer_001, time_plot/3600+6, transfer_01)
legend('ss = 0.001', 'ss = 0.01')
title('transfer')
%axis([6 30 -100 100])
grid on
saveas(gcf,'transfer_wc.pdf')

