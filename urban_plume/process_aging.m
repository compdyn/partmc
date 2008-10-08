load -ascii aging_a
load -ascii aging_p

time = aging_a(:,1);
activ = aging_a(:,4);
pure  = aging_p(:,4);
total = activ + pure;

lambda = 1.5e-5; % s^{-1}

activ_dot = (activ(60:end) - activ(1:(end-59))) ./ (time(60:end) - time(1:(end-59)));

activ_plot = activ(1:(end-59));
pure_plot = pure(1:(end-59));
time_plot = time(1:(end-59));

k_no_correct = activ_dot ./ pure_plot;
k = activ_dot ./ pure_plot + lambda * activ_plot ./ pure_plot;

tau_no_correct = 1.0 ./ k_no_correct;
tau = 1.0 ./ k;

plot(time/3600+6, activ, time/3600+6, pure, time/3600+6, total)
grid on
figure
plot(time_plot/3600+6, tau_no_correct/3600)
title('tau no correction')
axis([6 30 -100 100])
grid on
figure
plot(time_plot/3600+6, k*3600)
title('k')
grid on
figure
plot(time_plot/3600+6, tau/3600)
title('tau')
axis([6 30 -100 100])
grid on

