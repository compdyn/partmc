load -ascii aging_a
load -ascii aging_p

x = aging_a(:,1)
y1 = aging_a(:,2)

y1_dot = 

tau = aging_p ./ aging_a

plot(x, tau)
