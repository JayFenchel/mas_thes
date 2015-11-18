N = [5, 10, 20, 30, 40];
% durchschnittliche Rechenzeiten
% cond = [89.5, 500, 3586.65, 11193];
% sparse = [457, 1620 5273, 11358];
cond = [152, 668, 4313, 12946, 25132];
sparse = [531, 1633 5584, 11437, 19797];

% plot
figure(1)
subplot(2, 1, 1)
hold on;
plot(N, cond);
plot(N, sparse, 'r');
grid on
ylabel('Rechenzeit t [ms]')
subplot(2, 1, 2)
plot(N, sparse./cond)
ylabel('t_{sparse}/t_{cond} [ms/ms]')
xlabel('Praediktionshorizont T_p []')
grid on

N = [5, 10, 20, 30];
% vorher
cond1 = [100, 520, 4100, 13400];
sparse1 = [950, 4250, 27500, 80000];

% nachher
cond2 = [70, 295, 1200, 3750];
sparse2 = [585, 1620, 5500, 11200];

% % plot
% figure(1)
% subplot(2, 1, 1)
% hold on;
% plot(N, cond1);
% plot(N, sparse1, 'r');
% subplot(2, 1, 2)
% plot(N, sparse1./cond1)
% figure(2)
% subplot(2, 1, 1)
% hold on;
% plot(N, cond2);
% plot(N, sparse2, 'r');
% subplot(2, 1, 2)
% plot(N, sparse2./cond2)
figure(3)
subplot(2, 1, 1)
hold on;
plot(N, cond1);
plot(N, sparse2, 'r');
ylabel('Rechenzeit [ms]')
subplot(2, 1, 2)
plot(N, sparse2./cond1)
ylabel('sparse/cond []')
xlabel('Praediktionshorizont []')