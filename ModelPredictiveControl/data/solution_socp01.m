x4 = [-400.000,
        -399.789,
        -394.682,
        -382.518,
        -366.265,
        -347.807,
        -328.126,
        -307.756,
        -286.986,
        -266.041,
        -244.993,
        -223.838,
        -202.578,
        -181.221,
        -159.752,
        -138.172,
        -116.469,
        -94.650,
        -72.759,
        -50.927,
        -30.118,
        -12.032,
        1.855,
        10.862,
        15.160,
        15.572,
        13.270,
        9.493,
        5.329,
        1.587,
        -1.249,
        -2.993,
        -3.700,
        -3.577,
        -2.904,
        -1.963,
        -0.993,
        -0.163,
        0.435,
        0.774,
        0.880,
        0.808,
        0.625,
        0.397,
        0.174,
        -0.007,
        -0.131,
        -0.194,
        -0.206,
        -0.180,
        -0.133,
        -0.078,
        -0.028,
        0.011,
        0.036,
        0.048,
        0.048,
        0.040,
        0.028,
        0.015];
x2 = 180/pi*[0.000,
        0.226,
        0.408,
        0.404,
        0.399,
        0.398,
        0.397,
        0.381,
        0.377,
        0.371,
        0.372,
        0.368,
        0.367,
        0.366,
        0.358,
        0.342,
        0.333,
        0.292,
        0.202,
        0.097,
        0.005,
        -0.055,
        -0.080,
        -0.076,
        -0.054,
        -0.026,
        -0.010,
        -0.003,
        0.002,
        0.008,
        0.012,
        0.014,
        0.015,
        0.014,
        0.013,
        0.010,
        0.008,
        0.005,
        0.002,
        -0.000,
        -0.002,
        -0.003,
        -0.005,
        -0.006,
        -0.005,
        -0.004,
        -0.002,
        -0.001,
        0.000,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.001,
        0.000,
        -0.000,
        -0.000,
        -0.000,
        -0.000];
% load('/home/jayf/Arbeitsfläche/Ergebnisse QPcondensed/QPcond_n30_i6.mat')
load('/home/jayf/Downloads/acc2016/mpcat/mpcctl/systems/syspaper/QPCONDTEST500reg0001_n40_i7.mat')
t_plot = [0:0.5:29.5];
zwanni = 20.01*ones(60);
figure(19)
for i=1:10
% subplot(1, 2, 2)
%     hold on
%     grid on
%     plot(t_plot, x4(i, :), 'b')
%     ylabel('x_4 [m]')
%     xlabel('time [s]')
% subplot(1, 2, 1)
%     hold on
%     grid on
%     plot(t_plot, 180/pi*x2(i, :), 'b')
%     plot(t_plot, zwanni, 'r--')
%     xlabel('time [s]')
%     ylabel('x_2 [deg]')
end
% load('/home/jayf/Downloads/acc2016/mpcat/mpcctl/systems/syspaper/QPpce01.mat')
% t_plot = [0:0.5:29.5];
% zwanni = 20.01*ones(60);
% figure(8)
% for i=1:100
% subplot(3, 2, 4)
%     hold on
%     grid on
%     plot(t_plot, x4(i, :), 'b')
%     ylabel('x_4 [m]')
%     xlabel('time [s]')
% subplot(3, 2, 3)
%     hold on
%     grid on
%     plot(t_plot, 180/pi*x2(i, :), 'b')
%     plot(t_plot, zwanni, 'b--')
%     xlabel('time [s]')
%     ylabel('x_2 [deg]')
% end
% load('/home/jayf/Downloads/acc2016/mpcat/mpcctl/systems/syspaper/SOCPpce01.mat')
% t_plot = [0:0.5:29.5];
% zwanni = 20.01*ones(60);
% figure(8)
% for i=1:100
% subplot(3, 2, 6)
%     hold on
%     grid on
%     plot(t_plot, x4(i, :), 'b')
%     ylabel('x_4 [m]')
%     xlabel('time [s]')
% subplot(3, 2, 5)
%     hold on
%     grid on
%     plot(t_plot, 180/pi*x2(i, :), 'b')
%     plot(t_plot, zwanni, 'b--')
%     xlabel('time [s]')
%     ylabel('x_2 [deg]')
% end
sprintf('%16.8f', m_cost)
sprintf('%16.8f', m_viol)
sprintf('%16.8f', m_dt)