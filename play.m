%% RUNWRIGHTNET

clear;
net = getwrightnet( 1 );  % Condition 1
net.plot_every = 1;
net.neuron_to_plot = 2;
net.record_video = true;

n2_patts = find(mod(net.supplied_ts, 1000) == 8);
net.supplied_ts(n2_patts) = net.supplied_ts(n2_patts) + 3;
n2_sups = find(net.supplied_input == 4 & mod(net.supplied_ts, 1000) == 13);
net.supplied_input(n2_sups) = 5;

net.w_dend(2, :) = net.w_dend(1, :);
net.fgi = 4.0 + (43 * 0.2);

out = runsinglelayer(net);