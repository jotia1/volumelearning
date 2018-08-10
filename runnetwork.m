%% RUNNETWORK - Runs a default simulation

net = getdefaultnet();
net.input_source = 'D';
net.input_filename = 'moving_day.aedat';
net.sim_time_sec = 30;

% Lets make a fully connected network
net.num_axons = net.N_out;
net.num_dendrites = net.N_inp;
net.pre_dend = repmat(1:net.N_inp, net.N_hid, 1);
%net.delays_dend = ones(size(net.pre_dend));
net.post_axon = repmat(1:net.N_out, net.N_hid, 1) + net.N_inp + net.N_hid;
%net.delays_axon = ones(size(net.post_axon));


out = runsinglelayer(net);