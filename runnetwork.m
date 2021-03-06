%% RUNNETWORK - Runs a default simulation

% net = getdefaultnet();
% net.input_source = 'D';
% net.input_filename = 'moving_day.aedat';
% net.sim_time_sec = 600;
% 
% net.num_axons = 32;%net.N_out;
% net.num_dendrites = 32;%net.N_inp;
% net.pre_dend = repmat(1:32, net.N_hid, 1);
% %net.delays_dend = ones(size(net.pre_dend));
% 
% % TODO: Possible error below, 33:64 should be 33 + N_hid + 1 : 64 + N_hid
% net.post_axon = repmat(33:64, net.N_hid, 1) + net.N_inp + net.N_hid;
% %net.delays_axon = ones(size(net.post_axon));
% 
% out = runsinglelayer(net);


%% Network to run to test 1D data

net = getdefaultnet();
net.input_source = 'S';

net.sim_time_sec = 500;
net.N_inp = 8;
net.N_hid = 1;
net.N_out = net.N_inp;
net.N = net.N_inp + net.N_hid + net.N_out;
net.N_v = net.N_hid + net.N_out;
net.scaling_factor = 1;
net.w_init = 16.95;
net.w_max = 29;
net.syn_mean_thresh = 1;
net.weak_con_thres = 10;
net.plot_every =  25; % Never

net.Apre = 0.1;
net.Apost = -0.12;
net.taupre = 20;
net.taupost = 20;
net.synaptic_scaling_dend = true;
net.synaptic_scaling_axon = true;

net.num_axons = 8;
net.num_dendrites = 8;
[ net.supplied_input, net.supplied_ts ] = gen1Ddata(net.sim_time_sec, ...
                                            net.num_axons, 8, 0, 0);
net.num_dimensions_to_plot = 1;
                                        
% Fully connected network
net.pre_dend = repmat(1:net.num_axons, net.N_hid, 1);
net.post_axon = repmat(1:net.num_dendrites, ...
    net.N_hid, 1) + net.N_inp + net.N_hid;

net.neuron_to_plot = 1;

% SDVL off;
net.nu = 0;
net.nv = 0;
net.variance_max = 0.1;
net.variance_min = 0.1;
net.delay_max = 1;

out = runsinglelayer(net);


















