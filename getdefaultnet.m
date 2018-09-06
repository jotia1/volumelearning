function [ default_net ] = getdefaultnet( )

%% Set up a default network structure
default_net = struct();

% Network structure
default_net.inp_img_size = 8;
default_net.N_hid = 3;
default_net.N_inp = default_net.inp_img_size * default_net.inp_img_size;
default_net.N_out = default_net.N_inp;

% define connections
% Pass empty matrix to random initialise.
default_net.pre_dend = [];
default_net.delays_dend = [];
default_net.post_axon = [];
default_net.delays_axon = [];
default_net.w_dend = [];
default_net.w_axon = [];
default_net.variance_dend = [];
default_net.variance_axon = [];

% Simulation run parameters
default_net.rand_seed = 2;
default_net.sim_time_sec = 5;
default_net.delay_max = 1;
default_net.delay_min = 1;
default_net.num_dendrites = 32;
default_net.num_axons = 32;
default_net.scaling_factor = 40;
default_net.w_init = 0.65 * default_net.scaling_factor;
default_net.w_max = default_net.w_init * 1.1;
default_net.syn_mean_thresh = default_net.w_init * 0.6;
default_net.weak_con_thres = default_net.w_init * 0.1;
default_net.plot_every = 1;
default_net.lateral_inhibition_on = false;
default_net.synaptic_scaling_dend = true;
default_net.synaptic_scaling_axon = true;

% Neuron params
default_net.v_rest = -65;
default_net.v_reset = -70;
default_net.v_thres = -55;
default_net.neuron_tau = 20;
default_net.izhikevich_neurons = false;

% SDVL params
default_net.variance_max = 4;
default_net.variance_min = 0.1;
default_net.a1 = 1;
default_net.a2 = 3;
default_net.b1 = 8;
default_net.b2 = 5;
%default_net.k = 1;
default_net.nu = 0.0338;
default_net.nv = 0.0218;

default_net.taupre = 16;
default_net.taupost = 33;
default_net.Apre = 0.1 * default_net.scaling_factor;
default_net.Apost = -0.12 * default_net.scaling_factor;

default_net.num_dimensions_to_plot = 2;
default_net.neuron_to_plot = 1;

% supply inp and ts below
default_net.supplied_input = [4, 5, 6, 7, 8, 9 ];
default_net.supplied_ts = [1, 5, 10, 15, 20, 25];
end