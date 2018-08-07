function [ default_net ] = getdefaultnet( )

%% Set up a default network structure
default_net = struct();

% Network structure
default_net.inp_img_size = 8;
default_net.N_hid = 3;
default_net.N_inp = default_net.inp_img_size * default_net.inp_img_size;
default_net.N_out = default_net.N_inp;
default_net.N = default_net.N_inp + default_net.N_hid + default_net.N_out;
default_net.N_v = default_net.N_hid + default_net.N_out;

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
default_net.rand_seed = 1;
default_net.sim_time_sec = 5;
default_net.delay_max = 20;
default_net.num_dendrites = 32;
default_net.num_axons = 32;
default_net.scaling_factor = 150;
default_net.w_init = 0.65 * default_net.scaling_factor;
default_net.w_max = default_net.w_init * 1.5;
default_net.syn_mean_thresh = default_net.w_init * 0.8;
default_net.weak_con_thres = default_net.w_init * 0.1;
default_net.plot_every = 1;
default_net.lateral_inhibition_on = false;
default_net.synaptic_scaling_dend = true;

% Neuron params
default_net.v_rest = -65;
default_net.v_reset = -70;
default_net.v_thres = -55;
default_net.neuron_tau = 20;

% SDVL params
default_net.variance_max = 10;
default_net.variance_min = 0.1;
default_net.a1 = 1;
default_net.a2 = 3;
default_net.b1 = 8;
default_net.b2 = 5;
default_net.k = 1;
default_net.nu = 0.0338;
default_net.nv = 0.0218;

default_net.taupre = 20;
default_net.taupost = 20;
default_net.Apre = 0.1 * default_net.scaling_factor;
default_net.Apost = -0.12 * default_net.scaling_factor;

% Data params
% Options are 'G': Generated, 'D': DVS, 'S': Supplied
default_net.input_source = 'S';  
    % if using DVS, fill in below:
    default_net.input_filename = '';
    default_net.input_size = default_net.inp_img_size;
    default_net.input_max_size = 0; % 0 is all.
    % if supplying inp and ts fill in below
    default_net.supplied_input = [4, 5, 6, 7, 8, 9 ];
    default_net.supplied_ts = [1, 5, 10, 15, 20, 25];
end