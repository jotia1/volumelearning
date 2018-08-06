function [ default_net ] = getdefaultnet( )

%% Set up a default network structure
default_net = struct();

% Network sizes
default_net.inp_img_size = 8;
default_net.N_hid = 3;
default_net.N_inp = default_net.inp_img_size * default_net.inp_img_size;

% Simulation run parameters
default_net.sim_time_sec = 50;
default_net.delay_max = 20;
default_net.num_dendrites = 32;
default_net.num_axons = 32;
default_net.scaling_factor = 50;
default_net.w_init = 0.65 * default_net.scaling_factor;
default_net.w_max = default_net.w_init * 1.5;
default_net.syn_mean_thresh = default_net.w_init * 0.8;
default_net.weak_con_thres = default_net.w_init * 0.1;
default_net.plot_every = 1;

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

default_net.tau_pre = 20;
default_net.tau_post = 20;
default_net.Apre = 0.1 * default_net.scaling_factor;
default_net.Apost = -0.12 * default_net.scaling_factor;

% Data params
default_net.input_source = 'GENERATED';  % Or: DVS
    % if using DVS, fill in below:
    default_net.input_filename = '';
    default_net.input_size = default_net.inp_img_size;
    default_net.input_max_size = 0; % 0 is all.
end