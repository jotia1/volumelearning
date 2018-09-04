function [ net ] = getwrightnet( cond )
%% GETWRIGHTNET - settings for Paul Wrights different conditions
%   Cond is the condition (experiment number) to use in reference to the
%   paper: 
%
%       Wright, Wiles (2012) "Learning Transmission Delays in Spiking
%       Neural Networks: A Novel Approach to Sequence Learning Based on 
%       Spike Delay Variance". 
%
%   Example usage:
%       net = getwrightnet(1);


%% Set up a default network structure
net = struct();

% Simulation run parameters
net.rand_seed = 1;
net.scaling_factor = 1;
net.w_init = 1;
net.w_max = 1;
net.syn_mean_thresh = Inf;
net.weak_con_thres = Inf;
net.plot_every = Inf;
net.lateral_inhibition_on = true;
net.synaptic_scaling_dend = false;
net.synaptic_scaling_axon = false;

% Neuron params
net.v_rest = -65;
net.v_reset = -70;
net.v_thres = -55;
net.neuron_tau = 20;

% Turn off STDP
net.taupre = 0;
net.taupost = 0;
net.Apre = 0;
net.Apost = 0;

net.inp_img_size = Inf;
net.num_dimensions_to_plot = 1;
net.neuron_to_plot = 1;

if ~exist('cond', 'var')
    cond = 1;
    disp('No condition supplied, assuming condition 1.');
end

if cond == 1
    
    % Network structure
    net.N_hid = 2;
    net.N_inp = 3;
    net.N_out = 3;
    
    net.sim_time_sec = 60;
    net.delay_max = 20;
    net.delay_min = 1;
    net.num_dendrites = 3;
    net.num_axons = 3;
    
    net.variance_dend = ones(net.N_hid, net.num_dendrites) * 2;
    net.w_dend = [];
    net.pre_dend = [];
    net.delays_dend = [];
    
    % SDVL params
    net.variance_max = 4;
    net.variance_min = 0.1;
    net.a1 = 1;
    net.a2 = 3;
    net.b1 = 8;
    net.b2 = 5;
    net.k = 1;
    net.nu = 0.0338;
    net.nv = 0.0218;
    
    % supply inp and ts below
    net.supplied_input = [4, 5, 6, 7, 8, 9 ];
    net.supplied_ts = [1, 5, 10, 15, 20, 25];
    
elseif cond == 2
    error('Not yet implemented');
elseif cond == 3
    error('Not yet implemented');
end

net.post_axon = [];
net.delays_axon = [];
net.w_axon = [];
net.variance_axon = [];

end
