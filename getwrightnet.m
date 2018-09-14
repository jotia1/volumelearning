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
net.w_max = net.w_init;
net.syn_mean_thresh = Inf;
net.weak_con_thres = Inf;
net.lateral_inhibition_on = true;
net.synaptic_scaling_dend = false;
net.synaptic_scaling_axon = false;

% Neuron params
net.v_rest = -65;
net.v_reset = -65;
net.v_thres = 30;
net.neuron_tau = 20;
net.izhikevich_neurons = true;

% Turn off STDP
net.taupre = 0;
net.taupost = 0;
net.Apre = 0;
net.Apost = 0;

net.inp_img_size = Inf;
net.num_dimensions_to_plot = 1;
net.neuron_to_plot = 1;
net.plot_every = Inf;
net.record_video = false;

if ~exist('cond', 'var')
    cond = 1;
    disp('No condition supplied, assuming condition 1.');
end

if cond == 1
    
    % Network structure
    net.N_hid = 2;
    net.N_inp = 3;
    net.N_out = net.N_hid;
    net.num_dendrites = net.N_inp;
    net.num_axons = 1;
    
    net.sim_time_sec = 300;
    net.delay_max = 15;
    net.delay_min = 1;
    
    dend_conn_matrix_size = [net.N_hid, net.num_dendrites];
    net.variance_dend = ones(dend_conn_matrix_size) * 2;
    net.w_dend = ones(dend_conn_matrix_size) * 1;
    net.w_dend(2, :) = 0;
    net.pre_dend = repmat(1:net.N_inp, net.N_hid, 1);
    net.delays_dend = ones(dend_conn_matrix_size) * 5;
    
    % SDVL params
    net.variance_max = 10;
    net.variance_min = 0.1;
    net.a1 = 3;
    net.a2 = 2;
    net.b1 = 5;
    net.b2 = 5;
    net.nu = 0.03;
    net.nv = 0.01;
    net.fgi = 13;
    
    % supply inp and ts below
    seq = [0, 3, 7];
    net.supplied_input = repmat(1:3, 1, net.sim_time_sec * 2);
    net.supplied_ts = reshape(repmat(((0:(net.sim_time_sec * 2) -1 )' * 500)', 3, 1), 1, net.sim_time_sec * 3 * 2) + ...
        repmat(seq, 1, net.sim_time_sec * 2) + 1;
    
    % Set up supervision
    net.supervised_seconds = 100;
    net.supplied_input = [net.supplied_input, ones(1, net.supervised_seconds * 2) * net.N_inp + 1];
    net.supplied_ts = [net.supplied_ts, ((0:0.5:net.supervised_seconds - 0.5) * 1000) + 13];
    
    
elseif cond == 2
    error('Not yet implemented');
elseif cond == 3
    error('Not yet implemented');
end

net.post_axon = (net.N_inp + net.N_hid + 1 : net.N_inp + net.N_hid + net.N_out)';
net.delays_axon = ones(net.N_hid, net.num_axons);
net.w_axon = ones(net.N_hid, net.num_axons);
net.variance_axon = ones(net.N_hid, net.num_axons);

end
