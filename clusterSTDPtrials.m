%% CLUSTERSTDPTRIALS

net = getdefaultnet();
net.input_source = 'S';

net.sim_time_sec = 100;
net.N_inp = 8;
net.N_hid = 1;
net.N_out = net.N_inp;
net.N = net.N_inp + net.N_hid + net.N_out;
net.N_v = net.N_hid + net.N_out;
net.scaling_factor = 1;
% net.w_init = 16.95;
% net.w_max = 29;
% net.syn_mean_thresh = 1;
net.weak_con_thres = 10;
net.plot_every = net.sim_time_sec + 1; % Never

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

% Variables to tweak
net.w_init = 16.95;  % [5-25]
net.w_max = 29;      % [10-30]   
net.syn_mean_thresh = 1; % [0 - 20]

% Get a valid folder to dump results to
count = 0;
output_folder = 'out0';
while exist(output_folder, 'dir') == 7
    count = count + 1;
    output_folder = sprintf('out%02d', count);
end

mkdir(output_folder);

for w_init = 5:25
    for w_max = 10:30
        for smt = 0:20
            net.w_init = w_init;
            net.w_max = w_max;
            net.synaptic_mean_thresh = smt;

            out = runsinglelayer(net);
            save(sprintf('%s/w_init-%d_w_max-%d_smt-%d', output_folder, w_init, w_max, smt), 'net', 'out');
        end
    end
end

















