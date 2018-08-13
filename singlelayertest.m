% SINGLELAYERTEST - A suite of tests to verify functionality of networks
net = getdefaultnet();
% Tests to write still
%   - Test for STDP in axons/dendrites (both ways)
%   - Test for SDVL in axons/dendrites (both ways), both mean and variance
%   adjusting.
%   - Test for synaptic redistribution (on and off & axon and dendrite)
%   - Test for synaptic scaling (on and off & axon and dendrite)
%   - Test for weight bounding (over w_max and under min)
%   - Test for delay bounding (over delay_max and under 1)
%   - Test for variance bounding, above and below max/min

% preconditions
assert(validatenetwork(net), 'Default network has an error');

%%  Test 1: Dendrite connections
net.sim_time_sec = 1;
% Every connection connects between pixel 1 and all hidden neurons.
net.pre_dend = ones(net.N_hid, net.num_dendrites);
% Pixel 2 now connects to hid neuron 1 through connection 1
net.pre_dend(1, 1) = 2;
% Pixel 3 now connects to hid neuron 2 through last (num_dendrites)
% connection
net.pre_dend(2, net.num_dendrites) = 3;
% Pixel (N_Inp) connects to hid neuron 3 through last connection 
net.pre_dend(3, net.num_dendrites) = net.N_inp - 7;
net.pre_dend(3, 5:9) = sub2ind([8,8], [2, 3 ,4 5, 5], [2, 2, 2, 2, 3]);

% Make each spike a big deal (arrives after 1ms and a LOT at once).
net.variance_dend = ones(net.N_hid, net.num_dendrites) * net.variance_min;
net.delays_dend = ones(net.N_hid, net.num_dendrites); 

%
rows = repmat([2, 3, 1, 2, 3, 4, 5, 5], 1, 1);
cols = repmat([1, 1, 8, 2, 2, 2, 2, 3], 1, 1);

% Make specified pixels spike to induce spikes in hidden layer
net.supplied_input = sub2ind([net.inp_img_size net.inp_img_size], rows, cols);
net.supplied_ts = reshape(repmat(1:1, 8, 1) * 10, 1, []); 

assert(validatenetwork(net), 'Test 1 net invalid');
out = runsinglelayer(net);
assert(sum(find(out.spike_times_trace == 1 + net.N_inp)) > 0, 'No connections to hid neuron 1');
assert(sum(find(out.spike_times_trace == 2 + net.N_inp)) > 0, 'No connections to hid neuron 2');
assert(sum(find(out.spike_times_trace == 3 + net.N_inp)) > 0, 'No connections to hid neuron 3');


%% Test 2: Axon connections
% Check that neurons correctly feed into axons and to output.
net.sim_time_sec = 1;
% All axons go to pixel id 2
axons_start_idx = net.N_inp + net.N_hid + 1;
net.post_axon = ones(net.N_hid, net.num_axons) * (axons_start_idx + 1);
% neuron 1 connects to pixel id 1 through connection 1;
net.post_axon(1, 1) = axons_start_idx;
% neuron 2 connects to last pixel through connection 1
net.post_axon(2, 1) = net.N;
% neuron 3 connects to second last pixel through last connection
net.post_axon(3, net.num_axons) = net.N - 1;

% Make each spike a big deal (arrives after 1ms and a LOT at once).
net.variance_axon = ones(net.N_hid, net.num_axons) * net.variance_min;
net.delays_axon = ones(net.N_hid, net.num_axons); 

% Make neurons spike to induce spikes in output layer
net.supplied_input = repmat(net.N_inp + 1:net.N_inp + net.N_hid, 1, 1);
net.supplied_ts = reshape(repmat(1:1, 3, 1) * 10, 1, []); 

assert(validatenetwork(net), 'Test 2 net invalid');
out = runsinglelayer(net);
assert(sum(find(out.spike_times_trace == axons_start_idx)) > 0, 'No connections from hid neuron 1');
assert(sum(find(out.spike_times_trace == net.N)) > 0, 'No connections from hid neuron 2');
assert(sum(find(out.spike_times_trace == net.N - 1)) > 0, 'No connections from hid neuron 3');






