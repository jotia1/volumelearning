function [ output ] = runsinglelayer( net )
%% RUNSINGLELAYER - Takes a network data struct and runs network
%   Given a valid network description as defined by validatenetwork
%   description, run that network and calculate the output.

%% Connectivity explanations
% pre_dend [mat] - is the pixel ID which is presynaptic to a connection in 
%       the dendrite layer. Indexed by connection.
% post_dend [cell] - Neuron ID(s) that are postsynaptic to a given pixel
%       from the input layer. Indexed by input layer pixel ID.
% delays_dend [mat] - the delay along a given hidden connection. (Idexed by
%       connection id.
% pre_axon [cell] - Gives the axon connection ID(s) of neurons presynaptic 
%       to the given pixel ID. Indexed by pixel ID.
% post_axon [mat] - Gives the pixel ID postsynaptic of a connection.
%       Indexed by connection ID. 
% delays_axon [mat] - Gives the delay along a particular axon connection.
%       Indexed by connection ID. 



output = struct();
timing_info = struct();
output.timing_info.init_time = tic;
assert(validatenetwork(net), 'Error with network description');

layer_sizes = [net.N_inp, net.N_hid, net.N_out];
dend_conn_matrix_size = [net.N_hid, net.num_dendrites];
axon_conn_matrix_size = [net.N_hid, net.num_axons];

neuron_start_idx = net.N_inp + 1;
axon_start_idx = neuron_start_idx + net.N_hid;

v = net.v_rest * ones(net.N_v, 1);
w_dend = ones(dend_conn_matrix_size) * net.w_init;
w_axon = ones(axon_conn_matrix_size) * net.w_init;

% Connections
pre_dend = randi([1 net.N_inp], dend_conn_matrix_size);
delays_dend = rand(dend_conn_matrix_size) * net.delay_max;
post_dend = cell(net.N_inp, 1);
for n = 1 : net.N_inp
    post_dend{n} = find(pre_dend == n);
end

% output connections
delays_axon = rand(axon_conn_matrix_size) * net.delay_max;
post_axon = randi([axon_start_idx, net.N], axon_conn_matrix_size);
pre_axon = cell(net.N_out, 1);
for n = 1 : net.N_out
    n_idx = n + axon_start_idx - 1;
    pre_axon{n} = find(post_axon == n_idx);
end

ms_per_sec = 1e3;
last_spike_time = zeros(net.N, 1) * -Inf;

% Synapse dynamics parameters
g_dend = zeros(net.N_hid, 1);
g_axon = zeros(net.N_out, 1);
variance_range = (net.variance_max - net.variance_min) + net.variance_min;
variance_dend = rand(dend_conn_matrix_size) * variance_range;
variance_axon = rand(axon_conn_matrix_size) * variance_range;

% STDP variables
STDPdecaypre = exp(-1/net.taupre);
STDPdecaypost = exp(-1/net.taupost);
dApre_dend = zeros(dend_conn_matrix_size);
dApost_dend = zeros(dend_conn_matrix_size);
dApre_axon = zeros(axon_conn_matrix_size);
dApost_axon = zeros(axon_conn_matrix_size);
active_spikes = cell(net.delay_max, 1);  % To track when spikes arrive
active_idx = 1;

% Output data
spike_times_trace = [];
output.spike_times_trace = [];
vt = zeros(size(v, 1), ms_per_sec);
vt(:, 1) = v;

%% Load data
output.timing_info.load_data = tic;
%% Data
if net.input_source == 'G'
    [ inp, ts, patt_inp, patt_ts ] = embedPat( net.N_inp );
elseif net.input_source == 'D'
    % DVS data
    [xs, ys, ts, ps] = loadDVSsegment(128, false, 0);
    [ xs, ys, ts, ps ] = dvs2patch( xs, ys, ts, ps, 8, 50, 30 );
    inp = sub2ind([net.inp_img_size net.inp_img_size], xs, ys);
    ts = floor(ts /1000);
elseif net.input_source == 'S'
    inp = net.supplied_input;
    ts = net.supplied_ts;
else
    error(['Unknown input source: ', net.input_source, ' (try G).']);
end


output.timing_info.load_data = toc(output.timing_info.load_data);
output.timing_info.init_time = toc(output.timing_info.init_time);
output.timing_info.sim_sec_tics = uint64(zeros(net.sim_time_sec, 1));
output.timing_info.sim_sec_tocs = zeros(net.sim_time_sec, 1);
output.timing_info.plotting_tics = uint64([]);
output.timing_info.plotting_tocs = [];

for sec = 1 : net.sim_time_sec
    
    output.timing_info.sim_sec_times(sec) = tic;
    output.spike_times_trace = [output.spike_times_trace; spike_times_trace]; % TODO should speed this up (expanding list).
    spike_times_trace = [];
    vt = zeros(size(vt));
    vt(:, 1) = v;
    
    for ms = 1 : ms_per_sec
        
        time = (sec - 1) * ms_per_sec + ms;
        
        %% Caculate input at this step
        Iapp = zeros(size(v));
        % Dendrites
        t0_dend = time - reshape(last_spike_time(pre_dend), size(pre_dend));   % conn_matrix_size
        t0negu_dend = t0_dend - delays_dend;               % conn_matrix_size
        scale_dend = 1 ./ (variance_dend .* sqrt(2 * pi));
        g_dend = scale_dend .* exp((-1/2) .* ((t0negu_dend) ./ variance_dend) .^2 );
        g_dend(isnan(g_dend)) = 0;
        gaussian_values_dend = w_dend .* g_dend;
        Iapp(1:net.N_hid, :) = sum(gaussian_values_dend(1:net.N_hid, :), 2);
        % Axons
        t0_axon = time - last_spike_time(axon_start_idx:end);
        [~, active_conns] = idx2layerid(layer_sizes, post_axon);
        t0negu_axon = t0_axon(active_conns) - delays_axon;
        scale_axon = 1 ./ (variance_axon .* sqrt(2 * pi));
        g_axon = scale_axon .* exp((-1/2) .* ((t0negu_axon) ./ variance_axon) .^2 );
        g_axon(isnan(g_axon)) = 0;
        gaussian_values_axon = w_axon .* g_axon;

        [to_neurons, ~, to_neurons_idx] = unique(active_conns);  %TODO fix speed
        Iapp(net.N_hid + to_neurons) = accumarray(to_neurons_idx, gaussian_values_axon(:));
   
        %% A spike has arrived
        incoming = active_spikes{active_idx}; 
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            from_inp = incoming(:,1) <= net.N_inp;
            inp_spike_idxs = incoming(from_inp, 2);
            dApre_dend(inp_spike_idxs) = dApre_dend(inp_spike_idxs) + net.Apre;
            
            hid_spike_idxs = incoming(~from_inp, 2);
            dApre_axon(hid_spike_idxs) = dApre_axon(hid_spike_idxs) + net.Apre;
            
        end
        %% Update membrane voltages
        v = v + (net.v_rest + Iapp - v) / net.neuron_tau;
        vt(:, ms) = v;
        
        
        %% Deal with neurons that just spiked
        % TODO: Searching for ts == time is slow if events is large which
        % it is for real data (searching whole stream). Could trim only
        % feeding 1 sec to inp per second.
        fired_pixels = inp(ts == time);
        fired = [find(v >=net.v_thres) + net.N_inp; fired_pixels'; net.N_inp + net.N_hid + inp(ts == time  & inp <= net.N_inp)'];  % TODO Hack to inject spikes
        spike_times_trace = [spike_times_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time; 
        
    end
    output.timing_info.sim_sec_tocs(sec) = toc(output.timing_info.sim_sec_times(sec));
    
    %% Plotting
    if mod(sec, net.plot_every) == 0
        output.timing_info.plotting_tics(end + 1) = tic;
        %visualiseweights(); % TODO: turn into function
       
        output.timing_info.plotting_tocs(end + 1) = toc(output.timing_info.plotting_tics(end));
    end
    
end

%% Collect final state of variables
output.pre_dend = pre_dend;
output.post_dend = post_dend;
output.delays_dend = delays_dend;
output.variance_dend = variance_dend;
output.pre_axon = pre_axon;
output.post_axon = post_axon;
output.delays_axon = delays_axon;
output.variance_axon = variance_axon;

output.dApre_dend = dApre_dend;
output.dApost_dend = dApost_dend;
output.dApre_axon = dApre_axon;
output.dApost_axon = dApost_axon;

% Clean output
output.timing_info = rmfield(output.timing_info, {'sim_sec_tics', 'plotting_tics'});

end

















