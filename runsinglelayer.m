function [ output ] = runsinglelayer( net )
%% RUNSINGLELAYER - Takes a network data struct and runs network
%   Given a valid network description as defined by validatenetwork
%   description, run that network and calculate the output.

%% Connectivity explanations
% pre_dend [mat] - is the pixel ID which is presynaptic to a connection in 
%       the dendrite layer. Indexed by connection.
% post_dend [cell] - Dendrite connection idx(s) for a given pixel from the 
%       input layer. Indexed by input layer pixel ID.
% delays_dend [mat] - the delay along a given hidden connection. (Idexed by
%       connection id.
% pre_axon [cell] - Gives the axon connection ID(s) of neurons presynaptic 
%       to the given pixel ID. Indexed by pixel ID.
% post_axon [mat] - Gives the pixel index postsynaptic of a connection.
%       Indexed by connection ID. 
% delays_axon [mat] - Gives the delay along a particular axon connection.
%       Indexed by connection ID. 



output = struct();
timing_info = struct();
output.timing_info.init_time = tic;
assert(validatenetwork(net), 'Error with network description');
rand('seed', net.rand_seed);

layer_sizes = [net.N_inp, net.N_hid, net.N_out];
dend_conn_matrix_size = [net.N_hid, net.num_dendrites];
axon_conn_matrix_size = [net.N_hid, net.num_axons];

neuron_start_idx = net.N_inp + 1;
axon_start_idx = neuron_start_idx + net.N_hid;

v = net.v_rest * ones(net.N_v, 1);

%% Connections
w_dend = net.w_dend;
if isempty(w_dend)
    w_dend = ones(dend_conn_matrix_size) * net.w_init;
end
w_axon = net.w_axon;
if isempty(w_axon)
    w_axon = ones(axon_conn_matrix_size) * net.w_init;
end

pre_dend = net.pre_dend; 
if isempty(pre_dend) 
    pre_dend = randi([1 net.N_inp], dend_conn_matrix_size);
end
delays_dend = net.delays_dend;
if isempty(delays_dend)
    delays_dend = rand(dend_conn_matrix_size) * net.delay_max;
end
post_dend = cell(net.N_inp, 1);
for n = 1 : net.N_inp
    post_dend{n} = find(pre_dend == n);
end

% output connections
post_axon = net.post_axon;
if isempty(post_axon)
    post_axon = randi([axon_start_idx, net.N], axon_conn_matrix_size);
end
delays_axon = net.delays_axon;
if isempty(delays_axon)
    delays_axon = rand(axon_conn_matrix_size) * net.delay_max;
end
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
variance_dend = net.variance_dend;
if isempty(variance_dend)
    variance_dend = rand(dend_conn_matrix_size) * variance_range;
end
variance_axon = net.variance_axon;
if isempty(variance_axon)
    variance_axon = rand(axon_conn_matrix_size) * variance_range;
end

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

neuron_to_plot = net.neuron_to_plot;

%% Load data
output.timing_info.load_data = tic;
%% Data
if net.input_source == 'G'
    [ inp, ts, patt_inp, patt_ts ] = embedPat( net.N_inp );
elseif net.input_source == 'D'
    % DVS data
    [xs, ys, ts, ps] = loadDVSsegment(128, false, 0);
    [ xs, ys, ts, ps ] = dvs2patch( xs, ys, ts, ps, net.inp_img_size, net.sim_time_sec, 30 );
    rows = net.inp_img_size - ys + 1; cols = net.inp_img_size - xs + 1;
    inp = sub2ind([net.inp_img_size net.inp_img_size], rows, cols);
    ts = floor(ts / 1000);
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
        t0_axon = time - last_spike_time(neuron_start_idx:axon_start_idx - 1);
        [~, active_conns] = idx2layerid(layer_sizes, post_axon); % convert to id's
        t0negu_axon = repmat(t0_axon, 1, net.num_axons) - delays_axon;
        scale_axon = 1 ./ (variance_axon .* sqrt(2 * pi));
        g_axon = scale_axon .* exp((-1/2) .* ((t0negu_axon) ./ variance_axon) .^2 );
        g_axon(isnan(g_axon)) = 0;
        gaussian_values_axon = w_axon .* g_axon;

        [to_neurons, ~, to_neurons_idx] = unique(active_conns);  %TODO fix speed
        if mean(gaussian_values_axon(:)) > 0
            disp('');
        end
        Iapp(net.N_hid + to_neurons) = accumarray(to_neurons_idx, gaussian_values_axon(:));
   
        if mean(Iapp(:)) > 0
            disp('')
        end
        
        %% A spike has arrived do STDP
        incoming = active_spikes{active_idx}; 
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            from_inp = incoming(:,1) <= net.N_inp;
            inp_spike_idxs = incoming(from_inp, 2);
            dApre_dend(inp_spike_idxs) = dApre_dend(inp_spike_idxs) + net.Apre;
            
            hid_spike_idxs = incoming(:,1) < axon_start_idx & ~from_inp;
            dApre_axon(hid_spike_idxs) = dApre_axon(hid_spike_idxs) + net.Apre;
            if sum(dApre_axon(:)) > 0
                disp('');
            end
            
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
        
        for spike = 1 : length(fired)
            neuron_idx = fired(spike);
            [neuron_layer, neuron_id] = idx2layerid(layer_sizes, neuron_idx);
        
            if neuron_layer == 3 % output layer
                hid_conn_idxs = pre_axon{neuron_id};
                v(neuron_id + net.N_hid) = net.v_reset;
                
                % Update STDP

                w_axon(hid_conn_idxs) = w_axon(hid_conn_idxs) + dApre_axon(hid_conn_idxs);
                dApost_axon(hid_conn_idxs) = dApost_axon(hid_conn_idxs) + net.Apost;
                
                
                
                
                % Update SDVL
                [hid_ids, ~] = ind2sub(size(post_axon), hid_conn_idxs);
                pre_idxs = net.N_inp + hid_ids;  % TODO should really use layerid2idx here
                t0_axon = time - last_spike_time(pre_idxs);
                t0negu_axon = t0_axon - delays_axon(hid_conn_idxs);
                abst0negu_axon = abs(t0negu_axon);
                k = (variance_axon(hid_conn_idxs) + 0.9) .^ 2;
                shifts = sign(t0negu_axon) .* k .* net.nu;
                
                % Update SDVL mean
                du = zeros(size(hid_conn_idxs));
                du(t0_axon >= net.a2) = -k(t0_axon >= net.a2) .* net.nu;
                du(abst0negu_axon >= net.a1) = shifts(abst0negu_axon >= net.a1);% TODO: verify this line is correct, made an edit without checkign the maths.
                
                delays_axon(hid_conn_idxs) = delays_axon(hid_conn_idxs) + du;
                delays_axon = max(1, min(net.delay_max, delays_axon));
                
                % Update SDVL variance
                dv = zeros(size(hid_conn_idxs));
                dv(abst0negu_axon < net.b2) = -k(abst0negu_axon < net.b2);
                dv(abst0negu_axon >= net.b1) = k(abst0negu_axon >= net.b1);
                
                variance_axon(hid_conn_idxs) = variance_axon(hid_conn_idxs) + dv;
                variance_axon = max(net.variance_min, min(net.variance_max, variance_axon));
            
                
            elseif neuron_layer == 2  %hid neuron
                v([neuron_id, find(net.lateral_inhibition_on)*1:net.N_hid]) = net.v_reset;
                
                % Update STDP of dendrites
                w_dend(neuron_id, :) = w_dend(neuron_id, :) + dApre_dend(neuron_id, :);
                dApost_dend(neuron_id, :) = dApost_dend(neuron_id, :) + net.Apost;
                
                % Update STDP of axons
                w_axon(neuron_id, :) = w_axon(neuron_id, :) + dApost_axon(neuron_id, :);
                % Set up active spikes (to later adjust dApre)
                conn_delays = round(delays_axon(neuron_id, :));
                arrival_offsets = mod(active_idx + conn_delays - 1, net.delay_max) + 1;
                [to_offsets, ~, offset_idxs] = unique(arrival_offsets);
                for i = 1 : numel(to_offsets)
                    to_offset = to_offsets(i);
                    neuron_local_idxs = find(offset_idxs == i);
                    conn_idxs = sub2ind(size(delays_axon), ones(size(neuron_local_idxs)) * neuron_id, neuron_local_idxs);
                    active_spikes{to_offset} = [active_spikes{to_offset}; ones(size(conn_idxs)) * neuron_idx, conn_idxs];
                end
                
                % Update SDVL
                presyn_idxs = pre_dend(neuron_id, :);
                t0_dend = time - last_spike_time(presyn_idxs)';
                t0negu_dend = t0_dend - delays_dend(neuron_id, :);
                abst0negu_dend = abs(t0negu_dend);
                k = (variance_dend(neuron_id, :) + 0.9) .^ 2;
                shifts = sign(t0negu_dend) .* k .* net.nu;
                
                % Update SDVL mean
                du = zeros(size(presyn_idxs));              % Otherwise
                du(t0_dend >= net.a2) = -k(t0_dend >= net.a2) .* net.nu;             % t0 >= a2
                du(abst0negu_dend >= net.a1) = shifts(abst0negu_dend >= net.a1); % |t0-u| >= a1
                
                delays_dend(neuron_id, :) = delays_dend(neuron_id, :) + du;
                delays_dend = max(1, min(net.delay_max, delays_dend));
                
                % Update SDVL variance
                dv = zeros(size(presyn_idxs));               % Otherwise
                dv(abst0negu_dend < net.b2) = -k(abst0negu_dend < net.b2) .* net.nv;  % |t0-u| < b2
                dv(abst0negu_dend >= net.b1) = k(abst0negu_dend >= net.b1) .* net.nv; % |t0-u| >= b1

                variance_dend(neuron_id, :) = variance_dend(neuron_id, :) + dv;
                variance_dend = max(net.variance_min, min(net.variance_max, variance_dend));


            else  % First layer (input)
                conn_idxs = post_dend{neuron_id};
                % penalise connection(s) if the postsynaptic spiked recently.
                w_dend(conn_idxs) = w_dend(conn_idxs) + dApost_dend(conn_idxs);
                
                paths = post_dend{neuron_id};
                for path = 1 : numel(paths)
                    % Keep track of when this spike will arrive at each
                    % postsynaptic neuron and through which connection.
                    conn_idx = paths(path);
                    delay = round(delays_dend(conn_idx));
                    arrival_offset = mod(active_idx + delay - 1, net.delay_max) + 1;
                    active_spikes{arrival_offset} = [active_spikes{arrival_offset}; neuron_idx, conn_idx];
                end
            end
        end
        
        %% Update (decay) STDP variables
        dApre_dend = dApre_dend * STDPdecaypre;
        dApost_dend = dApost_dend * STDPdecaypost;
        dApre_axon = dApre_axon * STDPdecaypre;
        dApost_axon = dApost_axon * STDPdecaypost;
        active_idx = mod(active_idx, net.delay_max) + 1;
        
        % Dendritic synaptic scaling
        means = mean(w_dend, 2);
        to_scale = net.synaptic_scaling_dend & means < net.syn_mean_thresh;
        if sum(to_scale) > 0
            w_dend(to_scale, :) = w_dend(to_scale, :) .* (net.syn_mean_thresh ./ means(to_scale));
        end 
        
        % Axonal synaptic scaling
        means = mean(w_axon, 2);
        to_scale = net.synaptic_scaling_axon & means < net.syn_mean_thresh;
        if sum(to_scale) > 0
            w_axon(to_scale, :) = w_axon(to_scale, :) .* (net.syn_mean_thresh ./ means(to_scale));
            %        neuron_to_plot = 3;
            %visualiseweights;
            %disp(time);
        end         
        
        % Synaptic bounding - limit w to [0, w_max]
        w_dend = max(0, min(net.w_max, w_dend)); 
        w_axon = max(0, min(net.w_max, w_axon));
        %neuron_to_plot = 3;
        %visualiseweights;
        %disp('');
        
%         % Redistribute weak connections
%         % TODO - add in redistribution for output
%         weak_conns = find(net.synaptic_redistribution_on & w_hid < weak_con_thres);%  & delays_hid > 19.5);
%         delays_hid(weak_conns) = rand(size(weak_conns)) * delay_max;
%         variance_hid(weak_conns) = rand(size(weak_conns)) * (variance_max - variance_min) + variance_min;
%         w_hid(weak_conns) = w_init;
%         dApre_hid(weak_conns) = 0;
%         dApost_hid(weak_conns) = 0;
%         last_spike_time(weak_conns) = -Inf;
% 
%         for c = 1 : numel(weak_conns)
%             conn = weak_conns(c);
%             old_pre = pre_hid(conn);
%             post_hid{old_pre}(post_hid{old_pre} == conn) = [];
%             new_pre = randi([1 N_inp]);
%             pre_hid(conn) = new_pre;
%             post_hid{new_pre}(end + 1) = conn;
%         end  
    end
    output.timing_info.sim_sec_tocs(sec) = toc(output.timing_info.sim_sec_times(sec));
    
    %% Plotting
    if mod(sec, net.plot_every) == 0
        output.timing_info.plotting_tics(end + 1) = tic;
                
        if net.num_dimensions_to_plot == 2
            visualise2Dweights; 
        elseif net.num_dimensions_to_plot == 1
            visualise1Dweights;
        end
        

        output.timing_info.plotting_tocs(end + 1) = toc(output.timing_info.plotting_tics(end));
    end
    
    output.spike_times_trace = [output.spike_times_trace; spike_times_trace]; % TODO should speed this up (expanding list).

    
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

















