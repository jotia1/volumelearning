%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 

% Note on terminology: neuron_idx refers to the neurons index globally and
% each neuron has a unique index which is sequentially assigned. neuron_id
% refers to the neurons sequential index within its layer and thus several
% neurons may have the same id (and be from different layers) but all have
% unique index's. 
rand('seed', 1);
clear;

im_size = 16;
N_inp = im_size*im_size;
N_hid = 2;
N_out = N_inp;
N = N_inp + N_hid + N_out;
N_v = N_hid + N_out;
layer_sizes = [N_inp, N_hid, N_out];
sim_time_sec = 1;
delay_max = 20;
num_dendrites = 6;%floor(N_inp / N_hid);  % TODO - testing hack
num_axons = 8;%floor(N_inp / N_hid);        % TODO - testing hack
hid_conn_matrix_size = [N_hid, num_dendrites];
out_conn_matrix_size = [N_hid, num_axons];
scaling_factor = 150;                       % TODO - testing hack
w_init = 0.65 * scaling_factor;
w_max = w_init * 1.5;
syn_mean_thresh = w_init * 0.8;
weak_con_thres = w_init*0.1;
plot_every = 1;

% Constants and conversions
ms_per_sec = 1000;

% Neuron parameters
v_rest = -65;
v_reset = -70;
v_thres = -55;
neuron_tau = 20;
v = v_rest * ones(N_v, 1);
w_hid = ones(hid_conn_matrix_size) * w_init;
w_out = ones(out_conn_matrix_size) * w_init;
%w = [w_hid; w_out];

% Connections
% pre is which neurons are pre-synaptic and post is which are post-synaptic
pre_hid = randi([1 N_inp], hid_conn_matrix_size);
pre_hid(1, :) = [8, 7, 256, 6, 5, 4];                          % TODO - testing hack
delays_hid = rand(hid_conn_matrix_size) * delay_max;
delays_hid(1, :) = [1, 1, 1, 1, 1, 1];                      % TODO - testing hack
post_hid = cell(N_inp, 1);
for n = 1 : N_inp
    post_hid{n} = find(pre_hid == n);
end
% output connections
%pre_out = randi([1 N_hid] + N_inp, out_conn_matrix_size);
delays_out = rand(out_conn_matrix_size) * delay_max;
post_out = randi([N_inp + N_hid + 1, N], out_conn_matrix_size);
post_out(1, :) = [8, 7, 256, 6, 5, 4, 3, 2] + N_inp + N_hid;               % TODO - testing hack
pre_out = cell(N_out, 1);
for n = 1 : N_out
    n_idx = n + N_hid + N_inp;
    % Note: pre_out is the pre-synaptic connection index's 
    pre_out{n} = find(post_out == n_idx);
    %[layers, ids] = ind2sub(size(post_out), find(post_out == n_idx));
    %pre_out{n} = layers;
end
%pre_out = repmat([1:N_out]', 1, num_axons) + N_inp + N_hid; 
% post_out = cell(N_hid, 1);
% for n = 1 : N_hid
%     post_out{n} = find(pre_out == n);
% end
%delays = [delays_hid; delays_out];
%pre = [pre_hid; pre_out];
%post = [post_hid; post_out];
last_spike_time = zeros(N, 1) * -Inf;

% Synapse dynamics parameters
g_hid = zeros(N_hid, 1);
g_out = zeros(N_out, 1);
variance_max = 10;
variance_min = 0.1;
variance_hid = rand(hid_conn_matrix_size) * (variance_max - variance_min) + variance_min;
variance_hid(1, :) = [.4, .5, .6, 2, 2, 2]; 
variance_out = rand(out_conn_matrix_size) * (variance_max - variance_min) + variance_min;


% STDP variables
taupre = 20;
taupost = 20;
Apre = 0.1 * scaling_factor;
Apost = -0.12 * scaling_factor;
STDPdecaypre = exp(-1/taupre);
STDPdecaypost = exp(-1/taupost);
dApre_hid = zeros(hid_conn_matrix_size);
dApost_hid = zeros(hid_conn_matrix_size);
dApre_out = zeros(out_conn_matrix_size);
dApost_out = zeros(out_conn_matrix_size);
active_spikes = cell(delay_max, 1);  % To track when spikes arrive
active_idx = 1;

% SDVL variables
a1 = 1;         % 
a2 = 3;         % 
b1 = 8;         % 
b2 = 5;         % 
k = 1;          % Learning accelerator (scaling factor)
nu = 0.0338;    % Learning rate (for the mean)
nv = 0.0218;    % Learning rate (for the variance

%Logging
spike_times_trace = [];
vt = zeros(size(v, 1), ms_per_sec);
vt(:, 1) = v;
debug = [];

%% Data
inp = [8, 7, 256, 6, 5, 4 ];
ts = [1, 5, 10, 15, 20, 25];
%[ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp );

% DVS data
%[xs, ys, ts, ps] = loadDVSsegment(im_size, false, 0);
%inp = sub2ind([im_size im_size], xs, ys);
%ts = floor(ts /1000);

%% Main computation loop
for sec = 1 : sim_time_sec
    
    %% Reset second facilitating variables
    spike_times_trace = [];
    spike_arrival_trace = [];
    vt = zeros(size(vt));
    vt(:, 1) = v;
    
    %[ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp, patt_inp, patt_ts );
    %ts = ts + (sec-1) * 1000;
    inp_length_sec = floor(ts(end) - ts(1) / 1000);
    if mod(sec, inp_length_sec) == 0  %cheeky shift data to keep training
         ts = ts + inp_length_sec*1000;
    end
    tic; 
    %avgs = zeros(1000, 4);
    for ms = 1 : ms_per_sec
        %tic;
        time = (sec - 1) * ms_per_sec + ms;
        
        Iapp = zeros(size(v));
        % Hidden layer
        t0_hid = time - reshape(last_spike_time(pre_hid), size(pre_hid));   % conn_matrix_size
        t0negu_hid = t0_hid - delays_hid;               % conn_matrix_size
        scale_hid = 1 ./ (variance_hid .* sqrt(2 * pi));
        g_hid = scale_hid .* exp((-1/2) .* ((t0negu_hid) ./ variance_hid) .^2 );
        g_hid(isnan(g_hid)) = 0;
        gaussian_values_hid = w_hid .* g_hid;
        Iapp(1:N_hid, :) = sum(gaussian_values_hid(1:N_hid, :), 2);
        % Output layer
        t0_out = time - last_spike_time(N-N_out+1:end);
        [~, active_conns] = idx2layerid(layer_sizes, post_out);
        t0negu_out = t0_out(active_conns) - delays_out;
        scale_out = 1 ./ (variance_out .* sqrt(2 * pi));
        g_out = scale_out .* exp((-1/2) .* ((t0negu_out) ./ variance_out) .^2 );
        g_out(isnan(g_out)) = 0;
        gaussian_values_out = w_out .* g_out;
        
        %avgs(ms, 1) = toc;
        %tic; 
        
        % Collect input for each neuron based on synapses facing them
        % TODO can be optimised with a precalculated array, Yes this is
        % slow
        [to_neurons, conn_pre_from, to_neurons_idx] = unique(active_conns);  %TODO fix
        Iapp(to_neurons + N_hid) = accumarray(to_neurons_idx, gaussian_values_out(:));
        
        %avgs(ms, 2) = toc;
        %tic;
        
%         if time > 500
%            disp('')
%         end
                
        %% A spike has arrived
        incoming = active_spikes{active_idx}; 
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            inp_idxs = incoming(incoming(:,1) <= N_inp, 2);
            dApre_hid(inp_idxs) = dApre_hid(inp_idxs) + Apre;
            
            %out_idxs = incoming(incoming(:,1) > N_inp, 2);
            %dApre_out(out_idxs) = dApre_out(out_idxs) + Apre;
        end
        %avgs(ms, 3) = toc;
        %tic;
        
        %% Update membrane voltages
        v = v + (v_rest + Iapp - v) / neuron_tau;
        vt(:, ms) = v;

        %% Deal with neurons that just spiked
        % TODO: Searching for ts == time is slow if events is large which
        % it is for real data (searching whole stream). Could trim only
        % feeding 1 sec to inp per second.
        fired_pixels = inp(ts == time);
        fired = [find(v >=v_thres) + N_inp; fired_pixels'; N_inp + N_hid + inp(ts == time  & inp <= N_inp)'];  % TODO Hack to inject spikes
        spike_times_trace = [spike_times_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time; 
        
        %avgs(ms, 4) = toc;
        
        for spike = 1 : length(fired)
            neuron_idx = fired(spike);
            if neuron_idx == 516
               disp 
            end
            [neuron_layer, neuron_id] = idx2layerid(layer_sizes, neuron_idx);
            
            if neuron_layer ~= 3 % Anything thats not input
                
            end
            
            if neuron_layer == 3 % output layer
                %v(neuron_id) = v_reset; % No lateral inhibition
                % TODO - Split into hid and out
                %pre_idxs = find(post_out == n_idx);
                if time > 14
                   [[time;time], w_out];
                end
                hid_idxs = pre_out{neuron_id};  % The CONNECTION(s) pre to this
                % Maybe pre_out should give back the connection_idxs
                % Update STDP
                w_out(hid_idxs) = w_out(hid_idxs) + dApre_out(hid_idxs);
                dApost_out(hid_idxs) = dApost_out(hid_idxs) + Apost;
                
                % Output has spiked, reinforce any who spiked just before
                % TODO
                
                
            
            elseif neuron_layer == 2  %hid neuron
                v(neuron_id) = v_reset; % No lateral inhibition
                %v(1:N_hid) = v_reset; % Lateral inhibition
                w_hid(neuron_id, :) = w_hid(neuron_id, :) + dApre_hid(neuron_id, :);
                dApost_hid(neuron_id, :) = dApost_hid(neuron_id, :) + Apost;
                % Update STDP of output connections
                dApre_out(neuron_id, :) = dApre_out(neuron_id, :) + Apre;
                
                % Update SDVL
                presyn_idxs = pre_hid(neuron_id, :);
                t0_hid = time - last_spike_time(presyn_idxs)';
                t0negu_hid = t0_hid - delays_hid(neuron_id, :);
                abst0negu = abs(t0negu_hid);
                k = (variance_hid(neuron_id, :) + 0.9) .^ 2;
                shifts = sign(t0negu_hid) .* k .* nu;
                
                % Update SDVL mean
                du = zeros(size(presyn_idxs));              % Otherwise
                du(t0_hid >= a2) = -k(t0_hid >= a2) .* nu;             % t0 >= a2
                du(abst0negu >= a1) = shifts(abst0negu >= a1); % |t0-u| >= a1
                
                delays_hid(neuron_id, :) = delays_hid(neuron_id, :) + du;
                delays_hid = max(1, min(delay_max, delays_hid));
                
                % Update SDVL variance
                dv = zeros(size(presyn_idxs));               % Otherwise
                dv(abst0negu < b2) = -k(abst0negu < b2) .* nv;  % |t0-u| < b2
                dv(abst0negu >= b1) = k(abst0negu >= b1) .* nv; % |t0-u| >= b1

                variance_hid(neuron_id, :) = variance_hid(neuron_id, :) + dv;
                variance_hid = max(variance_min, min(variance_max, variance_hid));
                
                % STDP - hid (dendrites)
                % I spiked, penalise any post-synaptics that have spikes
                % just before this
                w_out(neuron_id, :) = w_out(neuron_id, :) + dApost_out(neuron_id, :);
                
                for d = 1 : delay_max
                    % TODO up to here, was selecting all connections with a
                    % given delay and then would append them all to
                    % acitve_spikes to arrive later. 
                    conn_delays = delays_out(neuron_id, delays_out(neuron_id, :) == d);
                    arrival_offset = mod(active_idx + delay - 1, delay_max) + 1;
                end  
                
            else  % First layer (input)
                conn_idxs = post_hid{neuron_idx};
                w_hid(conn_idxs) = w_hid(conn_idxs) + dApost_hid(conn_idxs);
                
                paths = post_hid{neuron_id};
                for path = 1 : numel(paths)
                    % Keep track of when this spike will arrive at each
                    % postsynaptic neuron and through which connection.
                    conn_idx = paths(path);
                    delay = round(delays_hid(conn_idx));
                    arrival_offset = mod(active_idx + delay - 1, delay_max) + 1;
                    active_spikes{arrival_offset} = [active_spikes{arrival_offset}; neuron_idx, conn_idx];
                end
            end
        end
        
        %% Update (decay) STDP variables
        dApre_hid = dApre_hid * STDPdecaypre;
        dApost_hid = dApost_hid * STDPdecaypost;
        dApre_out = dApre_out * STDPdecaypre;
        dApost_out = dApost_out * STDPdecaypost;
        active_idx = mod(active_idx, delay_max) + 1;
        
        % Synaptic scaling
        means = mean(w_hid, 2);
        to_scale = means < syn_mean_thresh;
        if sum(to_scale) > 0
            w_hid(to_scale, :) = w_hid(to_scale, :) .* (syn_mean_thresh ./ means(to_scale));
        end 
        % TODO - Synaptic scaling for output connections??
        %debug = [debug; means(:)'];
        
        % Limit w to between [0, w_max]
        w_hid = max(0, min(w_max, w_hid)); 
        
        % Redistribute weak connections
        % TODO - add in redistribution for output
        weak_conns = find(w_hid < weak_con_thres);%  & delays_hid > 19.5);
        delays_hid(weak_conns) = rand(size(weak_conns)) * delay_max;
        variance_hid(weak_conns) = rand(size(weak_conns)) * (variance_max - variance_min) + variance_min;
        w_hid(weak_conns) = w_init;
        dApre_hid(weak_conns) = 0;
        dApost_hid(weak_conns) = 0;
        last_spike_time(weak_conns) = -Inf;
        
        for c = 1 : numel(weak_conns)
            conn = weak_conns(c);
            old_pre = pre_hid(conn);
            post_hid{old_pre}(post_hid{old_pre} == conn) = [];
            new_pre = randi([1 N_inp]);
            pre_hid(conn) = new_pre;
            post_hid{new_pre}(end + 1) = conn;
        end
        
    end
    %floor(mean(avgs, 1) * 1000000)
    
    %% Plot results from this second of processing
    if mod(sec, plot_every) == 0
        %clf
        subplot(2, 2, 1);
        hold off
        offset = (sec - 1) * 1000;
        filter = find(spike_times_trace(:,2) <= N_inp + N_hid & spike_times_trace(:,2) > N_inp & spike_times_trace(:,1) > offset);
        l2_spike_idxs = spike_times_trace(filter, 2);
        l2_spike_times = spike_times_trace(filter, 1);
        filter = find(spike_times_trace(:,2) < N_inp + N_hid & spike_times_trace(:, 1) > offset);
        l1_idxs = spike_times_trace(filter, 2);
        l1_times = spike_times_trace(filter, 1);
        % Note this is the SPIKE TIME (not arrival time)
        plot(l1_times - offset, l1_idxs, '.k', 'MarkerSize', 8);
        hold on 
        plot(l2_spike_times - offset, l2_spike_idxs, '.r', 'MarkerSize', 8)
        ax = gca;
        i = 0;
        axis([0, 30, -30 N_v + 50]);
        while i < numel(l2_spike_times)
            i = i + 1;
            pos = l2_spike_times(i);
            c_idx = mod(l2_spike_idxs(i) - N_inp - 1, size(ax.ColorOrder, 1)) + 1;
            colour = ax.ColorOrder(c_idx, :);
            plot( [pos pos] - offset, get( gca, 'Ylim' ), '--', 'Color', colour, 'LineWidth',2);
        end

        xlabel('Time (ms)');
        ylabel('Neuron number');

        subplot(4, 2, 5);
        hist(pre_hid');
        %legend({'Blue - N1', 'Brown - N2', 'Yellow - N3'});
        title(sprintf('Synaptic distribution at second: %d', sec-1));
        xlabel('Pixel number');
        ylabel('Number of connections from each output in range');
        drawnow;
        
        subplot(4, 2, 2);
        hist(w_hid');
        title('weights\_hid');
        
        subplot(4, 2, 4);
        hist(delays_hid');
        title('delays\_hid');
        
        subplot(4, 2, 6);
        hist(variance_hid');
        title('variance\_hid')
        
        drawnow;
        
        subplot(4, 2, 7);
        plot(1:ms_per_sec, vt(1:N_hid,:));
        axis([0 30 -70 -50]);
        title(sprintf('second: %d', sec-1));
        %legend({'Neuron 1', 'Neuron 2', 'Neuron 3'});
        xlabel('Time (ms)');
        ylabel('Membrane potential (mV)');
        hold off;
        drawnow;
        
        %waitforbuttonpress;
    end
   
    toc;
    fprintf('Second: %d\n', sec);
end