%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 
rand('seed', 1);
clear;

N_inp = 2000;
N_out = 3;
N = N_inp + N_out;
layer_sizes = [N_inp, N_out];
sim_time_sec = 300;
delay_max = 20;
num_dendrites = 500;
connection_matrix_size = [N_out, num_dendrites];
w_init = 0.65*4;
w_max = w_init * 1.5;
syn_mean_thresh = 0.6;

% Constants and conversions
ms_per_sec = 1000;

% Neuron parameters
v_rest = -65;
v_reset = -70;
v_thres = -55;
neuron_tau = 20;
v = v_rest * ones(N_out, 1);
w = ones(connection_matrix_size) * w_init;

% Connections
% pre is which neurons are pre-synaptic and post is which are post-synaptic
pre = randi([1 N_inp], connection_matrix_size);
delays = rand(connection_matrix_size) * delay_max;
last_spike_time = zeros(N, 1) * -Inf;
post = cell(N_inp, 1);
for n = 1 : N_inp
    post{n} = find(pre == n);
end

% Synapse dynamics parameters
g = zeros(size(v));
sigma_max = 10;
sigma_min = 0.1;
sigma = rand(connection_matrix_size) * (sigma_max - sigma_min) + sigma_min;

% STDP variables
taupre = 20;
taupost = 20;
Apre = 0.1;
Apost = -0.12;
STDPdecaypre = exp(-1/taupre);
STDPdecaypost = exp(-1/taupost);
dApre = zeros(connection_matrix_size);
dApost = zeros(connection_matrix_size);
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
%inp = [5, 2002, 5, 2002, 5, 2002, 2002, 2002, 5, 5, 5];
%ts = [200, 500, 550, 560, 590, 600, 601, 602, 615, 617, 618];
[ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp );

%% Main computation loop
for sec = 1 : sim_time_sec
    %% Reset second facilitating variables
    spike_times_trace = [];
    spike_arrival_trace = [];
    vt = zeros(size(vt));
    vt(:, 1) = v;
    
    [ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp, patt_inp, patt_ts );
    ts = ts + (sec-1) * 1000;
    tic;
    for ms = 1 : ms_per_sec
        time = (sec - 1) * ms_per_sec + ms;
        
        Iapp = zeros(size(v));
        t0 = time - reshape(last_spike_time(pre), size(pre));   % conn_matrix_size
        t0negu = t0 - delays;               % conn_matrix_size
        scale = 1 ./ (sigma .* sqrt(2 * pi));
        g = scale .* exp((-1/2) .* ((t0negu) ./ sigma) .^2 );
        g(isnan(g)) = 0;
        gaussian_values = w .* g;
        Iapp = sum(gaussian_values, 2);
        
        %% An incoming spike has arrived
        incoming = active_spikes{active_idx}; 
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            conn_idxs = incoming(:, 2);
            
            dApre(conn_idxs) = dApre(conn_idxs) + Apre;
        end
        
        %% Update membrane voltages
        v = v + (v_rest + Iapp - v) / neuron_tau;
        vt(:, ms) = v;
        
        %% Deal with neurons that just spiked
        fired = [find(v >=v_thres) + N_inp; inp(ts == time)'];
        spike_times_trace = [spike_times_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time;
        
        for spike = 1 : length(fired)
            neuron_idx = fired(spike);
            [neuron_layer, neuron_id] = idx2layerid(layer_sizes, neuron_idx);
            
            if neuron_layer ~= 1  %Output neuron
                v(neuron_id) = v_reset; % No lateral inhibition
                %v(:) = v_reset; % Lateral inhibition
                w(neuron_id, :) = w(neuron_id, :) + dApre(neuron_id, :);
                dApost(neuron_id, :) = dApost(neuron_id, :) + Apost;
                
                % Update SDVL
                presyn_idxs = pre(neuron_id, :);
                t0 = time - last_spike_time(presyn_idxs)';
                t0negu = t0 - delays(neuron_id, :);
                abst0negu = abs(t0negu);
                k = (sigma(neuron_id, :) + 0.9) .^ 2;
                shifts = sign(t0negu) .* k .* nu;
                
                % Update SDVL mean
                du = zeros(size(presyn_idxs));              % Otherwise
                du(t0 >= a2) = -k(t0 >= a2) .* nu;             % t0 >= a2
                du(abst0negu >= a1) = shifts(abst0negu >= a1); % |t0-u| >= a1
                
                delays(neuron_id, :) = delays(neuron_id, :) + du;
                delays = max(1, min(delay_max, delays));
                
                % Update SDVL variance
                dv = zeros(size(presyn_idxs));               % Otherwise
                dv(abst0negu < b2) = -k(abst0negu < b2) .* nv;  % |t0-u| < b2
                dv(abst0negu >= b1) = k(abst0negu >= b1) .* nv; % |t0-u| >= b1

                sigma(neuron_id, :) = sigma(neuron_id, :) + dv;
                sigma = max(sigma_min, min(sigma_max, sigma));
                
            else  % First layer (input)
                conn_idxs = post{neuron_idx};
                w(conn_idxs) = w(conn_idxs) + dApost(conn_idxs);
                
                paths = post{neuron_id};
                for path = 1 : numel(paths)
                    % Keep track of when this spike will arrive at each
                    % postsynaptic neuron and through which connection.
                    conn_idx = paths(path);
                    delay = round(delays(conn_idx));
                    arrival_offset = mod(active_idx + delay - 1, delay_max) + 1;
                    active_spikes{arrival_offset} = [active_spikes{arrival_offset}; neuron_id, conn_idx];
                end
            end
        end

        %% Update (decay) STDP variables
        dApre = dApre * STDPdecaypre;
        dApost = dApost * STDPdecaypost;
        active_idx = mod(active_idx, delay_max) + 1;
        
        % Synaptic scaling
        means = mean(w, 2);
        to_scale = means < syn_mean_thresh;
        if sum(to_scale) > 0
            w(to_scale, :) = w(to_scale, :) .* (syn_mean_thresh ./ means(to_scale));
        end 
        %debug = [debug; means(:)'];

        % Limit w to between [0, w_max]
        w = max(0, min(w_max, w)); 
        
        % Redistribute weak connections
        weak_conns = find(w < 0.05  & delays > 19.5);
        delays(weak_conns) = rand(size(weak_conns)) * delay_max;
        sigma(weak_conns) = rand(size(weak_conns)) * (sigma_max - sigma_min) + sigma_min;
        w(weak_conns) = w_init;
        dApre(weak_conns) = 0;
        dApost(weak_conns) = 0;
        last_spike_time(weak_conns) = -Inf;
        
        for c = 1 : numel(weak_conns)
            conn = weak_conns(c);
            old_pre = pre(conn);
            post{old_pre}(post{old_pre} == conn) = [];
            new_pre = randi([1 N_inp]);
            pre(conn) = new_pre;
            post{new_pre}(end + 1) = conn;
        end
        
        
        
    end
    
    %% Plot results from this second of processing     
    clf
    subplot(2, 1, 1);
    offset = (sec - 1) * 1000;
    filter = find(spike_times_trace(:,2) > N_inp & spike_times_trace(:,1) > offset);
    l2_spike_idxs = spike_times_trace(filter, 2);
    l2_spike_times = spike_times_trace(filter, 1);
    filter = find(spike_times_trace(:, 1) > offset);
    l1_idxs = spike_times_trace(filter, 2);
    l1_times = spike_times_trace(filter, 1);
    % Note this is the SPIKE TIME (not arrival time)
    plot(l1_times - offset, l1_idxs, '.k', 'MarkerSize', 8);
    hold on 
    plot(l2_spike_times - offset, l2_spike_idxs, '.r', 'MarkerSize', 8)
    ax = gca;
    i = 0;
    while i < numel(l2_spike_times)
        i = i + 1;
        pos = l2_spike_times(i);
        c_idx = mod(l2_spike_idxs(i) - N_inp - 1, size(ax.ColorOrder, 1)) + 1;
        colour = ax.ColorOrder(c_idx, :);
        plot( [pos pos] - offset, get( gca, 'Ylim' ), '--', 'Color', colour, 'LineWidth',2);
    end
    axis([0, 1000, 0 N]);
    xlabel('Time (ms)');
    ylabel('Neuron number');
    
    subplot(2, 1, 2);
    hist(pre');
    legend({'Blue - N1', 'Brown - N2', 'Yellow - N3'});
    title(sprintf('Synaptic distribution at second: %d', sec-1));
    xlabel('Pixel number');
    ylabel('Number of connections from each output in range')
    drawnow;
% %     plot(1:ms_per_sec, vt);
% %     title(sprintf('second: %d', sec-1));
% %     legend({'Neuron 1', 'Neuron 2', 'Neuron 3'});
% %     xlabel('Time (ms)');
% %     ylabel('Membrane potential (mV)');
% %     drawnow;
% %     
%     hold off;
   
    toc;
    fprintf('Second: %d\n', sec);
end