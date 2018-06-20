%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 
rand('seed', 1);
clear;

% Often tweaked parameters
N_inp = 2000;
N_out = 3;
N = N_inp + N_out;
sim_time_sec = 300;
delay_max = 20;
num_connections = 3;
w_init = 0.65;
w_max = w_init * 1.5;
syn_mean_thresh = w_init * 0.65;

assert(w_init > syn_mean_thresh, 'SS will interfere');
assert(w_init < w_max, 'Careful synaptic scaling will limit weights');

% Constants and conversions
ms_per_sec = 1000;

% Neuron parameters
v_rest = -65;
v_reset = -70;
v_thres = -55;
neuron_tau = 20;
v = v_rest * ones(N, 1); 
w = ones(N, num_connections) * w_init;

% Synpatic links 
% Post is who is postsynaptic from a pre. N x M, M is number of connections
post = randi([N_inp + 1, N], N, num_connections);

% Synapse dynamics parameters
g = zeros(N, 1);
sigma_max = 10;
sigma_min = 0.1;
sigma = rand(N, num_connections) * (sigma_max - sigma_min) + sigma_min;
delays = ceil(rand(N, num_connections) * delay_max);
last_spike_time = zeros(N, 1) * -Inf;

% SDVL variables
a1 = 1;         % 
a2 = 3;         % 
b1 = 8;         % 
b2 = 5;         % 
k = 1;          % Learning accelerator (scaling factor)
nu = 0.0338;    % Learning rate (for the mean)
nv = 0.0218;    % Learning rate (for the variance

% STDP variables
taupre = 20;
taupost = 20;
Apre = 0.1;
Apost = -0.12;
STDPdecaypre = exp(-1/taupre);
STDPdecaypost = exp(-1/taupost);
% dApre/post represent the decayed activity of when pre/post synaptic
% neurons fired. It is the amount to change weights by.
dApre = zeros(N, num_connections);
dApost = zeros(N, 1);
active_spikes = cell(delay_max, 1);  % To track when spikes arrive
active_idx = 1;

% Synaptic Redistribution variables
% TODO

% Info logging variables
%spike_arrival_trace = [];
spike_times_trace = [];
vt = zeros(N, ms_per_sec);
vt(:, 1) = v;
debug = [];

%% DATA
[ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp );
% [xs, ys, ts, ps] = loadDVSsegment(16, false, 0);
% inp = sub2ind([16 16], xs, ys);
% ts = floor(ts /1000);

%% Main computation loop
for sec = 1 : sim_time_sec
    [ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp, patt_inp, patt_ts );
    ts = ts + (sec-1) * 1000;
%     if mod(sec+1, 32) == 0  %cheeky shift data to keep training
%         ts = ts + 32*1000;
%     end
    tic;
    for ms = 1 : ms_per_sec
        time = (sec - 1) * ms_per_sec + ms;  
        
        %% Calculate input to neurons at this time step
        Iapp = zeros(N, 1);
        t0 = time - last_spike_time;
        t0negu = t0 - delays;
        scale = 1 ./ (sigma .* sqrt(2 * pi));
        g = scale .* exp((-1/2) .* ((t0negu) ./ sigma) .^2 );
        g(isnan(g)) = 0;
        gaussian_values = w .* g;
      
        % Collect input for each neuron based on synapses facing them
        % TODO can be optimised with a precalculated array
        [to_neurons, conn_pre_from, to_neurons_idx] = unique(post);
        Iapp(to_neurons) = accumarray(to_neurons_idx, gaussian_values(:));

        %% An incoming spike arrived at a neuron
        incoming = active_spikes{active_idx}; 
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            % incoming(:,1) is pre-syn neuron, incoming(:, 2) is connection
            incoming_idxs = sub2ind(size(dApre), incoming(:, 1), incoming(:,2));
            
            % Update STDP variables
            dApre(incoming_idxs) = dApre(incoming_idxs) + Apre;
  
%             Log spike arrivals for plotting later
%             spike_arrival_trace = [spike_arrival_trace; 
%                                    time*ones(length(incoming),1), incoming];
        end

        
        %% Update membrane voltages
        v = v + (v_rest + Iapp - v) / neuron_tau;
        vt(:, ms) = v;
        
        %% Deal with neurons that just spiked
        fired = [find(v >=v_thres); inp(ts == time)'];
        spike_times_trace = [spike_times_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time;  % Used in synapse dynamics 
        
        for spike = 1 : length(fired)
            neuron = fired(spike);
            v(neuron) = v_reset;
            
            % Update STDP variables
            presynaptic_idxs = find(post == neuron);   
            % I spiked, reinforce any before me
            w(presynaptic_idxs) = w(presynaptic_idxs) + dApre(presynaptic_idxs);
            % weaken connections to any after me who spiked not long ago
            w(neuron, :) = w(neuron, :) + dApost(post(neuron, :))'; 
            
            dApost(neuron) = dApost(neuron) + Apost;
            
            % Update SVDL
            if neuron > N_inp
                v(N_inp + 1: end, :) = v_reset;
                [presyn_neurons, ~] = ind2sub(size(post), presynaptic_idxs);
                t0 = time - last_spike_time(presyn_neurons);
                t0negu = t0 - delays(presynaptic_idxs);
                abst0negu = abs(t0negu);
                k = (sigma(presyn_neurons) + 0.9) .^ 2;
                shifts = sign(t0negu) .* k .* nu;

                % Update SDVL mean
                du = zeros(size(presyn_neurons));              % Otherwise
                du(t0 >= a2) = -k(t0 >= a2) .* nu;             % t0 >= a2
                du(abst0negu >= a1) = shifts(abst0negu >= a1); % |t0-u| >= a1
                
                delays(presynaptic_idxs) = delays(presynaptic_idxs) + du;
                delays = max(1, min(delay_max, delays));

                % Update SDVL variance
                dv = zeros(size(presyn_neurons));               % Otherwise
                dv(abst0negu < b2) = -k(abst0negu < b2) .* nv;  % |t0-u| < b2
                dv(abst0negu >= b1) = k(abst0negu >= b1) .* nv; % |t0-u| >= b1

                sigma(presynaptic_idxs) = sigma(presynaptic_idxs) + dv;
                sigma = max(sigma_min, min(sigma_max, sigma));
              
            end
            
            for connection = 1:num_connections
                % Keep track of when this spike will arrive at each
                % postsynaptic neuron and through which connection.
                delay = round(delays(connection));
                arrival_offset = mod(active_idx + delay - 1, delay_max) + 1;
                active_spikes{arrival_offset} = [active_spikes{arrival_offset}; neuron, connection];
            end
        end
        
        %% Update (decay) STDP variables
        dApre = dApre * STDPdecaypre;
        dApost = dApost * STDPdecaypost;
        active_idx = mod(active_idx, delay_max) + 1;
        
        %% Apply synaptic scaling (SS) and weight bounding
        %output_means = accumarray(to_neurons_idx, w(:), [], @mean);
        %to_scale = output_means < syn_mean_thresh; 
        
        
        %w(conn_pre_from) = w(conn_pre_from)
        
        
        
        %debug = [debug; output_means'];
        
        w_mean = mean(w(:));
        % TODO this scales weights across ALL outputs not each individually
        if w_mean < syn_mean_thresh
            w = w .* (syn_mean_thresh / w_mean);
        end
        
        % Limit w to between [0, w_max]
        w = max(0, min(w_max, w)); 
        w(N_inp + 1 : end, :) = 0;
        
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
    xlabel('Time (ms)');
    ylabel('Neuron number');
    
    
    subplot(2, 1, 2);
    plot(1:ms_per_sec, vt(N_inp + 1:end, :));
    title(sprintf('second: %d', sec-1));
    legend({'Neuron 1', 'Neuron 2', 'Neuron 3'});
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    drawnow;
    
    hold off;
    
    %% Reset second facilitating variables
    spike_times_trace = [];
    spike_arrival_trace = [];
    vt = zeros(N, ms_per_sec);
    vt(:, 1) = v;
    toc;
    fprintf('Second: %d\n', sec);
    %waitforbuttonpress;
end