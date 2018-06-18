%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 

% Features
%   - 
%

rand('seed', 1);
%clear;
% Often tweaked parameters
N_inp = 2000;
N_out = 1;
N = N_inp + N_out;
sim_time_sec = 25;
delay_max = 1;
num_connections = 3;
w_init = 0.45;
w_max = 1;
syn_mean_thresh = w_init / 2;

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
post = ones(N, num_connections) * N; 
post(700, 1) = N-1; %TODO - Create actual connectivity, this is just a demo
post(800:804, :) = N-2;

% Synapse dynamics parameters
g = zeros(N, 1);
sigma_max = 7;
sigma_min = 2;  % TODO where did i get all these variables from...
sigma = rand(N, 1) * (sigma_max - sigma_min) + sigma_min;
delays = ceil(rand(N, num_connections) * delay_max);
last_spike_time = zeros(N, 1) * -Inf;

% STDP variables
taupre = 20;       % NOTE ABOUT STDP WITH DATA EXAMPLE:
taupost = 20;      % This example is somewhat cherry picked, small
Apre = 0.1;        % changes result in large destabilisations of learning
Apost = -0.12;     % Something is still wrong.
STDPdecaypre = exp(-1/taupre);
STDPdecaypost = exp(-1/taupost);
% dApre/post represent the decayed activity of when pre/post synaptic
% neurons fired. It is the amount to change weights by.
dApre = zeros(N, num_connections);
dApost = zeros(N, 1);
active_spikes = cell(delay_max, 1);  % To track when spikes arrive
active_idx = 1;
% TODO - later this is where I would put the pre variable

% Info logging variables
%spike_arrival_trace = [];
spike_times_trace = [];
vt = zeros(N, ms_per_sec);
vt(:, 1) = v;
debug = [];

%% DATA
%inp = [50, 100, 200, 500, 700, 800, 801, 802, 803, 804];
%ts = [5, 7, 50, 55, 100, 525, 505, 500, 500, 500];
inp = [800 1003 801 1003];
ts = [500 520 530 550];
[ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp );

%% Main computation loop
for sec = 1 : sim_time_sec
    [ inp, ts, patt_inp, patt_ts ] = embedPat( N_inp, patt_inp, patt_ts );
    ts = ts + (sec-1) * 1000;
    for ms = 1 : ms_per_sec
        time = (sec - 1) * ms_per_sec + ms;  
        
        %% Calculate input to neurons at this time step
        Iapp = zeros(N, 1);
        t0 = time - last_spike_time;
        t0negu = t0 - delays;
        scale = 1 ./ (sigma * sqrt(2 * pi));
        g = scale .* exp((-1/2) .* ((t0negu) ./ sigma) .^2 );
        g(isnan(g)) = 0;
        gaussian_values = w .* g;
      
        % Collect input for each neuron based on synapses facing them
        [to_neurons, ~, to_neurons_idx] = unique(post);
        Iapp(to_neurons) = accumarray(to_neurons_idx, gaussian_values(:));
        
        debug = [debug; mean(w(800, :)), mean(w(801, 1)), mean(w(802, 1))];
        if time == 2110
           disp(''); 
        end
        
        %% Update STDP
        incoming = active_spikes{active_idx}; % Get spikes arriving now
        if ~isempty(incoming)
            active_spikes{active_idx} = [];
            incoming_idxs = sub2ind(size(dApre), incoming(:, 1), incoming(:,2));

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
            
            presynaptic_idxs = find(post == neuron);   
            % I spiked, reinforce any before me
            w(presynaptic_idxs) = w(presynaptic_idxs) + dApre(presynaptic_idxs);
            % weaken connections to any after me who spiked not long ago
            w(neuron, :) = w(neuron, :) + dApost(post(neuron, :))'; 
            
            dApost(neuron) = dApost(neuron) + Apost;

            for connection = 1:num_connections
                % Keep track of when this spike will arrive at each
                % postsynaptic neuron and through which connection.
                delay = delays(connection);
                arrival_offset = mod(active_idx + delay - 1, delay_max) + 1;
                active_spikes{arrival_offset} = [active_spikes{arrival_offset}; neuron, connection];
            end
        end
        
        %% Update STDP variables
        dApre = dApre * STDPdecaypre;
        dApost = dApost * STDPdecaypost;
        active_idx = mod(active_idx, delay_max) + 1;
        
        %% Apply synaptic scaling (SS) and weight bounding
        w_mean = mean(w(:));
        % TODO - Only bounds against depression, should also bound above
        if w_mean < syn_mean_thresh
            w = w .* (syn_mean_thresh / w_mean);
        end
        
        % Limit w to between [0, w_max]
        w(1:N_inp, :) = max(0, min(w_max, w(1:N_inp,:))); 
        
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
    plot(l1_times - offset, l1_idxs, '.', 'MarkerSize', 8);
    hold on 
    plot(l2_spike_times - offset, l2_spike_idxs, '.r', 'MarkerSize', 14)
    if numel(l2_spike_times) > 0
        for i = l2_spike_times
            plot( [i i] - offset, get( gca, 'Ylim' ), '--r', 'LineWidth',2)
        end
    end
    
    
    subplot(2, 1, 2);
    plot(1:ms_per_sec, vt(N, :));
    title(sprintf('second: %d', sec-1));
    %axis([500 600 -Inf Inf]);
    legend({'Output neuron'});
    drawnow;
    
    hold off;
    
    %% Reset second facilitating variables
    spike_times_trace = [];
    spike_arrival_trace = [];
    vt = zeros(N, ms_per_sec);
    vt(:, 1) = v;
    
    fprintf('Second: %d\n', sec);
end