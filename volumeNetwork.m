%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 

% Features
%   - 
%


% Often tweaked parameters
N_inp = 1000;
N_out = 5;
N = N_inp + N_out;
sim_time_sec = 3;
delay_max = 10;
num_connections = 3;
w_init = 2;
w_max = 10;
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

% Spike propagation 
active_spikes = cell(delay_max, 1);
active_idx = 1;

% Synpatic links 
% Post is who is postsynaptic from a pre. N x M, M is number of connections
post = ones(N, num_connections) * N; 
post(700, 1) = N-1; %TODO - Create actual connectivity, this is just a demo
post(800, 1) = N-2;

% Synapse dynamics parameters
g = zeros(N, 1);
sigma_max = 30;
sigma_min = 2;  % TODO where did i get all these variables from...
sigma = rand(N, 1) * (sigma_max - sigma_min) + sigma_min;
delays = rand(N, num_connections) * delay_max;
last_spike_time = zeros(N, 1) * -Inf;

% Info logging variables
spike_arrival_trace = [];
vt = zeros(N, ms_per_sec);
vt(:, 1) = v;
debug = [];

%% DATA
inp = [50, 100, 200, 500, 700, 800];
ts = [2005, 2007, 2050, 2055, 2100, 2500];

%% Main computation loop
for sec = 1 : sim_time_sec
    for ms = 1 : ms_per_sec
        time = (sec - 1) * ms_per_sec + ms;  
        
        %% Calculate input to neurons at this time step
        Iapp = zeros(N, 1);
        t0 = time - last_spike_time;
        t0negu = t0 - u;
        scale = 1 ./ (sigma * sqrt(2 * pi));
        g = scale .* exp((-1/2) .* ((t0negu) ./ sigma) .^2 );
        g(isnan(g)) = 0;
        gaussian_values = w .* g;
      
        % Collect input for each neuron based on synapses facing them
        [to_neurons, ~, to_neurons_idx] = unique(post);
        Iapp(to_neurons) = accumarray(to_neurons_idx, gaussian_values(:));
        
        debug = [debug; sum(Iapp(:))];
        if time == 2110
           disp(''); 
        end
        
        %% Deal with spikes that have arrived
%         incoming = active_spikes{active_idx}; % Get spikes arriving now
%         active_spikes{active_idx} = [];
%         
%         for spike = 1 : length(incoming)
%             from_neuron = incoming(spike);
%             to_neuron = post(from_neuron);
%             % TODO - remove below, instantaneous application
%             %v(to_neuron) = v(to_neuron) + w(from_neuron); % TODO - conductance based?
%         end
%         
%         Log spike arrivals for plotting later
%         spike_arrival_trace = [spike_arrival_trace; 
%                                time*ones(length(incoming),1), incoming];

        
        %% Update membrane voltages
        v = v + (v_rest + Iapp - v) / neuron_tau;
        vt(:, ms) = v;
        
        %% Deal with neurons that just spiked
        fired = [find(v >=v_thres); inp(ts == time)'];
        last_spike_time(fired) = time;  % Used for calculating dynamics later
        
        for spike = 1 : length(fired)
            from_neuron = fired(spike);
            v(from_neuron) = v_reset;
        end
        
        
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
    plot(1:ms_per_sec, vt(N-2: end, :));
    title(sprintf('second: %d', sec-1));
    legend({'N-2','N-1','N'});
    drawnow;
    
    %% Reset second facilitating variables
    spike_times_trace = [];
    spike_arrival_trace = [];
    vt = zeros(N, ms_per_sec);
    vt(:, 1) = v;
    
    fprintf('Second: %d\n', sec);
end