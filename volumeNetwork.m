%% VOLUMENETWORK - A network for learning spatio-temporal volumes
% Based around the idea that temporal information plays a critical role in
% learning real world signals this network models conduction delays between
% neurons and will form synapses essentially sampling for useful places. 

% Features
%   - 
%


% Often tweaked parameters
w_init = 0.5;
N_inp = 1000;
N_out = 5;
N = N_inp + N_out;
sim_time_sec = 100;
delay_max = 20;

% Neuron parameters
v_rest = -65;
v_reset = -70;
v_thres = -55;
v = v_rest * ones(N, 1);
w = ones(N, 1) * w_init;

% Facilitating variables
activespikes = cell(delay_max, 1);
post = ones(N, 1) * N;


%% Main computation loop
for sec = 1 : sim_time_sec
    for ms = 1 : 1000
        time = (sec - 1) * 1000 + ms;  
        
        %% Deal with spikes that have arrived
        incoming = activespikes{active_idx}; % Get spikes arriving now
        activespikes{active_idx} = [];
        
        for spike = 1 : length(incoming)
            from_neuron = incoming(spike);
            
            w(from_neuron)
        end
        
        %% Update membrane voltages
        v = v + (v_rest + 
        
        %% Deal with neurons that just spiked
        
        
        
        
        
    end
end