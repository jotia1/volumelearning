function [ output ] = runsinglelayer( net )
%% RUNSINGLELAYER - Takes a network data struct and runs network
%   Given a valid network description as defined by validatenetwork
%   description, run that network and calculate the output.

output = struct();
assert(validatenetwork(net), 'Error with network description');

N = net.N_inp + net.N_hid + net.N_out;
layer_sizes = [N_inp, N_hid, N_out];

ms_per_sec = 1e3;

% Output data
output.spike_times_trace = [];


for sec = 1 : net.sim_time_sec
    
    for ms = 1 : ms_per_sec
        
        
    end
    
    
end

end