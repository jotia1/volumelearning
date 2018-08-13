function [ inp, ts ] = gen1Ddata( length_sec, num_pixels, ms_per_pixel, spatial_noise, temporal_noise_ms )
%% GEN1DDATA - Generate data of a dot moving along a 1D retina

ms_to_sec = 1000;

num_spikes = floor((length_sec * ms_to_sec) / ms_per_pixel);

%TODO: maybe make this gaussian?, would need to decide on a mean and 
% variance then though..
spatial_offsets = randi([-spatial_noise, spatial_noise], 1, num_spikes);  
inp = mod((1:num_spikes) + spatial_offsets, num_pixels) + 1;

temporal_offsets = randi([-temporal_noise_ms, temporal_noise_ms], 1, num_spikes);
ts = (1:num_spikes) * ms_per_pixel + temporal_offsets;

end