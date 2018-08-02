function [ nxs, nys, nts, nps ] = dvs2patch( xs, ys, ts, ps, patch_size, length_seconds, break_size_ms )
%% DVS2PATCH - Takes events and randomly rescales to something patch_size
%   Takes random regions (of patch_size by patch_size) and creates new
%   events that spatially fit into patch_size but for longer periods of
%   time.
%   It is designed to make one spatially smaller but temporally longer
%   sequence based on some original data. E.g. started with something
%   128*128*10secs make a sequence 8*8*500secs. Where the second sequence
%   is random regions of the first sequence seperated by break_size_ms of
%   quiet.
%
%   Example usage:
%       % Assumes xs, ys, ts, ps is preloaded.
%       [ nxs, nys, nts, nps ] = dvs2patch( xs, ys, ts, ps, 8, 100, 50 )
%
%   Assumptions;
%       - Assumes input is 128x128
%       - Assumes length_seconds > ts(end) - ts(1)


uS2SEC = 1e6;
MS2uS = 1e3;
INP_RESOLUTION = [128, 128];

rxmax = INP_RESOLUTION(1) - patch_size; rymax = INP_RESOLUTION(2) - patch_size;
rxmin = 1; rymin = 60;

inp_length_sec = (ts(end) - ts(1)) / uS2SEC;  
%assert(inp_length_sec < length_seconds, 'Input is temporally too long for given max length');

n_regions = ceil(length_seconds / inp_length_sec);
regions = [randi([rxmin, rxmax], n_regions, 1), randi([rymin, rymax], n_regions, 1)];

nxs = [];
nys = [];
nts = [];
nps = [];
last_time = 0;

for i = 1 : n_regions
    rx = regions(i, 1);
    ry = regions(i, 2);
    
    idxs = find(xs > rx & xs <= rx + patch_size ...
        & ys > ry & ys <= ry + patch_size);
    
    nxs = [nxs, xs(idxs) - rx];
    nys = [nys, ys(idxs) - ry];
    nts = [nts, ts(idxs) + last_time + (break_size_ms * MS2uS)];
    last_time = nts(end);
    
end

trim_idxs = find(nts > length_seconds * uS2SEC);
nxs(trim_idxs) = [];
nys(trim_idxs) = [];
nts(trim_idxs) = [];

end