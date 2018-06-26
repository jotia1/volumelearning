function [xs, ys, ts, ps] = loadDVSsegment(size, do_plot, eidx)
%LOADDVSSEGMENT Load and return trimmed DVS data
%   size - Size of region to sample (32, 16, 8)
%   do_plot - whether or not to plot
%   eidx - End index (max number of events to keep) 

% [xs, ys, ts, ps] = loadDVSsegment(16, true, 500)

addpath('../AedatTools/Matlab/');

aedat = struct;
aedat.importParams.filePath = '../data/lots_dots.aedat';

aedat = ImportAedat(aedat);
%PlotAedat(aedat,5, 'time')


all_data = aedat.data.polarity;

ts = all_data.timeStamp;
ys = all_data.y;
xs = all_data.x;
ps = all_data.polarity;
num = all_data.numEvents;

%plot3(xs, ys, ts, 'k.');

nts = [];
nys = [];
nxs = [];
nps = [];


% % 32x32 dots crossing
if size == 128
    xloc = floor(((190 - 128) / 2) + 50);
    yloc = floor(((180 - 128) / 2));
    nsize = 128;
elseif size == 32
    xloc = 105;
    yloc = 95;
    nsize = 32;
elseif size == 16
    % 16x16 dots crossing
    xloc = 112;
    yloc = 104;
    nsize = 16;
elseif size == 8
    % 8x8 dots crossing
    xloc = 117;
    yloc = 109;
    nsize = 8;
else 
    disp('NOT A VALID SIZE, TODO, cause ERROR'); 
end

idx = 1;
while idx < num
    
%     if xs(idx) <  nsize
%        disp(); 
%     end
    
    % if event is in the 32x32 rectangle at position (0, 104)
%     if xs(idx) > 240 - nsize - 1 && ...
%             ys(idx) > 104 && ys(idx) < 104 + nsize + 1
    
    % For collecting crossing dots
    if xs(idx) >= xloc && xs(idx) < xloc + nsize && ...
            ys(idx) >= yloc && ys(idx) < yloc + nsize
        % keep it
        nts(end + 1) = ts(idx);
        nys(end + 1) = ys(idx);
        nxs(end + 1) = xs(idx);
        nps(end + 1) = ps(idx);
        %disp('hi');
        
    end
    
    idx = idx + 1;
end

aedat.data.polarity.timeStamp = nts;
aedat.data.polarity.y = nys;
aedat.data.polarity.x = nxs;
aedat.data.polarity.polarity = nps;


% final cleaning touches
%eidx = 500;
if eidx == 0 || eidx > numel(nxs)
    eidx = numel(nxs);
end
xs = nxs(1:eidx) - (xloc - 1);   % Lets shift to 1 indexing now
ys = nys(1:eidx) - (yloc - 1);               % Same here
ts = nts(1:eidx) - nts(1);

%% Clean hot pixels
if nsize == 128
    idxs = [find(xs == 59 & ys == 104), ...
        find(xs == 101 & ys == 62), ...
        find(xs == 104 & ys == 5), ...
        find(xs == 34 & ys == 77), ...
        find(xs == 60 & ys == 18), ...
        find(xs == 2 & ys == 6), ...
        find(xs == 67 & ys == 50), ...
        find(xs == 44 & ys == 74)];
    xs(idxs) = [];
    ys(idxs) = [];
    ts(idxs) = [];
    
elseif nsize == 32
    idxs = [find(xs == 20 & ys == 5), find(xs == 10 & ys == 8)];
    xs(idxs) = [];
    ys(idxs) = [];
    ts(idxs) = [];
end

if do_plot
    subplot(2, 1, 1);
    plot(ts/1000000, sub2ind([nsize, nsize], xs, ys), '.k'); 
    axis([0 Inf 0 Inf]);

    subplot(2, 1, 2);
    plot3(xs, ys, ts/1000000, 'k.');
    axis([0 nsize 0 nsize -Inf Inf]);
    xlabel('xs'); ylabel('ys'); zlabel('ts (sec)');
    title('32x32 segment of lots dots datafile');

    figure
    pts = linspace(0, nsize, nsize);
    N = histcounts2(ys(:), xs(:), 0:nsize+1, 0:nsize+1 );
    imagesc(pts, pts, N);
    axis equal;
    set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
    colorbar
end

end

