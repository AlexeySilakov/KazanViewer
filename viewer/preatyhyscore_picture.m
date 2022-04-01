function hh = preatyhyscore_picture(x, y, z, levs, range, orientation, varargin)
%preatyhyscore_picture(x, y, z, levs, maxFreq, orientation)
% x, y - axes in MHz
% z - 2D matrix, size(z) = [numel(x), numel(y)]
% levs = [minLev, MaxLev]
% range - in MHz (if two numbers, one is lower limit other one is 
% orientation: 0 - (++) quadrant only
%              1 - (++) (+-), vertical
%              2 - (++) (-+), horisontal
%              3 - all quadrants

if length(range)==2
    minFreq = range(1);
    maxFreq = range(2);
    orientation = 0;
else
    minFreq = 0;
    maxFreq = range;
end

df = (maxFreq -minFreq)./[3, 4, 5, 6, 7]; 
idx = find( [df*10 - floor(df*10)]==0 );
if ~isempty(idx)
    df = df(max(idx));
else
    df = df(2);
end
switch orientation
    case 0 
        XLim = [minFreq, maxFreq];
        YLim = [minFreq, maxFreq];
        XTick = [minFreq:df:maxFreq];
        YTick = [minFreq:df:maxFreq];
    case 1
        XLim = [0, maxFreq];
        YLim = [-maxFreq, maxFreq];
        XTick = [0:df:maxFreq];
        YTick = [-maxFreq:df:maxFreq];
    case 2
        XLim = [-maxFreq, maxFreq];
        YLim = [0, maxFreq];
        XTick = [-maxFreq:df:maxFreq];
        YTick = [0:df:maxFreq];
    case 3
        XLim = [-maxFreq, maxFreq];
        YLim = [-maxFreq, maxFreq];
        XTick = [-maxFreq:df:maxFreq];
        YTick = [-maxFreq:df:maxFreq];
end

A = [1.0000    1.0000    1.0000; ...
    0         0    1.0000; ...
    0    1.0000    1.0000; ...
    1.0000    1.0000         0; ...
    1.0000         0         0; ...
    0.5000         0         0];

xx = 1:6;
yy = 1:0.1:6;
BB = interp1(xx, A, yy, 'linear');


hh(1) = imagesc(x, y, z, levs); hold on;
plot([0, maxFreq], [0, maxFreq], ':k');
if orientation==1, 
    plot([0, maxFreq], [0, 0], 'k');
    plot([0, maxFreq], [0, -maxFreq], ':k');
end
if orientation==2, 
    plot([0, 0], [0, maxFreq], ':k');
    plot([0, -maxFreq], [0, maxFreq], ':k');
end

axis equal; box on; grid on;
set(gca, 'XLim',XLim, 'YLim', YLim, ...
    'XTick', XTick, 'YTick', YTick, ...
    'TickDir', 'out', 'TickLength', [0.02, 0.01], ...
    'FontSize', 15, 'YDir', 'normal');
colormap(gca, BB);
if ~isempty(varargin) 
    if ~strcmp(varargin{1}, 'no labels')
        xlabel('Frequency, MHz', 'FontSize', 18);
        ylabel('Frequency, MHz', 'FontSize', 18);
    end
else
        xlabel('Frequency, MHz', 'FontSize', 18);
        ylabel('Frequency, MHz', 'FontSize', 18);
    
end
hold off;

