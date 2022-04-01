%   KVSKYLINE skyline projection 2D contour plot 
%   kvskyline(figure_number, ax, y, varargin)
%   'ax' parameters
%      ax.contour, ax.tc, ax.filt - see KAZAN viewer 
%   additional parameters (first value is default):
%      position                   - the position of the plot on the figure 
%      skyline   <0.1>            - width of skyline subplot
%      multyplot   <0>            - multiple plots per figure
%      skyspace  <0.1>            - space between skyline projections and axis
%      plot    <'image'>/'contour'  - plot type
%      hold    'on'/<'off'>       - hold previous plot  
%      style   <0>/1              - predefine styles
function hdl = kvskyline(n, ax, y, varargin)

if nargin > 3,
    if ~mod(nargin-3,2)
        for kk=1:2:nargin-3
            ax=setfield(ax, lower(varargin{kk}), varargin{kk+1});    
        end
    else
        error('Wrong amount of arguments')
    end
end
brd   = safeget(ax, 'skyline', .1);
ptype = safeget(ax, 'plot', 'image');
holdmode = lower(safeget(ax, 'hold', 'off'));

skylinex = safeget(ax,'skylinex',max(real(y), [], 1));
skyliney = safeget(ax,'skyliney',max(real(y), [], 2));
contr    = safeget(ax, 'contour', [0.1,0.2]);
skyspace = safeget(ax, 'skyspace', 0.2);
bounds   = safeget(ax, 'position', [0,0,1,1]);
multyplot = safeget(ax, 'multyplot', 0);
style = safeget(ax, 'style', 0);

figure(n)
hdl = get(n, 'UserData');
if isempty(hdl) | length(hdl) ~= 3,  holdmode = 'off'; end

AxLim = safeget(ax, 'Lim', [min(ax.x), max(ax.x), min(ax.y), max(ax.y)]);

% up plot
if multyplot,
    hdl(1) = axes('position',[.1*bounds(3)+bounds(1), (.9-brd)*bounds(4)+bounds(2),  ...
            (.8-brd)*bounds(3), brd*bounds(4)]);
elseif strcmp(lower(holdmode), 'off')
    clf;
    hdl(1) = axes('position',[.1*bounds(3)+bounds(1), (.9-brd)*bounds(4)+bounds(2),  ...
            (.8-brd)*bounds(3), brd*bounds(4)]);
else
    set(n, 'CurrentAxes', hdl(1));
    hold('on')
end

plot(ax.x, skyliney, 'Color', safeget(ax,'Color','b'));
set(hdl(1), 'XAxisLocation', 'top');
set(hdl(1), 'YTick', []);
set(hdl(1), 'XLimMode', 'manual', 'XLim', AxLim(1:2));
min1 = min(skyliney); max1 = max(skyliney);
set(hdl(1), 'YLimMode', 'manual', 'YLim', [min1,max1]+skyspace*(max1-min1)*[-1,1]/2);

% right plot 
if strcmp(lower(holdmode), 'off')
    hdl(3) = axes('position', [(.9-brd)*bounds(3)+bounds(1), .1*bounds(4)+bounds(2),...
            brd*bounds(3), (.8-brd)*bounds(4)]);
else
    set(n, 'CurrentAxes', hdl(3));
    hold('on')
end
plot(skylinex,ax.y, 'Color', safeget(ax,'Color','b'));
set(hdl(3), 'YAxisLocation', 'right');
set(hdl(3), 'XTick', []);
set(hdl(3), 'YLimMode', 'manual', 'YLim', AxLim(3:4));
set(hdl(3), 'XLimMode', 'manual', 'XLim', [min1,max1]+skyspace*(max1-min1)*[-1,1]/2);

% center plot
y = processing(y, ax);

maxyy=max(max(real(y))); minyy=min(min(real(y)));
if maxyy < 0, sh = minyy; else sh = 0; end
clim = real(contr*maxyy+sh);

if strcmp(lower(holdmode), 'off')
    hdl(2) = axes('position', [.1*bounds(3)+bounds(1), .1*bounds(4)+bounds(2),...
            (.8-brd)*bounds(3), (.8-brd)*bounds(4)]);
else
    set(n, 'CurrentAxes', hdl(2));
    hold('on')
end

if strcmp(lower(ptype), 'image')
    h = imagesc(ax.x, ax.y, real(y).', [min(clim), max(clim)]);
    set(gca, 'YDir', 'normal')
else
    contour(ax.x, ax.y, real(y+sh).', clim);
end
set(hdl(2), 'XLimMode', 'manual', 'XLim', AxLim(1:2), ...
    'YLimMode', 'manual', 'YLim', AxLim(3:4));
xlabel(safeget(ax, 'xlabel', ''),'UserData', 2)
ylabel(safeget(ax, 'ylabel', ''),'UserData', 3)

set(n, 'UserData', hdl);
set(n, 'CurrentAxes', hdl(2));

for ii=1:3
    set(hdl(ii),'UserData',ii);
end

switch style
    case 1,
        set(hdl(1),'Visible','off');
        set(hdl(3),'Visible','off');
        set(hdl(2),'Box','off','TickDir','out');
        if strcmp(lower(holdmode), 'off')
            axes('position', [.1*bounds(3)+bounds(1), .1*bounds(4)+bounds(2),...
                    (.8-brd)*bounds(3), (.8-brd)*bounds(4)],...
                    'XTick',[],'YTick',[],'Box','on','Color','none',...
                    'UserData',4);
        end
end