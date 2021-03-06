function varargout = spinbox(varargin)
% GUI object constructor
% AlSi 16.01.05
% handles of the objects are stored in push1;
% structure of the parameters are in push2
if nargin==0
   handleFig=gcf;
   action = 'new';
elseif nargin>=1 & ishandle(varargin{1}) & strcmp(get(varargin{1},'type'),'figure')
    action = 'new';
   handleFig=varargin{1};
elseif isa(varargin{1}, 'spinbox')
    hands = varargin{1};   
    action = 'set';
else
    return;
end
par = [];
% contruct structure of the new parameters
if nargin>1
    if mod(nargin-1, 2), 
        error('Number of input parameters must be ODD: set(handle, PARAM, VALUE)!!');
    end
    for ci = 2:2:nargin-1
        par = setfield(par, varargin{ci}, varargin{ci+1});
    end
end

switch action
    case 'new'
        prop = Local_getdeff;
               
        ed.edit = uicontrol(handleFig, 'Style', 'edit');
        ed.push1 = uicontrol(handleFig, 'Style', 'push');
        ed.push2 = uicontrol(handleFig, 'Style', 'push');
        
        builtin('set', ed.edit, 'buttondownfcn', ...
             ['spinboxcb(get(', sprintf('%20.15f',ed.push1),  ', ''userdata''), ''edit'')']);
         builtin('set', ed.edit, 'callback', ...
             ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''edit'')']);
%         builtin('set', ed.push1, 'buttondownfcn', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push1'')']);
         builtin('set', ed.push1, 'callback', ...
             ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push1'')']);
        
%         builtin('set', ed.push2, 'buttondownfcn', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push2'')']);
         builtin('set', ed.push2, 'callback', ...
             ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push2'')']);

        Local_resize(ed, prop);
        ed=class(ed,'spinbox');
        
        builtin('set',ed.push1,'userdata',ed);
        builtin('set',ed.push2,'userdata',prop);
        prop = Local_setparam(ed, prop);
        prop = Local_setparam(ed, par);
%        set(ed, 'Visible', 'on');
    case 'set'
        ed = hands;
        Local_setparam(ed, par);
end
if nargout == 1
    varargout{1} = ed;
end

% --------------------------------------------------
function Local_resize(ed ,prop)
buttonWidth = prop.Position(4);
switch lower(prop.Direction)
    case 'vertical'
        buttonWidth = prop.Position(4);
        pb_pos1 = [0, prop.Position(2), buttonWidth, 1/2*prop.Position(4)];
        pb_pos2 = [0, prop.Position(2) + prop.Position(4)/2, buttonWidth, 1/2*prop.Position(4)];
        multWidth = 1;
        CDATA_rot = 0;        
    case 'horizontal'
        buttonWidth = prop.Position(4)/2;
        pb_pos1 = [0, prop.Position(2), buttonWidth, prop.Position(4)];
        pb_pos2 = [0, prop.Position(2), buttonWidth, prop.Position(4)];
        CDATA_rot = 1; 
        multWidth = 2;        
end

switch prop.ButtonStyle
    case 'left'
        pb_pos1(1) = prop.Position(1);
        pb_pos2(1) = prop.Position(1) + buttonWidth*(multWidth-1);
        ed_posx = prop.Position(1)+buttonWidth*multWidth;
    case 'right'
        pb_pos1(1) = sum(prop.Position([1, 3])) - buttonWidth*multWidth;
        pb_pos2(1) = sum(prop.Position([1, 3])) - buttonWidth;
        ed_posx = prop.Position(1);
end
set(ed.edit, 'Position', [ed_posx, prop.Position(2), ...
        prop.Position(3)-buttonWidth*multWidth, prop.Position(4)]);
set(ed.push1, 'Position', pb_pos1);
set(ed.push2, 'Position', pb_pos2);

% --------------------------------------------------
function out = Local_getdeff
out.BackgroundColor = get(0, 'defaultuicontrolbackgroundcolor');
out.Callback = '';
out.Enable = 'on';
out.FontAngle = get(0, 'defaultuicontrolFontAngle');
out.FontName = get(0, 'defaultuicontrolFontName');
out.FontSize = get(0, 'defaultuicontrolFontSize');
out.FontUnits = get(0, 'defaultuicontrolFontUnits');
out.FontWeight = get(0, 'defaultuicontrolFontWeight');
out.ForegroundColor = get(0, 'defaultuicontrolForegroundColor');
out.HorizontalAlignment = get(0, 'defaultuicontrolHorizontalAlignment');
out.Max = 0;
out.Min = 0;
out.ButtonStyle = 'right';
out.Direction = 'vertical';
out.Position = [10, 10, 50, 21];
out.Step = 0.1;
out.TooltipString = '';
out.Units = get(0, 'defaultuicontrolunits');
out.Value = 0;
out.Tag = '';
out.ButtonDownFcn = '';
out.Clipping = 'on';
out.CreateFcn = '';
out.DeleteFcn = '';
out.BusyAction = 'queue';
out.HandleVisibility = 'on';
out.HitTest = 'on';
out.Interruptible = get(0, 'defaultuicontrolInterruptible');
out.SelectionHighlight = get(0, 'defaultuicontrolSelectionHighlight');
out.UIContextMenu = [];
out.UserData = [];
out.Visible = 'on';
out.EnableButton = 'on'
% ---------------------------------------------------
function varargout = Local_setparam(ed, newpar)

par = get(ed);
if isempty(newpar), 
    if nargout ==1, varargout{1} = par; end
    return; 
end
fnames = fieldnames(newpar);

% common parameters
commpar = {'BackgroundColor', 'ForegroundColor', ...
        'TooltipString', 'Units', 'Tag', 'Clipping', 'CreateFcn', ...
        'DeleteFcn', 'BusyAction', 'HandleVisibility', 'HitTest', ...
        'Interruptible', 'SelectionHighlight', 'UIContextMenu', ...
        'Visible', 'Enable'};

% separate 
seppar = {'Callback', 'Max', 'Min', 'Step', 'Position', ...
        'Site', 'Direction', 'Value', 'ButtonDownFcn',...
        'UserData'};
sep_method = [0, 0, 0, 0, 1,...
        1, 1, 2, 0, ...
        0]; % 0==do nothing, 1==Local_resize, 2==Local_value
% for edit
editpar = {'HorizontalAlignment', 'FontAngle', 'FontName', ...
        'FontSize', 'FontUnits', 'FontWeight'};
UpDownpar = {'EnableButtons'};
make_resize = 0;
for ci = 1:length(fnames)
    try
        newvalue = getfield(newpar, fnames{ci});
        par = setfield(par, fnames{ci}, newvalue);
        if sum(strcmp(commpar, fnames{ci}))
            builtin('set', ed.edit, fnames{ci}, newvalue);
            builtin('set', ed.push1, fnames{ci}, newvalue);
            builtin('set', ed.push2, fnames{ci}, newvalue);
        elseif sum(strcmp(seppar, fnames{ci}))
            num = find(strcmp(seppar, fnames{ci}));
            switch sep_method(num)
                case 1 % resize
                    make_resize = 1;
                case 2
                    builtin('set', ed.edit, 'String', newvalue);
            end
        elseif sum(strcmp(editpar, fnames{ci}))
            builtin('set', ed.edit, fnames{ci}, newvalue);
        elseif strcmp('EnableButtons', fnames{ci})
            builtin('set', ed.push1, 'Enable', newvalue);
            builtin('set', ed.push2, 'Enable', newvalue);
            
            if strcmp(newvalue, 'on')
                out = Local_ico;
                out(:, :, 1) = out(:, :, 1)*par.BackgroundColor(1);
                out(:, :, 2) = out(:, :, 2)*par.BackgroundColor(2);
                out(:, :, 3) = out(:, :, 3)*par.BackgroundColor(3);
            else
                out(1, 1, 1) = par.BackgroundColor(1);
                out(1, 1, 2) = par.BackgroundColor(2);
                out(1, 1, 3) = par.BackgroundColor(3);
            end
                builtin('set',ed.push1,'CDATA', out(end:-1:1, :, :));
                builtin('set',ed.push2,'CDATA', out);
        end
        
    catch
        error(['Wrong parameter name: ', fnames{ci}, '.'])
    end
end

set(ed.push2, 'Userdata', par);
if make_resize
    Local_resize(ed, par);
end
if nargout ==1
    varargout{1} = par;
end

function out = Local_ico
out(:, :, 1) = [1   1   1   1   1   1   1   1   1   1;
                1   1   1   1   1   1   1   1   1   1;
                1   1   1   1   0   0   1   1   1   1;
                1   1   1   0   0   0   0   1   1   1;
                1   1   0   0   0   0   0   0   1   1;
                1   0   0   0   0   0   0   0   0   1;
                0   0   0   0   0   0   0   0   0   0;
                0   0   0   0   0   0   0   0   0   0;
                1   1   1   1   1   1   1   1   1   1;
                1   1   1   1   1   1   1   1   1   1];
out(:, :, 2) = out(:, :, 1);
out(:, :, 3) = out(:, :, 1);