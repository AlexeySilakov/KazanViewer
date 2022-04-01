function varargout = spinbox(varargin)
% GUI object constructor
% AlSi 16.01.05
% handles of the objects are stored in push1;
% structure of the parameters are in push2
if nargin==0
   handleFig=gcf;
   action = 'new';
elseif isa(varargin{1}, 'char') %%% for ML>2014, callbacks need to be handeled diffently... uhhhhhh
%     builtin('gcbo')
%     builtin('gco')
    spinboxcb;
    action = '';
elseif nargin>=1 && ishandle(varargin{1}) && (strcmp(get(varargin{1},'type'),'figure')|| strcmp(get(varargin{1},'type'),'uipanel'))
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
    if isstruct(varargin{2})
        par = varargin{2};
    else
        if mod(nargin-1, 2),
            error('Number of input parameters must be ODD: set(handle, PARAM, VALUE)!!');
        end
        for ci = 2:2:nargin-1
            par = setfield(par, varargin{ci}, varargin{ci+1});
        end
    end
end

switch action
    case 'new'
        
        prop = Local_getdeff;
               
        ed.edit = uicontrol(handleFig, 'Style', 'edit');
        ed.push1 = uicontrol(handleFig, 'Style', 'push');
        ed.push2 = uicontrol(handleFig, 'Style', 'push');
        
%          builtin('set', ed.edit, 'buttondownfcn', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1),  ', ''userdata''), ''edit'')']);
%          
%          builtin('set', ed.edit, 'callback', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''edit'')']);
%          builtin('set', ed.push1, 'callback', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push1'')']);
%          builtin('set', ed.push2, 'callback', ...
%              ['spinboxcb(get(', sprintf('%20.15f',ed.push1), ',''userdata''), ''push2'')']);
% 
%          builtin('set', ed.edit, 'KeyPressFcn', ...
%               ['spinboxcb(get(', sprintf('%20.15f',ed.push1),  ', ''userdata''), ''keypress'')']);
%          builtin('set', ed.push1, 'KeyPressFcn', ...
%               ['spinboxcb(get(', sprintf('%20.15f',ed.push1),  ', ''userdata''), ''keypress'')']);
%          builtin('set', ed.push2, 'KeyPressFcn', ...
%               ['spinboxcb(get(', sprintf('%20.15f',ed.push1),  ', ''userdata''), ''keypress'')']);


    
        Local_resize(ed, prop);
        ed=class(ed,'spinbox');

        builtin('set', ed.edit, 'callback',      @(hObject,callbackdata)spinboxcb(ed, 'edit') );
        builtin('set', ed.edit, 'buttondownfcn', @(hObject,callbackdata)spinboxcb(ed, 'edit') );
        builtin('set', ed.push1, 'callback',     @(hObject,callbackdata)spinboxcb(ed, 'push1') );
        builtin('set', ed.push2, 'callback',     @(hObject,callbackdata)spinboxcb(ed, 'push2') );

        builtin('set', ed.edit, 'KeyPressFcn',   @(hObject,callbackdata)spinboxcb(ed, 'keypress') );
        builtin('set', ed.push1, 'KeyPressFcn',  @(hObject,callbackdata)spinboxcb(ed, 'keypress') );
        builtin('set', ed.push2, 'KeyPressFcn',  @(hObject,callbackdata)spinboxcb(ed, 'keypress') );

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
switch lower(prop.Direction(1))
    case 'v'
        buttonWidth = prop.Position(4);
        pb_pos1 = [0, prop.Position(2), buttonWidth, 1/2*prop.Position(4)];
        pb_pos2 = [0, prop.Position(2) + prop.Position(4)/2, buttonWidth, 1/2*prop.Position(4)];
        multWidth = 1;
        CDATA_rot = 0;        
    case 'h'
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
out.Precision = 8;
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
out.EditBackground = [1, 1, 1];
out.UserData = [];
out.Visible = 'on';
out.EnableButton = 'on';

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
        'TooltipString', 'Units', 'Tag', 'CreateFcn', ...
        'DeleteFcn', 'BusyAction', 'HandleVisibility', 'HitTest', ...
        'Interruptible', 'SelectionHighlight', 'UIContextMenu', ...
        'Visible', 'Enable'}; %'Clipping', 

% separate 
seppar = {'Callback', 'Max', 'Min', 'Step', 'Position', ...
        'Site', 'Direction', 'Value', 'ButtonDownFcn',...
        'String'};
sep_method = [0, 0, 0, 0, 1,...
        1, 1, 2, 3, ...
        4]; % 0==do nothing, 1==Local_resize, 2==Local_value
% for edit
editpar = {'HorizontalAlignment', 'FontAngle', 'FontName', ...
        'FontSize', 'FontUnits', 'FontWeight', 'EditBackground', 'UserData'};
UpDownpar = {'EnableButton'};
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
                    
                    builtin('set', ed.edit, 'String', sprintf('%f', newvalue));
%                     par.Value = newvalue;
                case 3 % callbacks
                    builtin('set', ed.edit, fnames{ci}, newvalue);
                    builtin('set', ed.push1, fnames{ci}, newvalue);
                    builtin('set', ed.push2, fnames{ci}, newvalue);   
                case 4
                    builtin('set', ed.edit, 'String', newvalue);
                    par.Value = str2num(newvalue);
            end
        elseif sum(strcmp(editpar, fnames{ci}))
            if strcmp(fnames{ci}, 'EditBackground')
               builtin('set', ed.edit, 'BackgroundColor', newvalue);
            else
               builtin('set', ed.edit, fnames{ci}, newvalue);
            end
        elseif strcmp('EnableButton', fnames{ci})
            builtin('set', ed.push1, 'Enable', newvalue);
            builtin('set', ed.push2, 'Enable', newvalue);
            Local_redraw(ed, par);
            
        end
        
    catch
        error(['Wrong parameter name: ', fnames{ci}, '.'])
    end
end

set(ed.push2, 'Userdata', par);
if make_resize
    Local_resize(ed, par);
    Local_redraw(ed, par);
end
if nargout ==1
    varargout{1} = par;
end

function Local_redraw(ed, par)
    switch lower(par.Direction(1))
        case 'v'
            ico = Local_icoV;
        case 'h'
            ico = Local_icoH;
    end
    if strcmp(par.EnableButton, 'on')
        out = ico;
        out(:, :, 1) = out(:, :, 1)*par.BackgroundColor(1);
        out(:, :, 2) = out(:, :, 2)*par.BackgroundColor(2);
        out(:, :, 3) = out(:, :, 3)*par.BackgroundColor(3);
    else
        out = ico*0.5+0.5;
        out(:, :, 1) = out(:, :, 1)*par.BackgroundColor(1);
        out(:, :, 2) = out(:, :, 2)*par.BackgroundColor(2);
        out(:, :, 3) = out(:, :, 3)*par.BackgroundColor(3);
    end
    builtin('set',ed.push1,'CDATA', out(end:-1:1, :, :));
    builtin('set',ed.push2,'CDATA', out(:, end:-1:1, :));
            

function out = Local_icoV
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

function out = Local_icoH
out(:, :, 1) = [1   1   1   1   1   1   0   0   1   1;
                1   1   1   1   1   0   0   0   1   1;
                1   1   1   1   0   0   0   0   1   1;
                1   1   1   0   0   0   0   0   1   1;
                1   1   0   0   0   0   0   0   1   1;
                1   1   0   0   0   0   0   0   1   1;
                1   1   1   0   0   0   0   0   1   1;
                1   1   1   1   0   0   0   0   1   1;
                1   1   1   1   1   0   0   0   1   1;
                1   1   1   1   1   1   0   0   1   1];
out(:, :, 2) = out(:, :, 1);
out(:, :, 3) = out(:, :, 1);

function spinboxcb(hands, obj)

if ~isa(hands, 'spinbox'), return; end

param = get(hands.push2, 'userdata');
if strcmp(param.Enable, 'off'), return; end
donothing = 0;
switch obj
    case 'edit'
            param.Value = str2num(get(hands.edit, 'String'));
            param.String = get(hands.edit, 'String');
            if ~isempty(param.Value)
                set(hands, 'EnableButtons', 'on');
            else
                
                set(hands, 'EnableButtons', 'off');
            end
            sss = builtin('get', hands.edit, 'String');
            
    case 'push2' % down
        if ~ischar(param.Value)
            num = param.Value + param.Step;
            param.Value = num;
            param.String = num2str(param.Value, param.Precision);
        end
        builtin('set', hands.edit, 'String', num2str(num, param.Precision));
    case 'push1' % up
        if ~ischar(param.Value)
            num = param.Value - param.Step;
            param.Value = num;
            param.String = num2str(param.Value, param.Precision);
        end
        builtin('set', hands.edit, 'String', num2str(num, param.Precision));
    case 'keypress'
        donothing = 1;
        key = double(builtin('get', gcf, 'CurrentCharacter'));
        if ~isempty(key)
            if key==30 %|| key==29 % up key or right
                if ~ischar(param.Value)
                    num = param.Value + param.Step;
                    param.Value = num;
                    param.String = num2str(num, param.Precision);
                end
                builtin('set', hands.edit, 'String', num2str(num, param.Precision));
                donothing = 0;
            elseif key==31 %|| key==28  % down key or left
                if ~ischar(param.Value)
                    num = param.Value - param.Step;
                    param.Value = num;
                    param.String = num2str(num, param.Precision);
                end
                builtin('set', hands.edit, 'String', num2str(num, param.Precision));
                donothing = 0;
            end
       
        end
end
if ~donothing  
builtin('set', hands.push2, 'userdata', param);
  
    if ~isempty(param.Callback)
        %     global ed;
        ed = hands;
        eval(param.Callback)
    end
end

function hhh = gcbo
% just for case, when in callback is 'gcbo' function
% global ed;
hhh = [];%builtin('gcbo');

function hhh = gco
% just for case, when in callback is 'gco' function
% global ed;
hhh = builtin('gco');