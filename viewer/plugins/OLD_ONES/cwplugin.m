function varargout = cwplugin(varargin)
% CWPLUGIN Application M-file for cwplugin.fig
%    FIG = CWPLUGIN launch cwplugin GUI.
%    CWPLUGIN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 15-Nov-2004 09:19:24

if nargin == 0  % LAUNCH GUI
    token = 'CWPLUGIN_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,n,e] = fileparts(which('kazan'));
    inifilename = [fpath '\kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try opc = ini.KazanViewer.MultyKazan; end
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end

	fig = openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    for kk=1:6
        handles.range(kk).pnts = [];
        handles.range(kk).bracket = [];
    end
    handles.Script = get(handles.lbPars, 'String');
    handles.Pages{1} = {'f=c'};
    handles.Pages{2} = handles.Script;
    handles.Pages{3} = {'f=c*x.^2+b*x+a'};
    handles.Pages{4} = {'f=b*lshape(x,x0,fw,0,1)+a'};
    handles.Pages{5} = {'f=b*lshape(x,x0,fw,1,1)+a'};
    handles.Pages{6} = {'f=c'};
    handles.ActivePage = 2;
    set(handles.pb1,'Value',1)
    handles.Color =   [1, 0, 0;...
            1, 1, 0;...
            0, 1, 0;...
            0, 0, 0;...
            1, 0, 1;...
            0, 1, 1;...
            1, 1, 0];
    handles.sim = [];
    guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

function varargout = figure1_DeleteFcn(h, eventdata, handles, varargin)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.figure1, handles.hh)']);

% --------------------------------------------------------------------
function FileUpdate(handles)
handles.area = [];
str = get(handles.lbPars, 'String');
val = 0;
hand = guidata(handles.hh);

% pnts = x*0;
for jj=1:6
    FitRng = [];
    str = handles.Pages{jj};
    for kk = 1:size(str, 1)
        if strfind(str{kk}, 'FitRng'), eval([str{kk}, ';']); end
    end
    [a,b] = getrange(FitRng, hand.src.ax.x);
    handles.range(jj).pnts = a; handles.range(jj).idx = b;
    handles.range(jj).bracket = FitRng;
end

guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = lbPars_Callback(h, eventdata, handles, varargin)
Res = get(handles.lbPars, 'String');
val = get(handles.lbPars, 'Value');

[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(2:end)), ' '];

if strcmp(trim(name), 'f')
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'function');
elseif str2(1)==''''
    primes = findstr(str2, '''');
    if size(primes, 2) < 2, primes(2) = size(str2, 2); end
    set(handles.ePars, 'String', str2(2:(primes(2)-1)));
    set(handles.ePars, 'UserData', 'string');
else
    dig = str2num(str2);
    set(handles.ePars, 'String', num2str(dig, 6));
    set(handles.ePars, 'UserData', 'value');
end


% --------------------------------------------------------------------
function varargout = ePars_Callback(h, eventdata, handles, varargin)
Res = get(handles.lbPars, 'String');
val = get(handles.lbPars, 'Value');
newval = trim(get(handles.ePars, 'String'));

[name, str1] = strtok(Res{val}, '=');
name = trim(name);

switch get(handles.ePars, 'UserData')
    case 'function'
        set(handles.ePars, 'String', newval);
        
        tokens = kv_autofit([name ' = ' newval]);
        
        %%%%%%%%%%%%% load script %%%%%%%%%%%%%%%%%%%%%%
        for k=1:size(Res,1),
            if ~strcmp(trim(strtok(Res{k},'=')),'f')
                try,eval([Res{k} ';']);end
            end
        end
        tempcell = Res;
        Res = {[name ' = ' newval]};
        for k = 1: size(tokens, 2)
            try, var = eval(tokens{k}); catch, var = 0; end
            Res{end+1} = [tokens{k},' = ', num2str(var)];
        end
        %%%% find strings with additional parameters %%%
        for count = 1:length(tempcell)
            if length(tempcell{count})>6,
                if strcmp(tempcell{count}(1:6), 'FitRng')
                    Res{end+1} = tempcell{count};
                end
            end
        end
        clear tempcell;
    case 'string'
        set(handles.ePars, 'string', newval);
        Res{val} = [name ' = ''', newval,''''];
    otherwise
        nums = str2num(newval);
        set(handles.ePars, 'String', num2str(nums));
        Res{val} = [name ' = ' num2str(nums)];
end
set(handles.lbPars, 'String', Res);
lbPars_Callback(h, eventdata, handles);
handles.Script = Res;
guidata(handles.figure1, handles);
hand = guidata(handles.hh);

if strcmp(get(handles.ePars, 'UserData'), 'value'),
    evallist(handles);

end

% --------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)
parname = get(handles.lbPars, 'String');
parval = get(handles.lbPars, 'Value');
if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;

str = get(handles.pmPars, 'String');

shift = get(handles.slPars, 'Value');
set(handles.slPars, 'Value', 0.5);
num = get(handles.pmPars, 'Value');

dim = str2num([str{num}]);

CurVal = str2double(get(handles.ePars, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.ePars, 'String', num2str(Val));
ePars_Callback(h, eventdata, handles, varargin);



% --------------------------------------------------------------------
function varargout = pmPars_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = pbIntegral_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
if ~isempty(hand.src.y)
    isShift = get(handles.cbShift,   'Value');
    isSquare = get(handles.cbSquare, 'Value');
    if get(handles.cbAntiderivative, 'Value')
        y = cumsum(hand.src.y - handles.sim.y);
    else
        y = hand.src.y - handles.sim.y;
    end
    if isSquare,  y=abs(y); end;
    xcc = (max(hand.src.ax.x)-min(hand.src.ax.x))/(length(hand.src.ax.x)-1);
    if get(handles.cbFullIntRange,'Value')
        res = sum(y,1)*xcc;
        handles.area = 1:size(y,1);
    else
        str1 = get(handles.eFirstPos, 'String');
        str2 = get(handles.eLastPos, 'String');
        p1=str2num(str1);    p2=str2num(str2);
        if isempty(p1), p1=hand.src.ax.x(1); end
        if isempty(p2), p2=hand.src.ax.x(size(handles.sim.y,1)); end
        [a,idx1]=min(abs(hand.src.ax.x-p1));
        [a,idx2]=min(abs(hand.src.ax.x-p2));
        set(handles.eFirstPos, 'String', num2str(hand.src.ax.x(idx1)));
        set(handles.eLastPos, 'String', num2str(hand.src.ax.x(idx2)));
        handles.area = idx1:idx2;
        if isShift, shift = y(min(handles.area));
        else shift = 0;
        end
        res = sum(y(handles.area)-shift,1)*xcc;
    end
    guidata(handles.figure1, handles);
    set(handles.eResInt, 'String', num2str(res,4));
end
plothis(handles)

% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = butFit_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
sliderval = get(hand.sl2Dplot, 'Value');

[prj, selec] = getslider_kazan(hand);

eval(['x  = hand.src.ax.', lower(prj), ';']);
eval(['y = hand.src.y', selec, ';']);
[ss, ind] = max(size(y));
if ind == 2, x = x.'; end
[strings,yres]= kv_autofit(get(handles.lbPars, 'String'), x, real(y), handles.range(handles.ActivePage).pnts);
set(handles.lbPars, 'String', strings);
handles.Script = strings;
guidata(handles.figure1, handles);

if strcmp(selec, '(:, :)'), sel = ''; else sel = selec; end
handles.sim.y = y;
eval(['handles.sim.y', sel,'= yres;']);
try
    set(handles.lbPars, 'String', strings);
catch
    set(handles.lbPars, 'String', strings, 'Value', 1);
end
lbPars_Callback(h, [], handles);

num = get(handles.pmPars, 'Value');
handles.Script{num} = strings;

guidata(handles.figure1, handles);
pbIntegral_Callback([], [], handles)
plothis(handles);
%--------------------------------------------------------------------

function plothis(handles)
handles = guidata(handles.figure1);
colors = handles.Color;

hand = guidata(handles.hh);
ax = hand.src.ax;
y = hand.src.y;

hand.out = {};

% Data
if 111==111
    hand.out{1}.ax = ax;
    hand.out{1}.ax.Color = [1 1 1]*0.4; % always gray ?
    hand.out{1}.ax.s = 0;
    hand.out{1}.y = y;
    hand.out{1}.title = 'Src';
    
    if handles.range(handles.ActivePage).bracket,
        FitRng = handles.range(handles.ActivePage).idx;
        
        Dy = (max(max(real(y))) - min(min(real(y))))/40; 
        eval(['axxx = ax.', lower(prj), ';']);
        for count = 1:size(FitRng, 1)  
            hand.out{end+1}.ax = ax;
            
            if FitRng(count, 1) >length(axxx), FitRng(count, 1) = length(axxx); end
            if FitRng(count, 2) >length(axxx), FitRng(count, 2) = length(axxx); end
            %%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eval(['hand.out{end}.ax.', lower(prj),' = ax.', lower(prj),...
                    '([FitRng(count, 1):FitRng(count, 2)], 1);']);
            hand.out{end}.ax.Color = 'b'; % always blue ?
            hand.out{end}.ax.s = 0;
            hand.out{end}.title = ['Bint', num2str(count)];
            
            switch prj
            case 'X'
                hand.out{end}.y = y([FitRng(count, 1):FitRng(count, 2)], :);                
            case 'Y'
                hand.out{end}.y = y(:, [FitRng(count, 1):FitRng(count, 2)]);
            end
            if 1==1 %handles.range.showbnds,
                %%%%%%%%%%%%%% Bra... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hand.out{end+1}.ax = ax;
                eval(['hand.out{end}.ax.', lower(prj),' = ax.', lower(prj),...
                        '(FitRng(count, [1, 1]), 1);']);
                hand.out{end}.ax.Color = 'b'; % always blue ?
                hand.out{end}.ax.LineWidth = 2;
                hand.out{end}.ax.s = 0;
                hand.out{end}.title = '[..';
                
                switch prj
                case 'X'
                    hand.out{end}.y(1, :) = y(FitRng(count, 1), :) - Dy;
                    hand.out{end}.y(2, :) = y(FitRng(count, 1), :) + Dy;
                case 'Y'
                    hand.out{end}.y(:, 1) = y(:, FitRng(count, 1)) - Dy;
                    hand.out{end}.y(:, 2) = y(:, FitRng(count, 1)) + Dy;
                end        
                %%%%%%%%%%%%%% cket... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hand.out{end+1}.ax = ax;
                eval(['hand.out{end}.ax.', lower(prj),' = ax.', lower(prj),...
                        '(FitRng(count, [2, 2]), 1);']);
                hand.out{end}.ax.Color = [0, 0, 1]; % always blue ?
                hand.out{end}.ax.LineWidth = 2;
                hand.out{end}.ax.s = 0;
                hand.out{end}.title = '..]';
                
                switch prj
                case 'X'
                    hand.out{end}.y(1, :) = y(FitRng(count, 2), :) - Dy;
                    hand.out{end}.y(2, :) = y(FitRng(count, 2), :) + Dy;
                case 'Y'
                    hand.out{end}.y(:, 1) = y(:, FitRng(count, 2)) - Dy;
                    hand.out{end}.y(:, 2) = y(:, FitRng(count, 2)) + Dy;
                end   
            end
        end
    else
        hand.out{1}.ax.Color = [0 0 1]; 
    end
    % Baseline
    hand.out{end + 1}.ax = ax;
    hand.out{end}.y = handles.sim.y;
    hand.out{end}.ax.Color = [1, 0, 0]; % always red
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Bsline';
end


if ~isempty(handles.area)
    area = (handles.area);
else 
    area = 1:length(ax.x);
end

if get(handles.cbAntiderivative, 'Value')
    yy = cumsum(y - handles.sim.y);
    shift = yy(min(area));
else
    yy = y - handles.sim.y;
    shift = 0;
end

% Data - result
hand.out{end+1}.ax = ax;
hand.out{end}.y = yy - shift;
hand.out{end}.ax.Color = [0, 1, 0]; % always green
hand.out{end}.ax.s = 0;
hand.out{end}.title = 'Diff';

% area
if get(handles.cbShowIntRange, 'Value')
    while length(area) > 150
        area = area(1:2:length(area));
    end
    hand.out{end+1}.ax = ax;
    hand.out{end}.ax.x=[];
    hand.out{end}.ax.x(1,:) = ax.x(area);
    hand.out{end}.ax.x(2,:) = ax.x(area);
    hand.out{end}.ax.y = [1;2];
    hand.out{end}.y(1,:) = yy(area)'-shift;
    hand.out{end}.y(2,:) = 0;
    hand.out{end}.ax.Color = [.2, .2, .2]; % always gray
    hand.out{end}.ax.diff = '';
    hand.out{end}.ax.filt = '';
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Area';
end

guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

%--------------------------------------------------------------------
function evallist(handles)
str = get(handles.lbPars, 'String');
val = 0;
hand = guidata(handles.hh);
FitRng = [];
eval(['x = hand.src.ax.', lower(hand.Prj), ';']);
eval(['y = hand.src.y', hand.sel, ';']);
[ss, ind] = max(size(y));
if ind == 2, x = x.'; end
for count = 1:size(str, 1)
    if ~strcmp(str{count}(1:3), 'f =')
        eval([str{count}, ';']);
    else
        val = count;
    end
end

[handles.range(handles.ActivePage).pnts, handles.range(handles.ActivePage).idx] = getrange(FitRng, x);
handles.range(handles.ActivePage).bracket = FitRng;

if ~val, disp('EVALLIST: Can not find function'); return; end
if findstr(str{val}, 'polyfit'), 
    butFit_Callback(handles.butFit, [], handles)
else
    eval([str{val}, ';']);
    if strcmp(hand.sel, '(:, :)'),
        switch hand.Prj
            case 'X'
                eval(['handles.sim.y = f(:, ones(size(hand.src.y, 2), 1));']);
            case 'Y'
                eval(['handles.sim.y = f(ones(size(hand.src.y, 1), 1), :);']);
        end
    else
        eval(['handles.sim.y', hand.sel, ' = f;']);
    end
    guidata(handles.figure1, handles);
    pbIntegral_Callback(handles.ePars, [], handles);
end

% --------------------------------------------------------------------
function [points,FitRng] = getrange(FitRngX, x)
numpoints=length(x);
pnts = zeros(1,numpoints);
sz = size(FitRngX);
allranges = reshape(FitRngX, 1, prod(sz));
[a,idx] = min(abs(allranges(ones(numpoints,1),:)-x(:,ones(1,prod(sz)))));
FitRng = reshape(idx, sz);
if ~isempty(FitRng)
    for i = 1:size(FitRng, 1)
        pnts([min(FitRng(i, 1), end):min(FitRng(i, 2), end)]) = 1;
    end
    points = find(pnts);
else
    points = [];    
end

% --------------------------------------------------------------------
function varargout = mSimContext_Callback(h, eventdata, handles, varargin)
% menu recreating
delete(findobj(h, 'UserData', 'ranges'));
if ~isempty(handles.range(handles.ActivePage).bracket)
    for k=1:size(handles.range(handles.ActivePage).bracket,1)
        uimenu(h, 'Label', ['Delete Range ', num2str(k)], 'UserData', 'ranges',...
            'Tag', 'mPluginsReposition', ...
            'Callback', ['cwplugin(''mContextDelRange_Callback'',gcbo,',num2str(k),',guidata(gcbo))']);
    end
end

% --------------------------------------------------------------------
function varargout = mContextAddRange_Callback(h, eventdata, handles, varargin)
str = get(handles.lbPars, 'String');
chkrng = 0;
for count = 1:length(str)
    if length(str{count})<6, continue; end
    chkrng = strcmp(str{count}(1:6), 'FitRng') + chkrng;
end
hand = guidata(handles.hh);
if isfield(hand, 'src'),
    eval(['numbegs = num2str(hand.src.ax.',lower(hand.Prj),'(1));']);
    eval(['numends = num2str(hand.src.ax.',lower(hand.Prj),'(end));']);
else
    numends  = '1'; numbegs = '1';
end
if ~chkrng, 
    str{end + 1} = ['FitRng(1, 1) = ', numbegs];
    str{end + 1} = ['FitRng(1, 2) = ', numends];
    handles.Script{end+1} = ['FitRng(1, 1) = ', numbegs];
    handles.Script{end+1} = ['FitRng(1, 2) = ', numends];
else
    str{end + 1} = ['FitRng(end+1, 1) = ',numbegs];
    str{end + 1} = ['FitRng(end, 2) = ', numends];
    handles.Script{end+1} = ['FitRng(end+1, 1) = ',numbegs];
    handles.Script{end+1} = ['FitRng(end, 2) = ', numends];
end
handles.Pages{handles.ActivePage} = handles.Script;
guidata(handles.figure1, handles);
set(handles.lbPars, 'String', str);
evallist(handles);

% --------------------------------------------------------------------
function varargout = mContextDelRange_Callback(h, del, handles, varargin)
str = get(handles.lbPars, 'String');
idx = 1:size(handles.range(handles.ActivePage).bracket,1);
idx = idx(idx~=del);
handles.range(handles.ActivePage).bracket = handles.range(handles.ActivePage).bracket(idx,:);

str1=[];
for kk=1:size(str,1)
    if strfind(str{kk},'FitRng'), break; end
    str1{kk}=str{kk};
end
str = str1;

for kk=1:size(handles.range(handles.ActivePage).bracket,1)
    str{end+1}=['FitRng(',num2str(kk),', 1) = ',num2str(handles.range(handles.ActivePage).bracket(kk,1))]; 
    str{end+1}=['FitRng(',num2str(kk),', 2) = ',num2str(handles.range(handles.ActivePage).bracket(kk,2))];
end
set(handles.lbPars, 'String',str, 'Value', 1);
handles.Script = str;

hand = guidata(handles.hh);
[handles.range(handles.ActivePage).pnts, handles.range(handles.ActivePage).idx] = ...
    getrange(handles.range.bracket, hand.src.ax.x);
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function varargout = mContextShowRange_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = eLastPos_Callback(h, eventdata, handles, varargin)
pbIntegral_Callback(h, eventdata, handles)
% --------------------------------------------------------------------
function varargout = eFirstPos_Callback(h, eventdata, handles, varargin)
pbIntegral_Callback(h, eventdata, handles)
% --------------------------------------------------------------------
function varargout = cbFullIntRange_Callback(h, eventdata, handles, varargin)
pbIntegral_Callback(h, eventdata, handles)
% --------------------------------------------------------------------
function varargout = cbShowIntRange_Callback(h, eventdata, handles, varargin)
plothis(handles)
% --------------------------------------------------------------------
function varargout = cbAntiderivative_Callback(h, eventdata, handles, varargin)
pbIntegral_Callback(h, eventdata, handles)

% --------------------------------------------------------------------
function varargout = pbPages_Callback(h, eventdata, handles, varargin)
handles.Pages{handles.ActivePage} = get(handles.lbPars,'String');
hh = [handles.pb0,handles.pb1,handles.pb2,handles.pb3,handles.pb4,handles.pb5];
set(hh, 'Value',0)
set(hh(eventdata+1), 'Value',1)
handles.ActivePage = eventdata+1;
set(handles.lbPars,'String', handles.Pages{eventdata+1},'Value',1);
guidata(handles.figure1, handles);
lbPars_Callback(h, eventdata, handles);

% --------------------------------------------------------------------
function varargout = cbShift_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = cbSquare_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function [prj, sel] = getslider_kazan(hand)
str = {'X', 'Y'};
prj = str{hand.Projection};
if hand.PlotType ==1
    switch hand.Projection
        case 1
            sel = ['(:, ','min(',  num2str(hand.Selection),', end))'];
        case 2
            sel = ['(min(', num2str(hand.Selection),', end), :)'];
    end
else
    sel = '(:, :)';
end

% --------------------------------------------------------------------
function pbLineWidth_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 1)']);
if length(x)~=1, return; end;

hand = guidata(handles.hh);
axx = hand.src.ax.x;
yy  = hand.src.y;

[mm,pidx] = min(abs(axx-x(1)));

for ii=pidx:-1:1, 
    if abs(yy(ii)) < abs(yy(pidx))/2,
        lpidx=ii; break;
    end
end
for ii=pidx:length(axx), 
    if abs(yy(ii)) < abs(yy(pidx))/2,
        rpidx=ii; break;
    end
end
disp(['LineWidth=',num2str(axx(rpidx)-axx(lpidx))])
set(handles.eResInt, 'String', num2str(axx(rpidx)-axx(lpidx)));
lpidx,rpidx

% --------------------------------------------------------------------
function pbLineWidthPP_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 1)']);
if length(x)~=1, return; end;

hand = guidata(handles.hh);
axx = hand.src.ax.x;
yy  = hand.src.y;

[mm,pidx] = min(abs(axx-x(1)));

cc= 0.02*(max(yy)-min(yy));

lpidx = pidx;
for ii=pidx:-1:1, 
    if yy(ii) >= yy(lpidx), lpidx=ii;
    elseif yy(ii) < yy(lpidx)-cc, break;
    end
end
rpidx = pidx;
for ii=pidx:length(axx), 
    if yy(ii) <= yy(rpidx), rpidx=ii;
    elseif yy(ii) > yy(rpidx)+cc, break;
    end
end
disp(['LineWidth=',num2str(axx(rpidx)-axx(lpidx))])
set(handles.eResInt, 'String', num2str(axx(rpidx)-axx(lpidx)));
