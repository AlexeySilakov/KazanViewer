function varargout = fittingbox(varargin)
% FittingBox Application M-file for fittingbox.fig
%    Intended to use only with 'kazan' viewer

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 02-May-2005 13:45:50
% Silakov Alexey
% AlSi 13.03.2004, 15:00

if nargin == 0  % LAUNCH GUI
    token = 'SPECFITBOX_open';
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
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    set(handles.figure1,'name','SPECFITBOX');    
    handles.tbMain = uitoolbar(handles.figure1);
    load([fpath, '\ico.mat']);
    
    handles.ptLoad = uipushtool(handles.tbMain, 'TooltipString', 'Load script from file', 'Tag', 'ptLoad', ...
        'CData', ico.load, 'ClickedCallback', 'SPECFITBOX(''mFile_Callback'',gcbo,guidata(gcbo))');
    handles.ptSave = uipushtool(handles.tbMain, 'TooltipString', 'Save current script', 'Tag', 'ptSave',...
        'CData', ico.save, 'ClickedCallback', 'SPECFITBOX(''mFile_Callback'',gcbo,guidata(gcbo))');
    handles.ptEdit = uipushtool(handles.tbMain, 'TooltipString', 'Change Script in editor', 'Tag', 'ptEdit', ...
        'CData', ico.edit, 'ClickedCallback', 'SPECFITBOX(''mFile_Callback'',gcbo,guidata(gcbo))');
    handles.ptClipboard = uipushtool(handles.tbMain, 'TooltipString', 'Copy to Clipboard', 'Tag', 'ptClipboard',...
        'CData', ico.clipboard,'ClickedCallback', 'SPECFITBOX(''mFile_Callback'',gcbo,guidata(gcbo))');    
   
    handles.hh = 0;
    handles.ploter = 0;
    handles.num = 1;
    handles.div = 0;
    handles.Data = [];
        doppar.Nsubs = 1;
        doppar.harm = 1;
        doppar.dist = 0;
        doppar.maxiter = 100000;
        doppar.fit_amps = 0;
        doppar.amp = 1;
        doppar.x0 = 0;
        doppar.fwhh = 1;
    handles.Script{1} = gen_script(doppar);

    set(handles.list, 'String', handles.Script{1});
    handles.sim.y = [];
    handles.par = doppar;
    handles.loadtype = 'onebyone'; % or 'allinone'  -all in one
    handles.Color =   [1, 0, 0;...
            1, 1, 0;...
            0, 1, 0;...
            0, 0, 0;...
            1, 0, 1;...
            0, 1, 1;...
            1, 1, 0];
    handles.process = 1;
    handles.range.pnts = [];
    handles.range.bracket = [];
    handles.range.showbnds = 0;
    handles.showranges = 1;
    handles.ShowLines = 0;
    handles.ShowSum = 1;
    %   set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
    guidata(fig, handles);
    figure1_ResizeFcn(handles.figure1, [], handles)
    list_Callback([], [], handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        %     [varargout{1:nargout}] = 
        feval(varargin{:}); % FEVAL switchyard
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
% --------------------------------------------------------------------
function varargout = figure1_DeleteFcn(h, eventdata, handles, varargin)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.figure1, handles.hh)']);

% --------------------------------------------------------------------
function varargout = figure1_CloseRequestFcn(h, eventdata, handles, varargin)
delete(handles.figure1);

% --------------------------------------------------------------------
function FileUpdate(varargin);
return;
% --------------------------------------------------------------------
function varargout = list_Callback(h, eventdata, handles, varargin)
% % Stub for Callback of the uicontrol handles.listbox1.
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');

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
function varargout = butLoadData_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
if isempty(hand),  
    msgbox('Can not find main box handles', 'Error');
    return;
end
[path,name,ext] = fileparts(hand.fname);
file_str = name;
sizey = size(hand.src.y, 2);
if sizey==1
    button = 'All in One';
else
    button = questdlg('How to load?',...
    '???','One by One','All in One','Cancel','One by One');        
end
popString = {};
switch button
    case 'One by One'
        if sizey > 10, 
            button = questdlg('Do you realy want to load data one by one?',...
                ['there are ' num2str(sizey), ' columns'],...
                'Oh, ja ja','Nie','Nie');
            if strcmp(button, 'Nie'), return; end
        end
        handles.loadtype = 'onebyone';
        for i = 1:sizey
            handles.Script{i} = get(handles.list, 'String');
            popString{end+1} = [file_str, '_col_ ', num2str(i)];
        end
    case 'All in One'
        handles.loadtype = 'allinone';
        handles.Script{1} = get(handles.list, 'String');
        popString = {file_str};
    otherwise
        return;
end
handles.sim.y = hand.src.y*0;
handles.src = hand.src;
% if ~isfield(handles, 'range'),
%     switch hand.Prj
%         case 'X'
%             handles.range.pnts = [1:size(hand.src.y, 1)];
%             handles.range.barcket = [1, size(hand.src.y, 1)];
%         case 'Y'
%             handles.range.pnts = [1:size(hand.src.y, 2)];
%             handles.range.barcket = [1, size(hand.src.y, 2)];
%     end
% end

set(handles.popSel, 'String', popString, 'Value', 1);
set(handles.list, 'String', handles.Script{1}, 'Value', 1);
list_Callback(handles.list, [], handles);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = ePars_Callback(h, eventdata, handles, varargin)
% % Stub for Callback of the uicontrol handles.ePars.
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');
newval = trim(get(handles.ePars, 'String'));


[name, str1] = strtok(Res{val}, '=');
name = trim(name);
if ~isempty(findstr(newval, 'lshapes')), lshapes(handles); return; end

    
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

        clear tempcell;
    case 'string'
        set(handles.ePars, 'string', newval);
        Res{val} = [name ' = ''', newval,''''];
    otherwise
        nums = str2num(newval);
        set(handles.ePars, 'String', num2str(nums));
        Res{val} = [name ' = ' num2str(nums)];
end
set(handles.list, 'String', Res);
if strcmp(Res{val}(1:3), 'N ='), lshapes(handles); end
list_Callback(h, eventdata, handles);
switch handles.loadtype
    case 'onebyone', 
        handles.Script{get(handles.popSel, 'Value')} = Res;
    case 'allinone'
        handles.Script{1} = Res;
end
handles.par = gen_par(Res);
guidata(handles.figure1, handles);
hand = guidata(handles.hh);

if strcmp(get(handles.ePars, 'UserData'), 'value'),
    evallist(handles);
end

%--------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)
parname = get(handles.list, 'String');
parval = get(handles.list, 'Value');
if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;

str = get(handles.popup, 'String');

shift = get(handles.slPars, 'Value');
num = get(handles.popup, 'Value');

dim = str2num([str{num}]);

CurVal = str2double(get(handles.ePars, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.ePars, 'String', num2str(Val));
ePars_Callback(h, eventdata, handles, varargin);

%--------------------------------------------------------------------
function plothis(handles)
handles = guidata(handles.figure1);
colors = handles.Color;

hand = guidata(handles.hh);
ax = handles.src.ax;
y = handles.src.y;

[prj, selec] = getslider_kazan(hand);
hand.out = {};
% ax = hand.src.ax;
% y = hand.src.y;

if handles.process, [y, ax] = processing(y, ax);
else
    ax.diff = '';
    ax.filt = '';
end
% Data
hand.out{1}.ax = ax;
hand.out{1}.ax.Color = [1 1 1]*0.4; % always gray ?
hand.out{1}.ax.s = 0;
hand.out{1}.y = y;
hand.out{1}.title = 'Src';

if handles.range.bracket,
    FitRng = handles.range.bracket;
    
    Dy = (max(max(real(y))) - min(min(real(y))))/50; 
    eval(['axxx = ax.', lower(prj), ';']);
    for count = 1:size(FitRng, 1)  
        hand.out{end+1}.ax = ax;
        
        if FitRng(count, 1) >length(axxx), FitRng(count, 1) = length(axxx); end
        if FitRng(count, 2) >length(axxx), FitRng(count, 2) = length(axxx); end
        %%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eval(['hand.out{end}.ax.', lower(prj),' = ax.', lower(prj),...
                '([FitRng(count, 1):FitRng(count, 2)], 1);']);
        hand.out{end}.ax.Color = [0, 0, 1]; % always blue ?
        hand.out{end}.ax.s = 0;
        hand.out{end}.title = ['Src', num2str(count)];
        
        switch prj
            case 'X'
                hand.out{end}.y = y([FitRng(count, 1):FitRng(count, 2)], :);                
            case 'Y'
                hand.out{end}.y = y(:, [FitRng(count, 1):FitRng(count, 2)]);
        end
        if handles.range.showbnds,
            %%%%%%%%%%%%%% Bra... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hand.out{end+1}.ax = ax;
            eval(['hand.out{end}.ax.', lower(prj),' = ax.', lower(prj),...
                    '(FitRng(count, [1, 1]), 1);']);
            hand.out{end}.ax.Color = 'b';
            hand.out{end}.ax.s = 0;
            hand.out{end}.title = ['br_', num2str(count)];
            
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
            hand.out{end}.ax.Color = 'b';
            hand.out{end}.ax.s = 0;
            hand.out{end}.title = ['er_', num2str(count)];
            
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
    hand.out{1}.ax.Color = 'b'; 
end

if handles.ShowLines
    eval(['axxx = ax.', lower(prj), ';']);
    yy = SepLines(axxx(:), handles.par);
    for i = 1:handles.par.Nsubs
        hand.out{end+1}.ax = ax;
        hand.out{end}.y = yy(:, i);
        hand.out{end}.ax.Color = [1, 1, 1]*0.5;
        hand.out{end}.ax.s = 0;
        hand.out{end}.title = ['Line_', num2str(i)];
    end
end
% Result
if handles.ShowSum
    hand.out{end + 1}.ax = ax;
    hand.out{end}.y = handles.sim.y;
    hand.out{end}.ax.Color = 'r';
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Out';
end
states = {'off', 'on'};
st = get(handles.mSimDif, 'Checked');
check = find(strcmp(states, st)) -1;
% diff = Data - result
if check, 
    hand.out{end+1}.ax = ax;
    hand.out{end}.y = y - handles.sim.y;
    hand.out{end}.ax.Color = 'g';
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Diff';
end

guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% ------------------------------------------------------------
function butCol_Callback(hObject, eventdata, handles)
num = get(handles.popSel, 'Value');
c = uisetcolor(handles.Color(num, :));
if size(c, 2) > 1
  
  handles.Color(num, :) = c;
  set(handles.butCol, 'BackgroundColor', c, 'ForegroundColor', [1 1 1] - c);
  guidata(handles.figure1, handles);
  plothis(handles);
end


% % ------------------------------------------------------------
function popSel_Callback(hObject, eventdata, handles)
num = get(handles.popSel, 'Value');
if strcmp(handles.loadtype, 'onebyone'),
    list = handles.Script{num};
    set(handles.list, 'String', list, 'Value', 1);
    
    hand = guidata(handles.hh);
    num = get(handles.popSel, 'Value');
    
    set(hand.sl2Dplot, 'Value', num);
    kazan('sl2Dplot_Callback',hand.sl2Dplot,[],hand);
else
    set(handles.list, 'String', handles.Script{num}, 'Value', 1);
end
list_Callback(handles.list, eventdata, handles);

% set(handles.butCol, 'BackgroundColor', handles.Color(num, :)); 

% --------------------------------------------------------------------
function varargout = nearme(X, xi)
% nearme(X, xi)
% k = nearme(X, xi)
% [k, pbshowhide] = nearme(X, xi)
% 
%   return value (pbshowhide) and index (k) 
%   of nearest to xi number of array X 

[c, k] = min(abs(X - ones(size(X))*xi));
v = X(k);
switch nargout 
    case 0
        varargout{1} = k;
    case 1
        varargout{1} = k;
    case 2
        varargout{1} = k;
        varargout{2} = v; 
end
%--------------------------------------------------------------------

function butFitt_Callback(h, eventdata, handles)
try
    hand = guidata(handles.hh);
    sliderval = get(hand.sl2Dplot, 'Value');
    
    [prj, selec, trans] = getslider_kazan(hand);
    
    if ~isfield(handles,'src'), error('Please load your data!'); end
    
    eval(['x  = handles.src.ax.', lower(prj), ';']);
    eval(['y = handles.src.y', selec, trans, ';']);
    [ss, ind] = max(size(y));
    if ind == 2, x = x.'; end
    num = get(handles.popSel, 'Value');
    
    isitlshapes = 0;
    for i =1:length(handles.Script{num})
        if ~isempty(findstr(handles.Script{num}{i}, 'lshapes')),
            isitlshapes = 1;
        end
    end
    
    if isitlshapes,
        set(handles.figure1,'name','SPECFITBOX: Busy');
        [par, yres] = specfit(handles.Script{num}, x, real(y), 1); %do fitt
        set(handles.figure1,'name','SPECFITBOX');
        strings = gen_script(par);
        handles.par = par;
        
    else 
        [strings,yres]= kv_autofit(handles.Script{num}, x, real(y), handles.range.pnts); 
    end
    if strcmp(selec, '(:, :)'), sel = ''; else sel = selec; end
    eval(['handles.sim.y', sel,'= yres',trans,';']);
    try
        set(handles.list, 'String', strings);
    catch
        set(handles.list, 'String', strings, 'Value', 1);
    end
    list_Callback(h, [], handles);
    
    handles.Script{num} = strings;
    
    guidata(handles.figure1, handles);
    plothis(handles);
catch
     set(handles.figure1,'name','SPECFITBOX: Error');
    if ~isempty(lastwarn), disp(['Warning: ', lastwarn]); end
    if ~isempty(lasterr), disp(['Error: ', lasterr]); end
end

%--------------------------------------------------------------------
function butFittAll_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
% sliderval = get(hand.sl2Dplot, 'Value');
[prj, selec] = getslider_kazan(hand);

switch hand.Projection
case 2
    x = handles.src.ax.y;
    y = real(handles.src.y)';
case 1
    x = handles.src.ax.x;
    y = real(handles.src.y);
end

numy = size(handles.src.y, ~(hand.Projection-1)+1);
switch handles.loadtype
case 'onebyone'
    for i = 1:numy
        set(handles.popSel, 'Value', i);
        set(handles.list, 'String', handles.Script{i});
        eval(['y = handles.src.y', selec, ';']);
        [ss, ind] = max(size(y));
        if ind == 2, x = x.'; end
        [par, yres] = specfit(handles.Script{1}, x, real(y), 0);
        script = gen_script(par);
        handles.par = par;
        eval(['handles.sim.y', selec,'= yres;']);
        set(handles.list, 'String', script);
        handles.Script{i} = script;
    end
case 'allinone'
    %%%%%%%% construct fitting function %%%%%%%
    strings = get(handles.list, 'String');
    
    [par, handles.sim.y] = specfit(strings, x, y, 0);
    script = gen_script(par);
    handles.par = par;
    if hand.Projection==2, handles.sim.y = handles.sim.y'; end
    
    set(handles.list, 'String', script);
    handles.Script{1} = script;
end

guidata(handles.figure1, handles);
plothis(handles);

set(handles.popSel, 'Value', 1);
popSel_Callback(handles.popSel, [], handles);
    %--------------------------------------------------------------------
function evallist(handles)
str = get(handles.list, 'String');
str_range = get(handles.listbox2, 'String');
val = 0;
hand = guidata(handles.hh);
FitRng = [];
% AlSi 02.12.04
[prj, selec] = getslider_kazan(hand);
eval(['x = handles.src.ax.', lower(prj), ';']);
eval(['y = handles.src.y', selec, ';']);

[ss, ind] = max(size(y));
if ind == 2, x = x.'; end
for count = 1:size(str, 1)
    if ~strcmp(str{count}(1:3), 'f =')
        eval([str{count}, ';']);
    else
        val = count;
    end
end
if ~ischar(str_range)
    for count = 1:size(str_range, 1)
        eval([str_range{count}, ';']);
    end
end
% pnts = x*0;
handles.range.pnts = getrange(FitRng, length(x));
handles.range.bracket = FitRng;
% if ~isempty(FitRng)
%     for i = 1:size(FitRng, 1)
%         pnts([FitRng(i, 1):min(FitRng(i, 2), end)]) = 1;
%     end
%     handles.range.pnts = find(pnts);
%     handles.range.bracket = FitRng;
% end

if ~val, disp('EVALLIST: Can not find function'); return; end
if findstr(str{val}, 'polyfit'), 
    butFitt_Callback(handles.butFitt, [], handles)
else
    if findstr(str{val}, 'lshapes'), 
        %doppar = gen_par(str);
        f = lshapes_func(x,handles.par); % don't fit, just calculate
    else
        eval([str{val}, ';']);
    end
    if hand.PlotType~=1,
        switch hand.Projection
            case 1
                handles.sim.y = f(:, ones(size(handles.src.y, 2), 1));
            case 2
                handles.sim.y = f(ones(size(handles.src.y, 1), 1), :);
        end
    else
        eval(['handles.sim.y', selec, ' = f;']);
    end
    guidata(handles.figure1, handles);
    plothis(handles);
end

%--------------------------------------------------------------------
function mSimDif_Callback(h, eventdata, handles)
states = {'off', 'on'};
st = get(handles.mSimDif, 'Checked');
check = find(strcmp(states, st))-1;

set(handles.mSimDif, 'Checked', states{~check + 1});
guidata(handles.figure1, handles);
plothis(handles);

%--------------------------------------------------------------------
function mSimFFunc_Callback(h, eventdata, handles)
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');

fittingbox_func(handles);
uiwait(handles.figure1);
str = get(handles.list, 'UserData');
if isempty(str), return; end
%%%% find strings with additional parameters %%%
for count = 1:length(Res)
    if length(Res{count})>6,
        if strcmp(Res{count}(1:6), 'FitRng')
            str{end+1} = Res{count};
        end
    end
end

set(handles.list, 'String', str, 'Value', 1);
list_Callback(handles.list, [], handles);

% --------------------------------------------------------------------
function mSimAgr_Callback(h, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimAgr, 'Checked'), 'on');
set(handles.mSimAgr, 'Checked', str{~stat + 1});
handles.process = ~stat;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function mSimAddRng_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
axx =hand.out{1}.ax.x;
axy =hand.out{1}.y;

[x,y]=eval([name '(''GetPoint'', hand, 1)']);
tx(1) = nearme(axx, x); % in points
Ylim = get(gca, 'YLim');
linehand = line([1, 1]*axx(max(tx(1), 1)), Ylim, 'LineStyle', ':', 'Parent', gca); 

[x,y]=eval([name '(''GetPoint'', hand, 1)']);
tx(2) = nearme(axx, x); % in points
delete(linehand);
str = get(handles.listbox2, 'String');
chkrng = 0;
if ischar(str)
    chkrng = 0;    
elseif length(str{1})<6
    chkrng = 0;    
elseif ~strcmp(str{1}(1:6), 'FitRng')
    chkrng = 0;    
else 
    chkrng = 1;    
end
hand = guidata(handles.hh);
[prj, sel] = getslider_kazan(hand);
tx = sort(tx);

if isfield(handles, 'src'),
    eval(['numends = (length(handles.src.ax.',lower(prj),'));']);
else
    numends  = 1;
end
if ~chkrng, 
    str{end + 1} = ['FitRng(1, 1) = ', num2str(max(1, tx(1)))];
    str{end + 1} = ['FitRng(1, 2) = ', num2str(min(numends, tx(2)))];
    set(handles.pbKillRange, 'Enable', 'on');
else
    str{end + 1} = ['FitRng(end+1, 1) = ' num2str(max(1, tx(1)))];
    str{end + 1} = ['FitRng(end, 2) = ', num2str(min(numends, tx(2)))];
end
guidata(handles.figure1, handles);
set(handles.listbox2, 'String', str, 'Value', length(str));
evallist(handles);

% --------------------------------------------------------------------
function mSimShowBds_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
if ~handles.range.bracket, disp('no range settings'); return; end
str = {'off', 'on'};
stat = get(handles.mSimShowBds, 'Checked');
num = find(strcmp(str, stat))-1;
set(handles.mSimShowBds, 'Checked', str{~num+1});
handles.range.showbnds = ~num;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function mSimDelRng_Callback(h, eventdata, handles)
if isempty(handles.range.bracket), disp('nothing to delete'); 
    set(handles.pbKillRange, 'Enable', 'off');
    return; 
end

list  = get(handles.listbox2, 'String');
num  = get(handles.listbox2, 'Value');
count = [];
FitRng = [];
% find strings with 'FitRng'
if ~ischar(list)
    for i = 1:length(list)
        eval([list{i}, ';']);
    end
end
if isempty(FitRng), 
    disp('nothing to delete'); 
    set(handles.pbKillRange, 'Enable', 'off');
    return; 
end

if length(list) ==2
    tlist = [];
    set(handles.pbKillRange, 'Enable', 'off');
elseif mod(length(list), 2) % just in case
    error('Odd number of stings ????');
else
    nn = (mod(num, 2)-0.5)*2; % if odd ->-1, else +1
    non = sort([num, num+nn]);
    I = [1:(non(1)-1), (non(2)+1):length(list)];
    tlist = list(I);

end
newsize = length(tlist);
%if (length(count) == 2)|(pos<0),
set(handles.listbox2, 'String', tlist, 'Value', min(num, newsize));
listbox2_Callback(h, [], guidata(handles.figure1));
handles.Script_Range = tlist;
guidata(handles.figure1, handles);
evallist(handles);
% --------------------------------------------------------------------
function points = getrange(FitRng, numpnts)
pnts = [1:numpnts]*0;
if FitRng
    for i = 1:size(FitRng, 1)
        pnts([min(FitRng(i, 1), end):min(FitRng(i, 2), end)]) = 1;
    end
    points = find(pnts);
else
    points = [];    
end



% --------------------------------------------------------------------
function varargout = popup_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = but1Step_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------

function [prj, sel, trans] = getslider_kazan(hand)
str = {'X', 'Y'};
trans = '';
prj = str{hand.Projection};
if hand.PlotType ==1
    switch hand.Projection
        case 1
            sel = ['(:, ','min(',  num2str(hand.Selection),', end))'];
        case 2
            sel = ['(min(', num2str(hand.Selection),', end), :)'];
            trans = '''';
    end
else
    sel = '(:, :)';
end

% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
brd = 5;
YPos = brd;
butw = 17;
sl = [1.2, 1]*butw;
figpos = get(handles.figure1, 'Position'); % points
inn = 0;
if figpos(3)<=150, figpos(3) = 150; inn = 1; end
if figpos(4)<=300, figpos(4) = 300; inn = 1; end
if inn,  set(handles.figure1, 'Position', figpos); end
% buttons
buth = butw;
% bFAp = [brd, brd,figpos(3)-2*brd, buth];
bFp = [brd, brd,figpos(3)-2*brd, buth];
% bFp = bFAp;
% bFp(2) = bFAp(2)+bFAp(4)+brd;
% b1Sp =bFAp; 
% b1Sp(2) = bFp(2)+bFp(4)+brd;
b1Sp = bFp;
%set(handles.butFittAll, 'position', bFAp);
set(handles.butFitt, 'position', bFp);
%set(handles.but1Step, 'position', b1Sp);

Ypos = b1Sp(2)+b1Sp(4)+brd;
% fitting range
if safeget(handles, 'showranges', 1);
    parth = 0.2; % 30% of the Figure height
    
    fRp = [brd, Ypos, figpos(3)-2*brd, figpos(4)*parth];
    slRSp = [fRp(1)+brd, fRp(2)+brd, sl];
    edRVp = [slRSp(1)+slRSp(3) + brd, fRp(2)+brd, fRp(3) - slRSp(3)-2*butw-4*brd, butw];
    pbARp = [edRVp(1) + edRVp(3)+brd, fRp(2)+brd, butw, butw];
    pbKRp = [pbARp(1) + butw, fRp(2)+brd, butw, butw];
    lb2p =  [fRp(1)+brd, fRp(2)+butw + 2*brd, fRp(3)-2*brd, fRp(4)-4*brd - butw];
    txtRp = [fRp(1)+brd, fRp(2)+fRp(4)-8, 56, 12];
    pbSHp =  [fRp(1)+fRp(3)-30, fRp(2)+fRp(4)-7.5, 30, 7.5];
    
    set(handles.frRange, 'position', fRp);
    set(handles.slRangeStep, 'position', slRSp, 'Visible', 'on');
    set(handles.edRangeVal, 'position', edRVp, 'Visible', 'on');
    set(handles.pbAddRange, 'position', pbARp, 'Visible', 'on');
    set(handles.pbKillRange, 'position', pbKRp, 'Visible', 'on');
    set(handles.listbox2, 'position', lb2p, 'Visible', 'on');
    set(handles.txtRange, 'position', txtRp);
    set(handles.pbShowHide, 'position', pbSHp);
else
    fRp = [brd, Ypos, figpos(3)-2*brd, 12];
    txtRp = [fRp(1)+brd, fRp(2)+fRp(4)-8, 56, 10];
    pbSHp =  [fRp(1)+fRp(3)-30, fRp(2)+fRp(4)-7.5, 30, 7.5];
    set(handles.frRange, 'position', fRp);
    set(handles.txtRange, 'position', txtRp);
    set(handles.pbShowHide, 'position', pbSHp);
    set(handles.slRangeStep,'Visible', 'off');
    set(handles.edRangeVal, 'Visible', 'off');
    set(handles.pbAddRange,'Visible', 'off');
    set(handles.pbKillRange, 'Visible', 'off');
    set(handles.listbox2, 'Visible', 'off');
end
% Fitting settings
Ypos = fRp(2)+fRp(4)+brd;
frMp = [brd, Ypos, figpos(3)-2*brd, figpos(4)-Ypos-butw-4*brd];
ePp = [frMp(1)+brd, frMp(2)+brd, frMp(3)-2*brd, butw];
slPp = [frMp(1)+brd, frMp(2)+butw+2*brd, sl];
ppp = [frMp(1)+ 2*brd+sl(1), frMp(2)+butw+2*brd, frMp(3)-3*brd-sl(1), butw];
lstp = [frMp(1)+brd, frMp(2)+3*brd+2*butw, frMp(3)-2*brd, frMp(4)-4*brd-2*butw];
txtMp = [frMp(1)+brd, frMp(2)+frMp(4)-5, 56, 10];

set(handles.frMain, 'position', frMp);
set(handles.ePars, 'position', ePp);
set(handles.slPars, 'position', slPp);
set(handles.popup, 'position', ppp);
set(handles.list, 'position', lstp);
set(handles.txtMain, 'position', txtMp);

Ypos = frMp(2)+frMp(4)+brd;
% Data selector
frChp = [brd, Ypos, figpos(3)-2*brd, butw+2*brd];
popSp = [frChp(1)+brd, frChp(2)+brd, frChp(3)-3*brd - 25, butw];
butCp = [frChp(1)+2*brd + popSp(3), frChp(2)+brd, 25, butw];
txtDp = [frChp(1)+brd, frChp(2)+frChp(4)-5, 56, 10];
set(handles.frChoose, 'position', frChp);
set(handles.popSel, 'position', popSp);
set(handles.butCol, 'position', butCp);
set(handles.txtDatas, 'position', txtDp);

% -----------------------------------------------------------------
function pbShoHide_Callback(h, eventdata, handles, varargin)
handles.showranges = ~handles.showranges;
if strcmp(get(handles.pbShowHide, 'String'), 'v')
    set(handles.pbShowHide, 'String', 'A');
else
    set(handles.pbShowHide, 'String', 'v');
end
guidata(handles.figure1, handles);
figure1_ResizeFcn(handles.figure1, [], handles)

% -----------------------------------------------------------------
function listbox2_Callback(h, eventdata, handles, varargin)
Res = get(handles.listbox2, 'String');
val = get(handles.listbox2, 'Value');
if val==0, set(handles.edRangeVal, 'String', 0); return; end
[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(2:end)), ' '];

dig = str2num(str2);
set(handles.edRangeVal, 'String', num2str(dig, 6));
set(handles.edRangeVal, 'UserData', 'value');

% -----------------------------------------------------------------
function edRangeVal_Callback(h, eventdata, handles, varargin)
% % Stub for Callback of the uicontrol handles.ePars.
Res = get(handles.listbox2, 'String');
val = get(handles.listbox2, 'Value');
newval = trim(get(handles.edRangeVal, 'String'));

[name, str1] = strtok(Res{val}, '=');
name = trim(name);

nums = str2num(newval);
set(handles.edRangeVal, 'String', num2str(nums));
Res{val} = [name ' = ' num2str(nums)];

set(handles.listbox2, 'String', Res);
listbox2_Callback(h, eventdata, handles);

handles.Script_Range = Res;

guidata(handles.figure1, handles);
hand = guidata(handles.hh);

evallist(handles);

%--------------------------------------------------------------------
function varargout = slRangeStep_Callback(h, eventdata, handles, varargin)
parname = get(handles.listbox2, 'String');
parval = get(handles.listbox2, 'Value');
if ~strcmp(get(handles.edRangeVal, 'UserData'),'value'), return; end;

shift = get(handles.slRangeStep, 'Value');
dim = 1;

CurVal = str2double(get(handles.edRangeVal, 'String'));
Val = (CurVal + dim*(shift - 0.5)*2);

set(handles.edRangeVal, 'String', num2str(Val));
edRangeVal_Callback(h, eventdata, handles, varargin);

%--------------------------------------------------------------------
function mFile_Callback(h, handles)
switch h
    case {handles.mFileLoad, handles.ptLoad}
        [fname,pname] = uigetfile('*.m','Open Script');
        if fname == 0, return; end
        ScriptName = [pname '\' fname];
        fid = fopen(ScriptName, 'r');
        [fpath,fname,ext] = fileparts(ScriptName); 
        if fid == -1, disp('File not found'); return; end
        str = {};
        while feof(fid) < 1,
            st = fgetl(fid);
            %if ~isempty(st) & st(1) ~= '%',
                str{end+1} = st;
                %end
        end
        fclose(fid);
        try
            handles.par = gen_par(str);
        catch
            error('Wrong input scipt');
        end
        set(handles.list, 'String', str, 'Value', 1);
        guidata(handles.figure1, handles);
        list_Callback(hObject, [], handles);
    case {handles.mFileSave, handles.ptSave}
      
        [fname,pname] = uiputfile('*.m','Save Script');
        if fname == 0, return; end
        [fpath,fname,ext] = fileparts([pname, fname]); 
        if strcmp(ext ,'') , ext = '.m'; end
        fullname = [fpath, '\', fname, ext];
        Scr = get(handles.list, 'String');
        fid = fopen(fullname, 'w');
        if fid > 0
            for kk = 1:length(Scr)
                fprintf(fid, '%s\r\n', Scr{kk});
            end
            fclose(fid);
        else
            error(['Can''t write to file ''',fullname,'''']);
        end
    case handles.ptEdit
        outscr = kv_scripteditor(get(handles.list,'String'));
        % open(handles.sim{num}.ScriptName);
        if isempty(outscr), return; end
        try
            handles.par = gen_par(outscr);
        catch
            error('Wrong input scipt');
        end
        set(handles.list, 'String', outscr, 'Value', 1);
        handles.sim.Script = outscr;
        guidata(handles.figure1, handles);
        list_Callback(h, [], handles);
    case handles.ptClipboard

        Scr = get(handles.list, 'String');
        out = '';
        for kk = 1:length(Scr)
            out = sprintf( '%s%s\n', out, Scr{kk});
        end
        clipboard('copy',out);
        disp(['SPECFITBOX: script have been copied to the clipboard.']);
    end

%--------------------------------------------------------------------
function [doppar, y]= specfit(script, x, data, fit)
% script must look like this:
% f=lshapes(N)
% N = 2
% fit_amps = 0
% ....
% amp(1) = 123
% x0(1) = 12
% fwhh(1) = 2
% ....
maxiter = 100000;
largescale = 'off';

for i = 1:length(script)
    if strcmp(script{i}(1:3), 'f ='), continue;end
    eval([script{i}, ';']);
end

doppar.Nsubs = N;
doppar.harm = harm;
doppar.dist =  distortion;
doppar.maxiter = maxiter;
doppar.fit_amps = fit_amps;
doppar.amp = amp;
options = optimset;
options.MaxIter = maxiter*fit; % if fit is 0, then it returns just start-point picture
option.LargeScale = largescale;

if fit_amps, 
    %linefit_fitamp(par, x, data, doppar)
    start_param = [x0, fwhh, amp];
    param = fminsearch(@linefit_fitamp, start_param, [], x, data, doppar);
   
    doppar.x0 = param(1:N);
    doppar.fwhh = param([1:N] + N);
    doppar.amp = param([1:N] + 2*N);
    y = lshapes_func(x, doppar);
else
    start_param = [x0, fwhh];
    param = fminsearch(@linefit_noamp, start_param, [], x, data, doppar);
    
    doppar.x0 = param(1:N);
    doppar.fwhh = param([1:N] + N);
    y = lshapes_func(x, doppar);
end


%--------------------------------------------------------------------
function f=linefit_noamp(par, x, data, doppar)
% par - struct: [xo[1], xo[2], fwhh[1], fwhh[2]]
% doppar - structure. doppar.amp, doppar.harm, etc.
% convention: all stuff here must be columns
% accept par!

% dist*Gaussian + (1-dist)*Lorentzian.
N = doppar.Nsubs;

if length(par)~=2*N, error('number of the parameters does not match to number of subfunctions'); end
doppar.x0 = par(1:N);
doppar.fwhh = par([1:N]+N);
y=lshapes_func(x, doppar);
f = sum((data - y).^2);

%--------------------------------------------------------------------
function f=linefit_fitamp(par, x, data, doppar)
% par - struct: [xo[1], xo[2], fwhh[1], fwhh[2]]
% doppar - structure. doppar.amp, doppar.harm, etc.
% convention: all stuff here must be columns
% accept par!

% dist*Gaussian + (1-dist)*Lorentzian.
N = doppar.Nsubs;

if length(par)~=3*N, error('number of the parameters does not match to number of subfunctions'); end
doppar.x0 = par(1:N);
doppar.fwhh = par([1:N]+N);
doppar.amp = par([1:N]+2*N);
y=lshapes_func(x, doppar);
f = sum((data - y).^2);

function y=lshapes_func(x, doppar)
% amplitudes are part of fitting parameters
N = doppar.Nsubs;
y = x(:)*0;
for i = 1:N
    xo = doppar.x0(i);
    fwhh = doppar.fwhh(i);
    amp = doppar.amp(i);
    if fwhh<=0, fwhh=1e-16; end
    tmpy = lshape(x, xo, fwhh, doppar.harm, doppar.dist(i));
    y = amp*(tmpy/max(tmpy)) + y;
end

%--------------------------------------------------------------------
function out_script = gen_script(doppar)
out_script = {'f = lshapes(N)',...
    ['N = ', num2str(doppar.Nsubs)],...
    ['harm = ', num2str(doppar.harm)],...
    ['maxiter = ', num2str(doppar.maxiter)],...
    ['fit_amps = ', num2str(doppar.fit_amps)]...
             };
         
for i = 1:doppar.Nsubs
    out_script{end+1} = ['%-----------'];
    out_script{end+1} = ['amp(', num2str(i), ') = ', num2str(doppar.amp(i))];    
    out_script{end+1} = ['x0(', num2str(i), ') = ', num2str(doppar.x0(i))];
    out_script{end+1} = ['fwhh(', num2str(i), ') = ', num2str(doppar.fwhh(i))];
    out_script{end+1} = ['distortion(', num2str(i), ') = ', num2str(doppar.dist(i))];
end

%--------------------------------------------------------------------
function doppar = gen_par(script)

for i = 1:length(script)
    if strcmp(script{i}(1:3), 'f ='), continue;end
    eval([script{i}, ';']);
end
doppar.Nsubs = N;
doppar.harm = harm;
doppar.dist =  distortion;
doppar.maxiter = maxiter;
doppar.fit_amps = fit_amps;
doppar.amp = abs(amp);
doppar.fwhh = abs(fwhh);
doppar.x0 = x0;

function lshapes(handles)

par = gen_par(get(handles.list, 'string'));

nn = length(par.x0);
amp = ones(1, par.Nsubs);
x0 = zeros(1, par.Nsubs);
fwhh = ones(1, par.Nsubs);
dist = zeros(1, par.Nsubs);

len = min(nn, par.Nsubs);
amp(1:len) = par.amp(1:len);
x0(1:len) = par.x0(1:len);
fwhh(1:len) = par.fwhh(1:len);
dist(1:len) = par.dist(1:len);

par.amp = amp;
par.x0 = x0;
par.fwhh = fwhh;
par.dist = dist;

handles.Script{1} = gen_script(par);
set(handles.list, 'String', handles.Script{1});
guidata(handles.figure1, handles);

%--------------------------------------------------------------------
function mSimShowLns_Callback(h, handles)
str = {'off', 'on'};
stat = get(handles.mSimShowLns, 'Checked');
num = find(strcmp(str, stat))-1;
set(handles.mSimShowLns, 'Checked', str{~num+1});

handles.ShowLines = ~num;
guidata(handles.figure1, handles);
plothis(handles);

%--------------------------------------------------------------------
function mSimShowSum_Callback(h, handles)
str = {'off', 'on'};
stat = get(handles.mSimShowSum, 'Checked');
num = find(strcmp(str, stat))-1;
set(handles.mSimShowSum, 'Checked', str{~num+1});

handles.ShowSum = ~num;
guidata(handles.figure1, handles);
plothis(handles);


%--------------------------------------------------------------------
function yy = SepLines(x, par)
% amplitudes are part of fitting parameters
N = par.Nsubs;
nn = length(x);
yy = zeros(nn, N);
for i = 1:N
    xo = par.x0(i);
    fwhh = par.fwhh(i);
    amp = par.amp(i);
    if fwhh<=0, fwhh=1e-16; end
    yy(:, i) = amp*(lshape(x, xo, fwhh, par.harm, par.dist(i)));
end


% --------------------------------------------------------------------
function mAbout_Callback(hObject, eventdata, handles)
% handles = guidata(handles.figure1);

% if hObject==handles.mAboutHelp
%     [fpath,n,e,v] = fileparts(which('kazan'));
%     filename = ['file:///',fpath,'\help\specfitbox.html'];
%     filename(filename=='\')='/';
%     stat = web(filename);
% else
    msgbox({'simplugin';
        'Interface for fitting spectra by';
        'lorenzian with gaussian distortion'
        'by Alexey Silakov, 21.09.2005';
        'Only for using with Kazan Viewer';
        'EasySpin (www.esr.ethz.ch) function lshape is used';'';
        'Parameters are use to be generated automaticaly,';
        'nevertheless there is possibility to load script from file';
        'Saving is also possible, but it is not recomended to change saved file manualy';'';
        'Parameters:';
        '   "N" - number of lines to fit';
        '       changing this parameter leads to reloading the script';
        '       according to the new number of lines';
        '   "harm" - [0] - absorbtion, [1]  - first derivative';
        '   "maxiter" - MAXimal number of ITERations';
        '   "fit_amps" - fit amplitudes (1) or not(0)';
        'All other parameters are different for each line';
        '   amps(i) - amplitude of line';
        '   x0(i)   - position of the line';
        '   fwhh(i) - full weight at half height';
        '   distortion(i) - ';
        '             dist..*gaussian + (1-dist..)*lorentzian';
        '';
        'parameters x0 and fwhh are always fitting variables'
    }, 'About', 'help');
        
% end

