function varargout = fittingbox(varargin)
% POWDERBOX Application M-file for powderbox.fig
%    Intended to use only with 'kazan' viewer

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 14-Dec-2003 16:23:19
% Boris Epel 5-dec-2003, MPI
% boep 21-dec-2003, MPI

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    handles.hh = 0;
    handles.ploter = 0;
    handles.num = 1;
    handles.div = 0;
    handles.Data = [];
    handles.Script{1} = get(handles.list, 'String');
    handles.sim.y = [];
    handles.loadtype = 'onebyone'; % or 'allinone'
    handles.Color =   [1, 0, 0;...
            1, 1, 0;...
            0, 1, 0;...
            0, 0, 0;...
            1, 0, 1;...
            0, 1, 1;...
            1, 1, 0];
    handles.script = {};
    
    %   set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
    guidata(fig, handles);
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
[path,name,ext,ver] = fileparts(hand.fname);
file_str = name;
sizey = size(hand.src.y, 2);

button = questdlg('How to load?',...
    '???','One by One','All in One','Cancel','One by One');        

popString = {};
switch button
    case 'One by One'
        if sizey > 10, 
            button = questdlg('Do you realy want to load data one by one?',...
                ['there are ' num2str(sizey), ' columns'],...
                'Oh, ja ja','May be later','May be later');
            if strcmp(button, 'May be later'), return; end
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
set(handles.popSel, 'String', popString, 'Value', 1);
set(handles.list, 'String', handles.Script{1}, 'Value', 1);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = ePars_Callback(h, eventdata, handles, varargin)
% % Stub for Callback of the uicontrol handles.ePars.
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');
newval = trim(get(handles.ePars, 'String'));

[name, str1] = strtok(Res{val}, '=');
name = trim(name);

switch get(handles.ePars, 'UserData')
    case 'function'
        set(handles.ePars, 'String', newval);
        tokens = fsplitter([name ' = ' newval]);
        % load script 
        for k=1:size(Res,1),
            if ~strcmp(trim(strtok(Res{k},'=')),'f')
                try,eval([Res{k} ';']);end
            end
        end
        Res = {[name ' = ' newval]};
        for k = 1: size(tokens, 2)
            try, var = eval(tokens{k}); catch, var = 0; end
            Res{end+1} = [tokens{k},' = ', num2str(var)];
        end
    case 'string'
        set(handles.ePars, 'string', newval);
        Res{val} = [name ' = ''', newval,''''];
    otherwise
        nums = str2num(newval);
        set(handles.ePars, 'String', num2str(nums));
        Res{val} = [name ' = ' num2str(nums)];
end
set(handles.list, 'String', Res);
list_Callback(h, eventdata, handles);
switch handles.loadtype
    case 'onebyone', 
        handles.Script{get(handles.popSel, 'Value')} = Res;
    case 'allinone'
        handles.Script = Res;
end
guidata(handles.figure1, handles);
hand = guidata(handles.hh);

if strcmp(get(handles.ePars, 'UserData'), 'value'),
    evallist(handles);
end

%--------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)

if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;
shift = get(handles.slPars, 'Value');
num = get(handles.popup, 'Value');
str = get(handles.popup, 'String');
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
ax.x = handles.src.ax.x;
ax.y = handles.src.ax.y;
y = handles.src.y;

prj = hand.Prj;

hand.out = {};

hand.out{1}.ax = handles.src.ax;
hand.out{1}.ax.c = [0, 0, 1]; % always blue ?
hand.out{1}.ax.s = 0;
hand.out{1}.y = handles.src.y;

hand.out{2}.ax.x = handles.src.ax.x;
hand.out{2}.ax.y = handles.src.ax.y;    
hand.out{2}.y = handles.sim.y;
hand.out{2}.ax.c = [1, 0, 0]; % always black
hand.out{2}.ax.s = 0;

states = {'off', 'on'};
st = get(handles.mSimDif, 'Checked');
check = find(strcmp(states, st)) -1;
if check, 
    hand.out{3}.ax.x = handles.src.ax.x;
    hand.out{3}.ax.y = handles.src.ax.y;    
    hand.out{3}.y = handles.src.y - handles.sim.y;
    hand.out{3}.ax.c = [0, 1, 0]; % always black
    hand.out{3}.ax.s = 0;
end

guidata(handles.hh, hand);
[fpath,name,ext,ver] = fileparts(get(handles.hh, 'FileName'));
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
set(handles.list, 'String', handles.Script{num}, 'Value', 1);
if strcmp(handles.loadtype, 'onebyone'), 
    hand = guidata(handles.hh);
    num = get(handles.popSel, 'Value');
    
        set(hand.sl2Dplot, 'Value', num);
        kazan('sl2Dplot_Callback',hand.sl2Dplot,[],hand);
end
list_Callback(handles.list, eventdata, handles);

% set(handles.butCol, 'BackgroundColor', handles.Color(num, :)); 

% --------------------------------------------------------------------
function varargout = nearme(X, xi)
% nearme(X, xi)
% k = nearme(X, xi)
% [k, v] = nearme(X, xi)
% 
%   return value (v) and index (k) 
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
% --------------------------------------------------------------------
function varargout = autofitt(handles, x, y, steps)
handles = guidata(handles.figure1);

%%%%%%%% construct fitting function %%%%%%%
strings = get(handles.list, 'String');
val = 0;
num = [];
valstr = {};
for k = 1:size(strings, 1)
    tempstr = strings{k};
    if strcmp(tempstr(1:3), 'f ='), val = k; 
    else
        valstr{end + 1} = strtokstr(tempstr, ' =');
        num(end + 1) = k;     
    end
    clear tempstr;
end
if val == 0, 
    msgbox('can not find function'); 
    return;
end

% recognising variables names
[in, funstr] = convfunc(strings{val}, {strings{num}});

%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%
resstr = funstr(4:end);
strChi = ['sum((', resstr,' - real(y)).^2)/1e-6'];
f = inline(strChi, 'in', 'x', 'y');

opt = [];
if ~isempty(steps), opt = optimset('MaxIter', steps); end
out = fminsearch(f, in, opt, x, y);

in = out;
eval(['Yres = ', resstr, ';']);
for m = 1:length(num)
    eval(['finalstr = [valstr{m},'' = '', num2str(in(', num2str(m),'))];']);
    strings{num(m)} = finalstr;
end
if nargout == 1, varargout{1} = Yres;
else varargout{1} = Yres; varargout{2} = strings;
    return;
end

% splitter of formulae
% --------------------------------------------------------------------
function [tokens] = fsplitter(formulae)
[name, f] = strtok(formulae, '='); f = f(2:end);
idx = f=='+' | f=='-' | f=='/' | f=='*' | f=='^' | f==','| f=='''';
funlist = {'./'; '.*';'.^';'(';')';'exp';'log';'cos';'sin';'abs'; 'lshape'};
for k=1:size(funlist, 1)
    idxexp = strfind(f, funlist{k});
    for l = 1:size(funlist{k}, 2)
        idx(idxexp + l - 1)=1;
    end
end
f(idx) = ' ';
tokens = {};
while 1
    [tok, f] = strtok(trim(f));
    if isempty(tok), break; end;
    % Is token an x axis ?
    if strcmp(tok, 'x'), continue; end;
    % Is token an i or j ?
    if strcmp(tok, 'i') | strcmp(trim(tok), 'j'), continue; end;
    % Is token a number ?
    if isempty(str2num(tok))
        tokens{end+1} = tok;
    end
end

%--------------------------------------------------------------------
function [funpar, funstr] = convfunc(strFunc, cellVar)
% stFunc - string with function, like 'f = n*i*x + u*j*a';
% cellVar - cell array of strings of parameters like {'n = 0'_;_ 'i = 0'}
% strFunc = 'f =n*x + b*a';
% cellVar = {'n = 1'; 'b = 2'; 'a = 5'};
string = strFunc(1:end); % cut 'f ='
string(end+1) = ' ';

funlist = {'lshape'; 'exp';'log';'cos';'sin';'abs'; './'; '.*';'.^'; ...
        '+'; '-'; '/'; '*'; '^'; ','; ''''; '='; '('; ')'};

% AlSi comment: funlist must be sorted by decreasing length(funlist{i})
% otherwise additional "variables" can occure: 'acos' -> 'a   ' (if 'cos' will checked at first)
% now it is not so important, but in future ....

tempstr = string;
for j = 1:length(funlist)
    pos = findstr(tempstr, funlist{j});
    if isempty(pos), continue; end
    for i = 1: length(pos)
        len = length(funlist(j));
        for q = (1:len)
            tempstr(q + pos-1) = ' ';
        end
    end
end
shiftpos = 0;
for k = 1:length(cellVar)
    varstr = ['in(' num2str(k), ')'];
    varname = strtokstr(cellVar{k}, ' ='); % name of variable
    varlen = length(varname);
    
    varpos = findstr(tempstr, [' ', varname, ' ']);
    
    stpos = varpos(1)+shiftpos;
    endpos = varpos(1)+shiftpos + varlen + 1;
    
    str1 = string(1:stpos);
    str2 = string(endpos:end);
    string = [str1, varstr, str2];
    
    eval([cellVar{k}, ';']);
    eval([varstr, '=', varname, ';']);
    
    shiftpos = length(varstr) - varlen + shiftpos;
end
funpar = in;
funstr = string;

%--------------------------------------------------------------------
function butFitt_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
sliderval = get(hand.sl2Dplot, 'Value');
prj = strcmp({'X', 'Y'}, hand.Prj);

switch h
    case handles.butFitt
        steps = [];    
    case handles.but1Step
        steps = 1;
    otherwise
        steps = [];
end

eval(['x  = handles.src.ax.', lower(hand.Prj), ';']);
eval(['y = handles.src.y', hand.sel, ';']);

[yres strings]= autofitt(handles, x, y, steps);
eval(['handles.sim.y', hand.sel,'= yres;']);
set(handles.list, 'String', strings, 'Value', 1);
list_Callback(h, [], handles);
if strcmp(handles.loadtype, 'onebyone'), 
    num = get(handles.popSel, 'Value');
    handles.Script{num} = strings;
end
guidata(handles.figure1, handles);
plothis(handles);

%--------------------------------------------------------------------
function butFittAll_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
sliderval = get(hand.sl2Dplot, 'Value');
prj = find(strcmp({'Y', 'X'}, hand.Prj));

num = size(handles.src.y, prj);
switch prj
    case 2
        str = '(:, i)';
    case 1
        str = '(i, :)';
end

eval(['x  = handles.src.ax.', lower(hand.Prj), ';']);

switch handles.loadtype
    case 'onebyone'
        for i = 1:num
            set(handles.popSel, 'Value', i);
            set(handles.list, 'String', handles.Script{i});
            eval(['y = handles.src.y', str, ';']);
            [yres script] = autofitt(handles, x, y, []);
            eval(['handles.sim.y', str,'= yres;']);
            set(handles.list, 'String', script);
            handles.Script{i} = script;
        end
    case 'allinone'
        %%%%%%%% construct fitting function %%%%%%%
        strings = get(handles.list, 'String');
        val = 0;
        num = [];
        valstr = {};
        for k = 1:size(strings, 1)
            tempstr = strings{k};
            if strcmp(tempstr(1:3), 'f ='), val = k; 
            else
                valstr{end + 1} = strtokstr(tempstr, ' =');
                num(end + 1) = k;     
            end
            clear tempstr;
        end
        if val == 0, 
            msgbox('can not find function'); 
            return;
        end
        % recognising variables names
        [in, funstr] = convfunc(strings{val}, {strings{num}});
        % fitting
        for i = 1:num
            eval(['y = handles.src.y', str, ';']);
            resstr = funstr(4:end);
            strChi = ['sum((', resstr,' - real(y)).^2)/1e-6'];
            f = inline(strChi, 'in', 'x', 'y');
            out = fminsearch(f, in, [], x, y);
            in = out;
            eval(['handles.sim.y', str,' = ', resstr, ';']);
        end
        for m = 1:length(num)
            eval(['strings{num(m)}  = [valstr{m},'' = '', num2str(in(', num2str(m),'))];']);
        end
        set(handles.list, 'String', strings);
end

guidata(handles.figure1, handles);
plothis(handles);

set(handles.popSel, 'Value', 1);
popSel_Callback(handles.popSel, [], handles);

%--------------------------------------------------------------------
function evallist(handles)
str = get(handles.list, 'String');
val = 0;
hand = guidata(handles.hh);
eval(['x = handles.src.ax.', lower(hand.Prj), ';']);
for count = 1:size(str, 1)
    if ~strcmp(str{count}(1:3), 'f =')
        eval([str{count}, ';']);
    else
        val = count;
    end
end
if ~val, disp('EVALLIST: Can not find function'); return; end
eval([str{val}, ';']);
eval(['handles.sim.y', hand.sel, ' = f;']);
guidata(handles.figure1, handles);
plothis(handles);
%--------------------------------------------------------------------
function figure1_close(h, eventdata, handles)
try
    hand = guidata(handles.hh);
    hand.rx = {};
    hand.ry = {};
    hand.rc = {};
    hand.rs = {};
    guidata(handles.hh, hand);
    eval([name '(''SetDataSource'', 1, hand)']);
end
delete(handles.figure1);

%--------------------------------------------------------------------
function mSimDif_Callback(h, eventdata, handles)
states = {'off', 'on'};
st = get(handles.mSimDif, 'Checked');
check = find(strcmp(states, st))-1;

set(handles.mSimDif, 'Checked', states{~check + 1});
guidata(handles.figure1, handles);
plothis(handles);