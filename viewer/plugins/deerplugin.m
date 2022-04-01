function varargout = deerplugin(varargin)
if nargin<1
    
    token = 'DEERPLUGIN_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,n,e] = fileparts(which('kazan'));
    inifilename = [fpath, filesep, 'kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try opc = ini.KazanViewer.MultyKazan; end
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end
% 
% 	fig = openfig(mfilename,'new');
%     setappdata(0, token, [oldfig, fig]);
    
    fig = figure; %openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    set(fig, 'Units', 'pixels', 'Name', 'DEER', 'NumberTitle', 'off');
    set(0, 'Units', 'pixels');
    vv = get( 0, 'ScreenSize');
    set(fig, 'Position', [319 100 250 min(800, vv(4)-200)], 'Tag', 'MainFigure', 'Visible', 'off', 'MenuBar', 'none'); %%% hide until everything is done
    set(fig, 'ResizeFcn', 'deerplugin(''resizeFigure'', guidata(gcf))', 'Units', 'pixels');
     
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    
    handles.hh = 0;
    handles.figStruct = guiDB;
    handles.BL{1}.type = 'b*exp(-a*t)+c';
    handles.BL{1}.a = 4;
    handles.BL{1}.b = 1e4;
    handles.BL{1}.c = 100;
    handles.BL{1}.x0_ns = 0;
    handles.BL{1}.x1_ns = 100;
    handles.BL{1}.x2_ns = 100;
    
    handles.BL{2}.type = 'loc_polyn(y, t, np, x1_ns, x2_ns)';
    handles.BL{2}.np = 2;
    handles.BL{2}.x0_ns = 0;
    handles.BL{2}.x1_ns = 100;
    handles.BL{2}.x2_ns = 100;
    
    handles.BL{3}.type = ''; % none
%     handles.BL{3}.x0_ns = 0;
%     handles.BL{3}.x1_ns = 100;
%     handles.BL{3}.x2_ns = 100;
%     
    handles.Pake{1}.type = 'C0*exp(-(r-r0).^2/dr^2/4)';
    handles.Pake{1}.r0 = 2.5;
    handles.Pake{1}.dr = 0.2;
    handles.Pake{1}.C0 = 1;
    handles.Pake{1}.ShiftY = 0;
    handles.Pake{1}.Amp = 1;
    handles.PakeFit{1} = [0, 0, 0, 1, 1];
    
    handles.Pake{2}.type = 'C01*exp(-(r-r01).^2/dr1^2/4)+C02*exp(-(r-r02).^2/dr2^2/4)';
    handles.Pake{2}.r01 = 2.5;
    handles.Pake{2}.dr1 = 0.2;
    handles.Pake{2}.C01 = 1;
    handles.Pake{2}.r02 = 4.5;
    handles.Pake{2}.dr2 = 0.2;
    handles.Pake{2}.C02 = 1;
    handles.Pake{2}.ShiftY = 0;
    handles.Pake{2}.Amp = 1;
    handles.PakeFit{2} = [0, 0, 0, 0, 0, 0, 0, 1, 1];
    
    load('pake_base40.mat'); % Pake(r, t), t ns, r nm
    handles.kern = base;
    handles.kern_r = r;
    handles.kern_t = t;
    
    handles.src = [];
    
    handles = createGUI(handles, handles.figStruct);
    guidata(handles.MainFigure, handles);
    set(fig, 'Visible', 'on');
    resizeFigure(handles);
    pmFitMethod_Callback(handles);
    pmBGFitMethod_Callback(handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        %     [varargout{1:nargout}] = 
        feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
%         dispstatus(gcf, 2); % Error
    end
end

function Struct = guiDB

Struct{1} = struct('Tag', 'pbLoad', 'String', 'Load File', 'Enable','on', 'Style', 'pushbutton',...
                   'Callback', 'deerplugin(''LoadFile_Callback'', guidata(gcf))' , ...
                   'Size', [-1, 30], 'Units', 'pixels');
Struct{1}.Children = {};  
% ---------------------------------------------- %
strsub{1} = struct('Tag', 'pmFitMethod', ...
                   'Style', 'popupmenu', ...
                   'Callback', 'deerplugin(''pmFitMethod_Callback'', guidata(gcf))' , ...
                   'Size', [-1, 30], 'Units', 'pixels');
strsub{1}.String = {'Single Gaussian', 'Double Gaussian'};          
strsub{1}.Children = {}; 

strsub{2} = struct('Tag', 'lbFitParams', 'String', {'nothing loaded yet'}, ...
                   'Style', 'listbox', 'BackgroundColor', [1, 1, 1], ...
                   'Callback', 'deerplugin(''lbFitParams_Callback'', guidata(gcf))' ,  ...
                   'Size', [-1, -1], 'Units', 'pixels');     
strsub{2}.Children = {};   

strsub{3} = struct('Tag', 'spEditParams', 'Value', 0, ...
                   'Style', 'spinbox', ...
                   'Callback', 'deerplugin(''spEditParams_Callback'', guidata(gcf))' ,  ...
                   'Size', [-1, 30], 'Units', 'pixels');    
strsub{3}.Children = {};  

% >>>>>>
grsub{1} = struct('Tag', 'edFitNIter', 'String', '300', ...
                   'Style', 'edit', ...
                   'Callback', 'deerplugin(''edFitNIter_Callback'', guidata(gcf))' ,  ...
                   'Size', [-1, 30], 'Units', 'pixels');  
grsub{1}.Children = {};  

grsub{2} = struct('Tag', 'pbFit', 'String', 'Fit', ...
                   'Style', 'pushbutton', ...
                   'Callback', 'deerplugin(''pbFit_Callback'', guidata(gcf))' ,...
                   'Size', [-1, 30], 'Units', 'pixels');    
grsub{2}.Children = {}; 

strsub{4} = struct('Tag', 'fr1', 'Value', 0, ...
                   'Style', 'subgrid', ...
                   'Callback', '' ,  ...
                   'Size', [-1, 40], 'ChildPosition', 'H', 'Units', 'pixels');        
strsub{4}.Children = grsub;                
%<<<<<    

Struct{2} = struct('Tag', 'panFit', 'Title', 'Pake fit', 'Style', 'panel', ...
                     'ChildPosition', 'V', 'Size', [-1, -3], 'Units', 'pixels');
Struct{2}.Children=strsub;                
% ----------------------------------------------
strsub1{1} = struct('Tag', 'chkShowBG', 'String', 'Show Background', ...
                   'Style', 'checkbox', ...
                   'Callback', 'deerplugin(''chkShowBG_Callback'', guidata(gcf))' , ...
                   'Size', [-1, 30], 'Units', 'pixels');
strsub1{1}.Children = {};  

strsub1{2} = struct('Tag', 'pmBGFitMethod',  ...
                   'Style', 'popupmenu', ...
                   'Callback', 'deerplugin(''pmBGFitMethod_Callback'', guidata(gcf))' ,  ...
                   'Size', [-1, 30], 'Units', 'pixels');
strsub1{2}.String = {'exponential', 'polynom', 'none'};          
strsub1{2}.Children = {};     

strsub1{3} = struct('Tag', 'lbBGFitParams', 'String', {'nothing loaded yet'}, ...
                   'Style', 'listbox','BackgroundColor', [1, 1, 1], ...
                   'Callback', 'deerplugin(''lbBGFitParams_Callback'', guidata(gcf))' , ...
                   'Size', [-1, -1], 'Units', 'pixels');           
strsub1{3}.Children = {};

strsub1{4} = struct('Tag', 'spBGEditParams', 'Value', 0, ...
                   'Style', 'spinbox', ...
                   'Callback', 'deerplugin(''spBGEditParams_Callback'', guidata(gcf))' , ...
                   'Size', [-1, 30], 'Units', 'pixels');  
strsub1{4}.Children = {};      

% >>>>>>
grsub1{1} = struct('Tag', 'edBGFitNIter', 'String', '300', ...
                   'Style', 'edit', ...
                   'Callback', 'deerplugin(''edBGFitNIter_Callback'', guidata(gcf))' , ...
                   'Size', [-1, 30], 'Units', 'pixels');  
grsub1{1}.Children = {};   

grsub1{2} = struct('Tag', 'pbBGFit', 'String', 'Fit', ...
                   'Style', 'pushbutton', ...
                   'Callback', 'deerplugin(''pbBGFit_Callback'', guidata(gcf))' ,  ...
                   'Size', [-1, 30], 'Units', 'pixels'); 
grsub1{2}.Children = {};  

strsub1{5} = struct('Tag', 'fr2', 'Value', 0, ...
                   'Style', 'subgrid', ...
                   'Callback', '' ,  ...
                   'Size', [-1, 40], 'ChildPosition', 'H', 'Units', 'pixels');  
strsub1{5}.Children=grsub1;
% <<<<<
Struct{3} = struct('Tag', 'panBG', 'Title', 'Baseline', 'Style', 'panel', ...
                    'ChildPosition', 'V', 'Size', [-1, -2], 'Units', 'pixels');
Struct{3}.Children=strsub1;

function handles = createGUI(handles, HStruct)     
for ii = 1:length(HStruct)
    handles = createObj(handles, HStruct{ii});
end


function handles = createObj(handles, Obj)
    parent = handles.MainFigure;
    tObj = Obj;
    tObj = rmfield(tObj, 'Children');
    tObj = rmfield(tObj, 'Size');
    if isfield(tObj, 'ChildPosition')
        tObj = rmfield(tObj, 'ChildPosition');
    end
    switch tObj.Style
        case {'edit', 'listbox', 'pushbutton', 'checkbox', 'popupmenu'}
            hh = uicontrol(parent, tObj);
        case 'spinbox' 
            if isfield(handles, 'uiCM')
                lastMenu = length(handles.uiCM)+1;
            else
                lastMenu = 1;
            end
            handles.uiCM{lastMenu} = uicontextmenu;
            for ii = 1:17
                cb = sprintf('h=guidata(gcf);set(h.%s,''Step'',1e%d);gg=get(h.uiCM{%d}, ''Children'');set(gg, ''Checked'', ''off'');set(gcbo, ''Checked'', ''on'');',tObj.Tag,ii-8,lastMenu);
                uimenu(handles.uiCM{lastMenu}, 'Label', sprintf('1e%d', ii-8), 'Callback', cb);
            end
            tObj.UIContextMenu = handles.uiCM{lastMenu};
            hh = spinbox(parent, tObj);
        case 'panel'
            tObj = rmfield(tObj, 'Style');
            hh = uipanel(parent, tObj);
            
        case 'subgrid'
            hh = [];% do nothing
    end

    if ~isempty(hh)
        handles = setfield(handles, Obj.Tag, hh);
    end
    
    if ~isempty(Obj.Children)
        handles = createGUI(handles,  Obj.Children);
    end

function resizeFigure(handles)
handles = guidata(handles.MainFigure);

parPos = get(handles.MainFigure, 'Position');
posStr = getPositions({}, [0, 0, parPos(3:4)], 'V', handles.figStruct);

for ii = 1:length(posStr)
    hh = getfield(handles, posStr{ii}.Tag);
    set(hh, 'Position', posStr{ii}.Position);
end

function posStr = getPositions(posStr, parPos, arr, HStruct)

space = 5;

if strcmp(arr, 'V')
    ntofit = 2;
    stretch = 1;
    range = length(HStruct):-1:1;
else 
    ntofit = 1;
    stretch = 2;
    range = 1:length(HStruct);
end

% find all
fixed = 0;
portions = 0;
for ii = 1:length(HStruct)
    ss = HStruct{ii}.Size(ntofit);
    if ss>0, fixed = fixed+ss;
    elseif ss<0
        portions = portions-ss;
    else
        disp('Size == 0, not possible');
    end        
end

Val = parPos(ntofit+2)-fixed-space*(length(HStruct)+3);

dV = Val/portions;

iniVar = parPos(ntofit);
for ii = range
    
    pos(stretch) = parPos(stretch)+space;
    if HStruct{ii}.Size(stretch)<0
        pos(stretch+2) = parPos(stretch+2)-2*space;
    else
        pos(stretch+2) = HStruct{ii}.Size(stretch);
    end
    
    pos(ntofit) = iniVar+space;
    
    ss = HStruct{ii}.Size(ntofit);
    if ss>0, 
        pos(ntofit+2) = ss;
    elseif ss<0
        pos(ntofit+2) = dV*(-ss);
    else
        disp('Size == 0, not possible');
    end    
    iniVar = pos(ntofit)+pos(ntofit+2);
    if ~strcmp(HStruct{ii}.Style, 'subgrid')
        
        posStr{end+1}.Position = pos;
        posStr{end}.Tag = HStruct{ii}.Tag;
    end
    if ~isempty(HStruct{ii}.Children)
        posStr = getPositions(posStr, pos, HStruct{ii}.ChildPosition, HStruct{ii}.Children);
    end
end

function pmFitMethod_Callback(handles)
num = get(handles.pmFitMethod, 'Value');
Str = handles.Pake{num};
fitPars = handles.PakeFit{num};
ListStr = getList(Str, fitPars);
set(handles.lbFitParams, 'String', ListStr, 'Value', 1);
lbFitParams_Callback(handles);
plotThings(handles);

function pmBGFitMethod_Callback(handles)
num = get(handles.pmBGFitMethod, 'Value');
Str = handles.BL{num};
ListStr = getList(Str);
set(handles.lbBGFitParams, 'String', ListStr, 'Value', 1);
lbBGFitParams_Callback(handles);
plotThings(handles);

function ListStr = getList(Struc, varargin)

names = fieldnames(Struc);
if nargin>1
    cols = varargin{1};
else
    cols = zeros(length(names), 1);
end
ListStr = {};

str{1} = '';
str{2} = '<HTML><FONT color="red">';
str{3} = '<HTML><FONT color="green">';

for ii = 1:length(names)
    val = getfield(Struc, names{ii});
    if ischar(val)
        if strcmp(names{ii}, 'type')
            continue;
        else
            ListStr{end+1} = sprintf('%s %s = %s', str{cols(ii-1)+1}, names{ii}, val);
        end
    end
    if length(val)==1
         ListStr{end+1} = sprintf('%s %s = %f', str{cols(ii-1)+1}, names{ii}, val);
    else
        for jj = 1:length(val)
            ListStr{end+1} = sprintf('%s %s(%d) = %f', str{cols(ii-1)+1}, names{ii}, jj, val(jj));
        end
    end
end

function lbFitParams_Callback(handles)
% % Stub for Callback of the uicontrol handles.listbox1.
Res = get(handles.lbFitParams, 'String');
val = get(handles.lbFitParams, 'Value');

[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(2:end)), ' '];

if str2(1)==''''
    primes = findstr(str2, '''');
    if size(primes, 2) < 2, primes(2) = size(str2, 2); end
    set(handles.spEditParams, 'String', str2(2:(primes(2)-1)));
    set(handles.spEditParams, 'UserData', 'string');
else
    dig = str2num(str2);
    set(handles.spEditParams, 'String', num2str(dig, 6));
    set(handles.spEditParams, 'UserData', 'value');
end

function lbBGFitParams_Callback(handles)
% % Stub for Callback of the uicontrol handles.listbox1.
Res = get(handles.lbBGFitParams, 'String');
val = get(handles.lbBGFitParams, 'Value');

[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(2:end)), ' '];

if str2(1)==''''
    primes = findstr(str2, '''');
    if size(primes, 2) < 2, primes(2) = size(str2, 2); end
    set(handles.spBGEditParams, 'String', str2(2:(primes(2)-1)));
    set(handles.spBGEditParams, 'UserData', 'string');
else
    dig = str2num(str2);
    set(handles.spBGEditParams, 'String', num2str(dig, 6));
    set(handles.spBGEditParams, 'UserData', 'value');
end

function spEditParams_Callback(handles)
Res = get(handles.lbFitParams, 'String');
val = get(handles.lbFitParams, 'Value');
newval = trim(get(handles.spEditParams, 'String'));
PType = get(handles.pmFitMethod, 'Value');

[name, str1] = strtok(Res{val}, '=');
%name = trim(name); opt_mode=name(1:2); name=name(3:end);
commOK = 1;
switch get(handles.spEditParams, 'UserData')
case 'string'
    set(handles.spEditParams, 'String', newval);
    Res{val} = [name ' = ''', newval,''''];
    handles.Pake{PType}  = setfield(handles.Pake{PType}, name, newval);
otherwise
    nums = str2num(newval);
    set(handles.spEditParams, 'String', trim(newval));
    Res{val} = [name ' = ' trim(newval)];
    vv = get(handles.spEditParams, 'Value');
    handles.Pake{PType}  = setfield(handles.Pake{PType}, name, vv);
end

set(handles.lbFitParams, 'String', Res);
guidata(handles.MainFigure, handles);
plotThings(handles);
% list_Callback(h, eventdata, handles);

function spBGEditParams_Callback(handles)
Res = get(handles.lbBGFitParams, 'String');
val = get(handles.lbBGFitParams, 'Value');
newval = trim(get(handles.spBGEditParams, 'String'));
BGType = get(handles.pmBGFitMethod, 'Value');

[name, str1] = strtok(Res{val}, '=');
name = trim(name);
commOK = 1;
switch get(handles.spBGEditParams, 'UserData')
case 'string'
    set(handles.spBGEditParams, 'String', newval);
    Res{val} = [name ' = ''', newval,''''];
    handles.BL{BGType} = setfield(handles.BL{BGType}, name, newval);
otherwise
    nums = str2num(newval);
    set(handles.spBGEditParams, 'String', trim(newval));
    Res{val} = [name ' = ' trim(newval)];
    vv = get(handles.spBGEditParams, 'Value');
    handles.BL{BGType} = setfield(handles.BL{BGType}, name, vv);
end
set(handles.lbBGFitParams, 'String', Res);

guidata(handles.MainFigure, handles);
plotThings(handles);
% list_Callback(h, eventdata, handles);

function LoadFile_Callback(handles)
hand = guidata(handles.hh);
if isempty(hand),  
    msgbox('Can not find main box handles', 'Error');
    return;
end
[path,name,ext] = fileparts(hand.fname);
file_str = name;
set(handles.pbLoad, 'String', file_str, 'BackgroundColor', [0.5, 0.8, 0.5]);
handles.sim.y = hand.src.y*0;
handles.src = hand.src;

guidata(handles.MainFigure, handles);
plotThings(handles);

function plotThings(handles)
if isempty(handles.src), return; end

src = handles.src;
hand = guidata(handles.hh);
% [prj, selec] = getslider_kazan(hand);
hand.out = {};
% ax = hand.src.ax;
% y = hand.src.y;

% if handles.process, [y, ax] = processing(y, ax);
% else
%     ax.diff = '';
%     ax.filt = '';
% end

showBG = get(handles.chkShowBG, 'Value');
PakeModel = get(handles.pmFitMethod, 'Value');
BGType = get(handles.pmBGFitMethod, 'Value');

%%%%%% get BG
type = handles.BL{BGType}.type;
listStr = get(handles.lbBGFitParams, 'String');
outBG = src;
outBG.y = src.y*0;
idx1 = 1; idx2 = length(src.y);
x0_ns = 0;
if ~isempty(type)
    
    for ii = 1:length(listStr)
        eval([listStr{ii}, ';']);
    end
    
    [~, idx1] = min(abs(src.ax.x-x1_ns));
    [~, idx2] = min(abs(src.ax.x-x2_ns));
    
    src.ax.x = src.ax.x-x0_ns;
    outBG.ax.x = src.ax.x;
    
    t = src.ax.x;
    y = real(src.y);
    eval(['outBG.y =', type, ';']);
end

DEERy = calcDEER(handles);
% BG = getEval(handles.src, type, get(handles.lbBGFitParams, 'String'));
if showBG
    % Data
    hand.out{1} = src;
    hand.out{1}.ax = src.ax;
    hand.out{1}.ax.Color = [0 0 1];
    hand.out{1}.ax.s = 0;
%     hand.out{1}.y = y;
    hand.out{1}.title = 'Src';
     % BG
    hand.out{end+1} = src;
    hand.out{end}.ax.x = src.ax.x(idx1:idx2);
    hand.out{end}.ax.Color = [0 0 1];
    hand.out{end}.ax.Marker = '.';
    hand.out{end}.ax.s = 0;
    hand.out{end}.y = src.y(idx1:idx2, :);
    hand.out{end}.title = 'Src_forBG';
%%%%%%%%% BG
    hand.out{end+1} = outBG;
    hand.out{end}.ax.Color = [1 0 0]; % always gray ?
    hand.out{end}.ax.s = 0;
    hand.out{end}.y = outBG.y;
    hand.out{end}.title = 'BG';
    
    hand.out{end+1}.ax = src.ax;
    hand.out{end}.ax.Color = [1 0 0];
    hand.out{end}.ax.s = 0;
    hand.out{end}.y = DEERy+outBG.y;
    hand.out{end}.title = 'DEER+BG';
else
    hand.out{1}.ax = src.ax;
    hand.out{1}.ax.Color = [0 0 1];
    hand.out{1}.ax.s = 0;
    hand.out{1}.y = src.y-outBG.y;
    hand.out{1}.title = 'Src';
    
    
    hand.out{end+1}.ax = src.ax;
    hand.out{end}.ax.Color = [1 0 0]; 
    hand.out{end}.ax.s = 0;
    hand.out{end}.y = DEERy;
    hand.out{end}.title = 'DEER';
    
end


%%%%%% send things to Kazan Viewer for plot
guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

function out = getBGEval(src, type, listStr)
out = src;
if isempty(type)
    out.y = out.y*0;
    return;
end
for ii = 1:length(listStr)
    eval([listStr{ii}, ';']);
end
idx1 = min(abs(t-x1_ns));
idx2 = min(abs(t-x2_ns));
t = src.ax.x;
eval(['out.y =', type, ';']);

function chkShowBG_Callback(handles)
plotThings(handles);

function bl = loc_polyn(y, t, np, x1_ns, x2_ns)
    [~, idx1] = min(abs(t-x1_ns));
    [~, idx2] = min(abs(t-x2_ns));
pp = polyfit(t(idx1:idx2), y(idx1:idx2), np);
bl = polyval(pp, t);

function [yy, Pars] = doFitting(Struct, Data, whattofit)
ss = rmfield(Struct, 'type');
names = fieldnames(ss);
funct = Struct.type;

% strfit = {'%f', 'inp(%i)'};
nVars = 0;
funStr = ['(', funct, ')'];
iniInp = [];
for ii= 1:(length(names)-2)
    %     if ~(strcmp(names, 'ShiftY')|strcmp(names, 'Amp'))
    idx1 = findstr(funStr, names{ii})-1;
    idx2 = idx1+length(names{ii})+2;
    
    if whattofit(ii)
        funStr = sprintf([funStr(1:idx1), 'inp(%i)', funStr(idx2:end)], ii);
        nVars = nVars+1;
        iniInp(end+1) = getfield(Struct, names{ii});
    else
        funStr = sprintf([funStr(1:idx1), '%f', funStr(idx2:end)], getfield(Struct, names{ii}));
    end
    
%     end
end


% eval(['distr =', handles.Pake{PakeModel}.type, ';']);
% distr = 0.01*distr/sum(distr);

% ['y = Amp*(exp(', funStr, '*kernel)-1)+ShiftY;

% eval(funStr);
f = inline(funStr, 'inp', 'x', 'y');

res = fminsearch(fitfun, iniInp);

yy = fitfun(res);

function y = calcDEER(handles)

%     load('Pake_arrays.mat'); % Pake(r, t), t ns, r nm
% handles.kern = Pake;
% handles.kern_r = r;
% handles.kern_t = t;

BGType = get(handles.pmBGFitMethod, 'Value');
x0 = handles.BL{BGType}.x0_ns;

src = handles.src;
x = src.ax.x-x0;
xmax = max(x); %
xmin = min(x);
[~, idx2] = min( abs(handles.kern_t-xmax)); idx2 = min(length(handles.kern_t), idx2+1);
[~, idx1] = min( abs(handles.kern_t-xmin)); idx1 = max(1, idx1-1);
kernel=handles.kern(:, idx1:idx2)-ones(size(handles.kern(:, idx1:idx2)));
kernel = interp1(handles.kern_t(idx1:idx2), kernel.',  x);
kernel = kernel.';
lastx = x;

PakeModel = get(handles.pmFitMethod, 'Value');
listStr = get(handles.lbFitParams, 'String');
for ii = 1:length(listStr)
    eval([listStr{ii}, ';']);
end
r = handles.kern_r;
eval(['distr =', handles.Pake{PakeModel}.type, ';']);
distr = 0.01*distr/sum(distr);

y = Amp*(exp(distr.'*kernel)-1)+ShiftY;
y = y(:);


function pbFit_Callback(handles)
num = get(handles.pmFitMethod, 'Value');
Str = handles.Pake{num};
fitPars = handles.PakeFit{num};

[Res, Pars] = doFitting(firPars, Data, handles.PakeFit);

ListStr = getList(Str, fitPars);
set(handles.lbFitParams, 'String', ListStr, 'Value', 1);
lbFitParams_Callback(handles);
plotThings(handles);

function pbBGFit_Callback(handles)
    