function varargout = fftplugin(varargin)
% POWDERBOX Application M-file for powderbox.fig
%    Intended to use only with 'kazan' viewer

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.0 21-Dec-2004 13:28:22
% Boris Epel 5-dec-2003, MPI
% boep 21-dec-2003, MPI

if nargin == 0  % LAUNCH GUI
    token = 'FFTPLUGIN_open';
    oldfig = getappdata(0, token)
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,~,e] = fileparts(which('kazan'));
    inifilename = [fpath filesep 'kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try opc = ini.KazanViewer.MultyKazan; end
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end

	fig = figure; %openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    set(fig, 'Units', 'pixels', 'Name', 'FFT plugin', 'NumberTitle', 'off');
    set(0, 'Units', 'pixels');
    vv = get( 0, 'ScreenSize');
    set(fig, 'Position', [319 100 250 min(800, vv(4)-200)], 'Tag', 'figure1', 'MenuBar', 'none'); %%% hide until everything is done
    set(fig, 'ResizeFcn', 'fftplugin(''figure1_ResizeFcn'',[],[],guidata(gcf))');
    
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  handles = loadGUI(handles);
  
  handles.hh = 0;
  handles.ploter = 0;
  handles.num = 1;
  handles.div = 0;
  handles.Spec = {};
  handles.sim{1}.Script = get(handles.list, 'String');
  handles.sim{1}.x = [];
  handles.sim{1}.x1 = [];
  handles.sim{1}.y = [];
  handles.sim{1}.stype = 'sim';
  handles.Color =   [1, 0, 0;...
      1, 1, 0;...
      0, 1, 0;...
      0, 0, 0;...
      1, 0, 1;...
      0, 1, 1;...
      1, 1, 0];
  handles.SumCol = [1, 0, 0];
  handles.script = {};

  set(handles.LoSc, 'Accelerator', 'L');
  set(handles.SvSc, 'Accelerator', 'S');
  set(handles.mCalc, 'Accelerator', 'C');
  set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
  guidata(fig, handles);
  list_Callback([], [], handles);
  
  if nargout > 0
    varargout{1} = fig;
  end
  
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  
  try
%     
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
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
function FileUpdate(varargin);
return;

function handles = loadGUI(handles)

handles.File = uimenu(handles.figure1, 'Tag', '', 'Label', 'File', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.Untitled_1 = uimenu(handles.File, 'Tag', 'Untitled_1', 'Label', 'Add baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', '', 'Accelerator', '' ); 
handles.mFileDelete = uimenu(handles.File, 'Tag', 'mFileDelete', 'Label', 'Delete', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''butDel_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.LoSc = uimenu(handles.File, 'Tag', 'LoSc', 'Label', 'Load', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''LoSc_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'L' ); 
handles.SvSc = uimenu(handles.File, 'Tag', 'SvSc', 'Label', 'Save', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''SvSc_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'S' ); 
handles.mFileCopy = uimenu(handles.File, 'Tag', 'mFileCopy', 'Label', 'Copy', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mFileCopy_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mFileScript = uimenu(handles.File, 'Tag', 'mFileScript', 'Label', 'Show script', 'Checked', 'off', 'Separator', 'on', 'Callback', 'fftplugin(''mFileScript_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mFileExecScript = uimenu(handles.File, 'Tag', 'mFileExecScript', 'Label', 'Execute script', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mFileScript_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 

handles.Options = uimenu(handles.figure1, 'Tag', '', 'Label', 'Options', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
handles.mSimSubBl = uimenu(handles.Options, 'Tag', 'mSimSubBl', 'Label', 'Subtract baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mSimShow_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
handles.mCalc = uimenu(handles.Options, 'Tag', 'mCalc', 'Label', 'Calculate', 'Checked', 'off', 'Separator', 'on', 'Callback', 'fftplugin(''check_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'C' ); 

handles.mAbout = uimenu(handles.figure1, 'Tag', 'mAbout', 'Label', 'About', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mAbout_Callback'',gcbo,[],guidata(gcbo))' ); 

handles.frame1 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame1', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 
% handles.frame3 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame3', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', ''); 

handles.butFFT = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'butFFT', 'String', 'FFT', 'Value', 0, 'Units', 'pixels', 'Callback', 'fftplugin(''butFFT_Callback'',gcbo,[],guidata(gcbo))'); 
handles.butCalc = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'butCalc', 'String', 'Calculate', 'Value', 0, 'Units', 'pixels', 'Callback', 'fftplugin(''butCalc_Callback'',gcbo,[],guidata(gcbo))'); 
handles.butCol = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butCol', 'String', '', 'Value', 0, 'Units', 'pixels', 'Callback', 'fftplugin(''butCol_Callback'',gcbo,[],guidata(gcbo))'); 
handles.popSel = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popSel', 'String', '1', 'Value', 1, 'Units', 'pixels', 'Callback', 'fftplugin(''popSel_Callback'',gcbo,[],guidata(gcbo))'); 
handles.popup = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popup', 'String', ...
    {'1e-6', ...
     '1e-5', ...
     '1e-4', ...
     '1e-3', ...
     '1e-2', ...
     '0.1', ...
     '1', ...
     '10', ...
     '100', ...
     '1e+3', ...
     '1e+4', ...
     '1e+5', ...
     '1e+6'}, 'Value', 6, 'Units', 'pixels', 'Callback', 'fftplugin(''popup_Callback'',[],[], guidata(gcf))'); 

% handles.slPars = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slPars', 'String', '{, ''}', 'Value', 0, 'Units', 'pixels', 'Callback', 'fftplugin(''slPars_Callback'',gcbo,[],guidata(gcbo))'); 
% handles.ePars = uicontrol(handles.figure1, 'Style', 'edit', 'Tag', 'ePars', 'String', 'ham', 'Value', 0, 'Units', 'pixels', 'Callback', 'fftplugin(''ePars_Callback'',gcbo,[],guidata(gcbo))'); 
    handles.ePars = spinbox(handles.figure1, 'Tag', 'ePars', 'Value', 0.5, 'Step', 0.1, ...
        'Units', 'pixels', 'Callback', 'fftplugin(''ePars_Callback'',[],[], guidata(gcf))');
    
handles.list = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'list', 'String', ...
    {'fft.awin = ''ham''', ...
     'fft.awidth = 1', ...
     'fft.aalpha = 0.6', ...
     'fft.ashift = 0', ...
     'fft.zerofill = 1', ...
     'fft.rshift = 0', ...
     'fft.lshift = 0', ...
     'fft.opt = ''real''', ...
     'fft.xscale = 1', ...
     'fft.cta = 100', ...
     'fft.phase0 = 0', ...
     'fft.bsline = ''none'''},...
     'BackgroundColor', [1, 1, 1], ...
     'Value', 1, 'Units', 'pixels', 'Callback', 'fftplugin(''list_Callback'',gcbo,[],guidata(gcbo))'); 

function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)

FigaPos = get(handles.figure1, 'Position');

if FigaPos(3)<130, FigaPos(3) = 130; set(handles.figure1, 'Position', FigaPos); end
if FigaPos(4)<210, FigaPos(4) = 210; set(handles.figure1, 'Position', FigaPos); end

brd = 5;

hhh  = brd;
butH = 30;
butW = 60;
pos = [brd, hhh, FigaPos(3)-2*brd, butH];
set(handles.butCalc, 'Position', pos);

hhh = hhh+pos(4)+brd;

pos = [brd, hhh, FigaPos(3)-2*brd, butH];
set(handles.butFFT, 'Position', pos);

hhh = hhh+pos(4)+brd*2;
hhF1 = hhh-brd;

pos = [brd*2, hhh, FigaPos(3)-4*brd, butH];
set(handles.ePars, 'Position', pos);

hhh = hhh+pos(4)+brd;
pos = [brd*2, hhh, FigaPos(3)-4*brd, butH];
set(handles.popup, 'Position', pos);

% on the other hand. ...

zzz = FigaPos(4) - brd*2-butH;
pos = [brd*2, zzz, FigaPos(3)-6*brd - butW, butH];
set(handles.popSel, 'Position', pos);

pos = [brd*3+pos(3), zzz, butW, butH];
set(handles.butCol, 'Position', pos);

hhF2 = pos(2)-3*brd;
pos = [brd*2, hhF1+5*brd+2*butH, FigaPos(3)-4*brd, hhF2-hhF1-2*butH-4*brd];
set(handles.list, 'Position', pos);

pos = [brd*1, hhF1, FigaPos(3)-2*brd, hhF2-hhF1+brd];
set(handles.frame1, 'Position', pos);

% --------------------------------------------------------------------
function varargout = popup_Callback(h, eventdata, handles, varargin)
num = get(handles.popup, 'Value');
str = get(handles.popup, 'String');
dim = str2num([str{num}]);
set(handles.ePars, 'Step', dim);

% --------------------------------------------------------------------
function varargout = list_Callback(h, eventdata, handles, varargin)
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');
if size(Res)<2, 
   val = 1; 
end

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
% Stub for Callback of the uicontrol handles.ePars.
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
      delim1 = '';
      delim2 = '';
    if strfind(newval, '[')
        delim1 = '[';
        delim2 = ']';
    end
    nums = str2num(newval);
    set(handles.ePars, 'String', [delim1, num2str(nums), delim2]);
    Res{val} = [name ' = ' num2str(nums)];
end
set(handles.list, 'String', Res);
list_Callback(h, eventdata, handles);
handles.sim{get(handles.popSel, 'Value')}.Script = Res;
guidata(handles.figure1, handles);


if strcmp(get(handles.mCalc, 'Checked'), 'on'), 
   stype = get(handles.popSel, 'Value');
   switch handles.sim{stype}.stype
      case 'sim'
         plothis(handles);
      case 'baseline'
         handles = bscalc(get(handles.popSel, 'Value'), handles);
         plothis(handles);
   end
end

% --------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)
 
if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;
shift = get(handles.slPars, 'Value');
set(handles.slPars, 'Value', 0.5);
num = get(handles.popup, 'Value');
str = get(handles.popup, 'String');
dim = str2num([str{num}]);

CurVal = str2double(get(handles.ePars, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.ePars, 'String', num2str(Val));
ePars_Callback(h, eventdata, handles, varargin);

% --------------------------------------------------------------------
function handles = bscalc(num, handles)

hand = guidata(handles.hh);
x = hand.src.ax.x;
x1 = hand.src.ax.y;
y = hand.src.y;
func = ';';
% load script 
str = get(handles.list, 'String');
for k=1:size(str,1),
  if strcmp(trim(strtok(str{k},'=')),'f')
    func = str{k};
  else
    try,eval([str{k} ';']);end
  end
end
try,eval([func ';']);end

handles.sim{num}.x = x;
handles.sim{num}.x1 = x1;
handles.sim{num}.amp = 1; %abs(max(f)-min(f));
handles.sim{num}.y = f;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function plothis(handles)
handles = guidata(handles.figure1);
colors = handles.Color;

hand = guidata(handles.hh);
ax = hand.src.ax;
y = hand.src.y;

hand.out = {};
    hand.out{1}.ax = ax;
    hand.out{1}.y = y; 
    hand.out{1}.ax.Color = [0, 0, 1]; % always blue
    hand.out{1}.ax.s = 0;
    hand.out{1}.title = 'Src';
% baseline correction
basln = strcmp(get(handles.mSimSubBl, 'Checked'), 'on');

% load script 
fft = [];
str = get(handles.list, 'String');
handles.script={};
for k=1:size(str,1),
    try,eval([str{k} ';']);end
    handles.script{end+1} = [str{k} ';'];
end

bline = zeros(size(y));
for i = 2:handles.num
    if strcmp(handles.sim{i}.stype, 'baseline')
        bline = bline + handles.sim{i}.y;
    end
end
opt = safeget(fft, 'opt', 'real');
if strcmp(opt, 'reim43D'), 
% see below   
    [sx, sy] = size(y);
    n = floor(sy/sx);
    if n-sy/sx ~=0
        error('for reim43D sx ~= N x sy');
    end
    bline = [];
    subspec = safeget(fft, 'subspec', 0);
    if subspec<=0
        reim43Drange = 1:n;
    else
        reim43Drange = subspec;
    end
    
    disp(['reim43D: selected range: [', num2str(reim43Drange),']']);
    ty = [];
    for ii = reim43Drange
        if ii>n
            error(['Specified range of data is wrong. Max N = ', num2str(n)]);
        end
        cutY = y(:, [1:sx] + (ii-1)*sx, :);
        
        switch fft.bsline
            case 'none',
            case 'polyfit1D',
                warning off MATLAB:polyfit:RepeatedPointsOrRescale;
                tbline = polyfit1(ax.x, cytY, n, 1);
            case 'polyfit2D',
                warning off MATLAB:polyfit:RepeatedPointsOrRescale;
                tbline = polyfit2(ax.x, ax.x, cutY, n);
        end
        bline = [bline, tbline];
        ty = [ty, cutY-tbline];
    end
    y = ty;
else
    if isfield(fft,'bsline')
        n = safeget(fft, 'n_pol', 2);
        switch fft.bsline
            case 'none',
            case 'polyfit1D',
                warning off MATLAB:polyfit:RepeatedPointsOrRescale;
                bline = polyfit1(ax.x, y, n, 1) + bline;
            case 'polyfit2D',
                warning off MATLAB:polyfit:RepeatedPointsOrRescale;
                bline = polyfit2(ax.x, ax.y, y, n) + bline;
        end
    end
    y = y - bline;
end




if get(handles.butFFT, 'Value')
    %%%%%%%%%%%%%%%% FFT %%%%%%%%%%%%%%%%%%%%%%%
    handles = guidata(handles.figure1);
    
    fft.fft = 1;
    
    if strcmp(opt, 'reim2'), % 2D FFT: fft(fft(real(y)))
        fft2 = 1;
        fft.opt = 'real';
        [tax, ty] = fftprocessing(ax, y, fft);
        fft.opt = 'imag';
        tax.y = tax.x;
        tax.x = ax.y;
        [rax, ry] = fftprocessing(tax, ty.', fft);
        %%%%%%% alsi %%%%%%%%%%%%% turn back 29.05.07
        tax.y = rax.x;
        tax.x = rax.y;
        ty = ry.';
        rax = tax;
        ry=ty;
        
    elseif strcmp(opt, 'creim2'),
        % try CTA for
        avpar = safeget(fft, 'cta', 1);
            fft2 = 1;
            rry = 0;
        for p=1:avpar,
            yy = y(p:end, :);
            yy = yy(:, p:end);
            ax1 = ax;
            ax1.x = ax.x(p:end, 1);
            ax1.y = ax.y(p:end, 1);            
            fft.opt = 'real';
            [tax, ty] = fftprocessing(ax1, yy, fft);
            fft.opt = 'imag';
            tax.y = tax.x;
            tax.x = ax1.y;
            [rax, ry] = fftprocessing(tax, ty.', fft);
            
            rry = rry + abs(ry);
        end
        %%%%%%% alsi %%%%%%%%%%%%% turn back 29.05.07
        tax.y = rax.x;
        tax.x = rax.y;
        rax = tax;
        ry=rry.';
        %%%%%%% alsi %%%%%%%%%%%%% 29.05.07
    elseif strcmp(opt, 'reim4'), % 2D FFT: fft(fft(real(y)))
        fft2 = 1;
        fft.opt = 'imag';
        [tax, ty] = fftprocessing(ax, real(y), fft);
        fft.opt = 'imag';
        tax.y = tax.x;
        tax.x = ax.y;
        [rax, ry] = fftprocessing(tax, ty.', fft);
        %%%%%%% alsi %%%%%%%%%%%%% turn back 29.05.07
        tax.y = rax.x;
        tax.x = rax.y;
        ty = ry.';
        rax = tax;
        ry=ty;
        %%%%%%% alsi %%%%%%%%%%%%% 29.05.07   
        %%%%%%% alsi %%%%%%%%%%%%% 05.09.2011
        elseif strcmp(opt, 'reim4abs'), % 2D FFT: fft(fft(real(y)))
        fft2 = 1;
        fft.opt = 'imag';
        [tax, ty] = fftprocessing(ax, real(y), fft);
        fft.opt = 'imag';
        tax.y = tax.x;
        tax.x = ax.y;
        [rax, ry] = fftprocessing(tax, ty.', fft);
        %%%%%%% alsi %%%%%%%%%%%%% turn back 29.05.07
        tax.y = rax.x;
        tax.x = rax.y;
        ty = ry.';
        rax = tax;
        ry=abs(ty);
        %%%%%%% alsi %%%%%%%%%%%%% 29.05.07   
        %%%%%%% alsi %%%%%%%%%%%%% 05.09.2011    
    elseif strcmp(opt, 'reim43D')
        %% this is a trickery for 3D HYSCORE
        fft2 = 1;
        [sx, sy] = size(y);
        n = floor(sy/sx);
        if n-sy/sx ~=0
            error('for reim43D sx ~= N x sy');
        end
        data = [];
        axy = [];
        for ii = 1:n
            cutY = y(:, [1:sx] + (ii-1)*sx);
            tax = ax;
            tax.y = tax.x;
            fft.opt = 'imag';
            [tax, ty] = fftprocessing(tax, real(cutY), fft);
            fft.opt = 'imag';
            tax.y = tax.x;
            tax.x = ax.x;
            [rax, ry] = fftprocessing(tax, ty.', fft);

            tax.y = rax.x;
            tax.x = rax.y;
            ty = ry.';
            
            maxy = max(tax.x);
            data = [data, ty];
            axy = [axy; rax.y+2*maxy*(ii-1)];
            
        end
        
        rax = tax;
        rax.y = axy;
        ry = data;
        %%%%%%% alsi %%%%%%%%%%%%% 05.09.2011         
    elseif strcmp(lower(opt), 'ftir'), % 2D FFT: fft(fft(real(y)))
        fft2 = 1;
        fft.opt = 'real';
        midpoint = safeget(fft, 'midpoint', floor(length(y)/2));
        ty1 = real(y(midpoint:-1:1));
        ty2 = real(y((midpoint):end));
        if length(ty1) > length(ty2)
            ty1 = ty1(1:length(ty2));
        else
            ty2 = ty2(1:length(ty1));
        end
        ax.x = ax.x(1:length(ty1));
        [ttax, tty] = fftprocessing(ax, ty1, fft);
        [ttax, tty2] = fftprocessing(ax, ty2, fft);
        
        rax = ttax;
        ry=abs(tty + tty2)/2;

    else
        fft2 = 0;
        if hand.Projection == 1
            [rax, ry] = fftprocessing(ax, y, fft);
        else
            ax1 = ax;
            ax1.x      = ax.y; ax1.xlabel = ax.ylabel;
            ax1.y      = ax.x; ax1.ylabel = ax.xlabel;
            [ax1, ry] = fftprocessing(ax1, y', fft);
            ry = ry'; 
            rax = ax;
            rax.x = ax1.y; rax.xlabel = ax1.ylabel;
            rax.y = ax1.x; rax.ylabel = ax1.xlabel;
        end
    end
    hand.out{1} = hand.src;
    hand.out{1}.ax = rax;
    hand.out{1}.ax.dx = 0;
    hand.out{1}.y = ry;
    hand.out{1}.ax.Color = [0.3, 0.3, 1]; % always blue
    hand.out{1}.ax.s = 0;
    hand.out{1}.title = 'FFT';
else
    opt = safeget(fft, 'opt', 'real');
    src.ax = hand.src.ax;
    if strcmp(opt, 'ftir')
        midpoint = safeget(fft, 'midpoint', floor(length(y)/2));
        ty1 = real(y(midpoint:-1:1));
        ty2 = real(y((midpoint):end));
        if length(ty1) > length(ty2)
            ty1 = ty1(1:length(ty2));
        else
            ty2 = ty2(1:length(ty1));
        end
        y = [ty1, ty2];
        src.ax.x = ax.x(1:length(ty1));
    end
    
    
    hand.out{end+1}.ax = src.ax;
    hand.out{end}.y = y; 
    hand.out{end}.ax.Color = [0, 0, 1]; % always blue
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Src-BL';
    
    handles = guidata(handles.figure1);
    
    % load script 
    str = get(handles.list, 'String');
    fft = [];
    for k=1:size(str,1),
        try,eval([str{k} ';']);end
    end
    
    
    fft.fft = 0; % do no FFT
    if strcmp(opt, 'reim2'), % if 2D FFT
        fft.opt = 'real';
        [tax, ty, tawin] = fftprocessing(src.ax, y, fft);
        fft.opt = 'imag';
        tax.y = tax.x;
        tax.x = ax.y;
        [rax, rty, ttawin] = fftprocessing(tax, ty.', fft);
        ry = rty.';
        rawin = ttawin.*tawin.';
        % now data is rotated, not good
    else
        if hand.Projection == 1
            [rax, ry, rawin] = fftprocessing(src.ax, y, fft);
        else
            ax1 = ax;
            ax1.x      = ax.y; ax1.xlabel = ax.ylabel;
            ax1.y      = ax.x; ax1.ylabel = ax.xlabel;
            [ax1, ry, rawin] = fftprocessing(ax1, y', fft);
            ry = ry'; rawin = rawin';
            rax = ax;
            rax.x      = ax1.y; rax.xlabel = ax1.ylabel;
        end
    end
    nawin = max(rawin);    % try to cause error
    nyr = max(max(real(hand.out{1}.y)));
    nyi = max(max(imag(hand.out{1}.y)));
    hand.out{end+1}.ax.x = src.ax.x;
    hand.out{end}.ax.y = src.ax.y;
    hand.out{end}.y = rawin*(nyr+j*nyi);
    hand.out{end}.ax.Color = [0, 1, 0]; % always green
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Awin';
    % application of window
    hand.out{end+1}.ax = src.ax;
    hand.out{end}.y = ry;
    hand.out{end}.ax.Color = [0, 0, 0]; % always black
    hand.out{end}.ax.s = 0;
    hand.out{end}.title = 'Out';
    % show baseline 
    if ~basln,
        hand.out{end+1}.ax = hand.src.ax;
        hand.out{end}.ax.y = hand.src.ax.y;    
        hand.out{end}.y = bline;
        hand.out{end}.ax.Color = [1, 0, 0]; % red
        hand.out{end}.ax.s = 0;
        hand.out{end}.title = 'Bline';
    end
end


guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function LoSc_Callback(hObject, eventdata, handles)
full = which('kazan');
[fpath,name,ext] = fileparts(full); 
[fname,pname] = uigetfile([fpath, filesep, 'pepper_scr',filesep,'*.m'],'Open Script');
if fname == 0, return; end

if ischar(get(handles.popSel, 'String')), 
  str = {get(handles.popSel, 'String')};
else
  str = get(handles.popSel, 'String');
end
val = get(handles.popSel, 'Value');

[name, scr] = LoadScript([pname filesep fname], handles.sim{val}.stype);
set(handles.list, 'String', scr, 'Value', 1);
list_Callback(hObject, eventdata, handles);

str{val} = [num2str(val), ': ', fname];
set(handles.popSel, 'String', str);
handles.sim{val}.Script = scr;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = LoadScript(fullname, stype)
fid = fopen(fullname, 'r');
[fpath,fname,ext] = fileparts(fullname); 
if fid == -1, disp('File not found'); return; end
str = {};
while feof(fid) < 1,
  st = fgetl(fid);
  if ~isempty(st) & st(1) ~= '%',
    str{end+1} = st;
  end
end
fclose(fid);

if nargout > 0 , 
  varargout = {fname , str.'};
end

% --------------------------------------------------------------------
function SvSc_Callback(hObject, eventdata, handles)

[fname,pname] = uiputfile([filesep, '*.m'],'Save Script');
if fname == 0, return; end
[fpath,fname,ext] = fileparts([pname, fname]); 
if strcmp(ext ,'') , ext = '.m'; end
fid = fopen([pname, fname, ext], 'w');
str = get(handles.list, 'String');
for i = 1:size(str, 1)
  fprintf(fid, '%s\r\n', str{i});
end
fclose(fid);

% ------------------------------------------------------------
function check_Callback(hObject, eventdata, handles)
c = {'off', 'on'};
col = {[1 0.7 0.7], [0.7, 1, 0.7]};
stat = strcmp(get(handles.mCalc, 'Checked'), 'on');
set(handles.mCalc, 'Checked', [c{~stat + 1}]);
set(handles.butCalc, 'Value', ~stat,'BackgroundColor', [col{~stat + 1}]);

if ~stat,
    ePars_Callback(hObject, eventdata, handles)
else
    set(handles.butFFT, 'Value',0);
end

% ------------------------------------------------------------
function figure1_CloseRequestFcn(hObject, eventdata, handles)
try
  hand = guidata(handles.hh);
  hand.rx = {};
  hand.ry = {};
  hand.rc = {};
  hand.rs = {};
  guidata(handles.hh, hand);
  eval([name '(''SetDataSource'', 2, hand)']);
end
delete(handles.figure1);

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

% --- Executes on button press in butDel.
function butDel_Callback(hObject, eventdata, handles)
num = handles.num;
val = get(handles.popSel, 'Value');
if num == 1, return; end

switch val
  case 1,
    ss = 2:num;
  case num
    ss = 1:num-1;
  otherwise
    ss = [1:val-1, val+1:num];
end

handles.sim = {handles.sim{ss}};
handles.num = handles.num - 1;
str = get(handles.popSel, 'String');
set(handles.popSel, 'String', str{ss}, 'Value', 1);
popSel_Callback(hObject, eventdata, handles);
guidata(handles.figure1, handles);

% ------------------------------------------------------------
function popSel_Callback(hObject, eventdata, handles)
num = get(handles.popSel, 'Value');
set(handles.list, 'String', handles.sim{num}.Script, 'Value', min(2, size(handles.sim{num}.Script, 1)));
list_Callback(hObject, eventdata, handles);
set(handles.butCol, 'BackgroundColor', handles.Color(num, :)); 

% --------------------------------------------------------------------
function AddSim_Callback(hObject, eventdata, handles)
full = get(handles.figure1, 'FileName');
[fpath,name,ext] = fileparts(full); 
[fname, pname] = uigetfile([fpath ,filesep,'pepper_scr',filesep,'*.m'],'Open Script');
if fname == 0, return; end

switch hObject
case handles.mFileBaseline
  stype = 'baseline';
otherwise
  stype = 'sim';
end

[name, script] = LoadScript([pname fname], stype);
set(handles.list, 'String', script, 'Value', 1);
list_Callback(hObject, eventdata, handles);

if ischar(get(handles.popSel, 'String')),
  str = {get(handles.popSel, 'String')};
else
  str = get(handles.popSel, 'String');
end

str{handles.num + 1} = [num2str(handles.num + 1), ': ', name];
set(handles.popSel, 'String', str, 'Value', handles.num + 1);

hand = guidata(handles.hh);
handles.sim{handles.num + 1}.Script = script;
handles.sim{handles.num + 1}.x = hand.src.ax.x;
handles.sim{handles.num + 1}.y = zeros(size(hand.src.ax.x));
handles.sim{handles.num + 1}.stype = stype;
handles.num = handles.num + 1;
guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function mSimShow_Callback(hObject, eventdata, handles)
c = {'off', 'on'};
stat = strcmp(get(hObject, 'Checked'), 'on');
set(hObject, 'Checked', [c{~stat + 1}]);
if ~stat, 
  plothis(handles);
end

% --------------------------------------------------------------------
function varargout = mFileSalt_Pepper_Callback(h, eventdata, handles, varargin)
func = {'Pepper', 'Salt'}; % For future, when Num. of Functions will be > 2
Name = get(h, 'Label');
for i = 1:size(func, 2)
  if strcmp(func{i}, Name), 
    set(h, 'Checked', 'on');
    handles.CalcType = lower(Name);
  else
    eval(['set(handles.',['mFile', func{i}], ', ''Checked'', ''off'');']);
  end
end
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = mFileCopy_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function mColSum_Callback(hObject, eventdata, handles)
num = get(handles.popSel, 'Value');
c = uisetcolor(handles.Color(num, :));
if size(c, 2) > 1
  handles.SumCol = c;
  %set(handles.butCol, 'BackgroundColor', c, 'ForegroundColor', [1 1 1] - c);
  guidata(handles.figure1, handles);
  plothis(handles);
end

% --------------------------------------------------------------------
function mFixAxis_Callback(h, eventdata, handles)
c = {'off', 'on'};
stat = strcmp(get(handles.mFixAxis, 'Checked'), 'on');
set(handles.mFixAxis, 'Checked', [c{~stat + 1}]);
if ~stat, 
  plothis(handles);
end

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

% --- Executes on button press in butCalc.
function butCalc_Callback(hObject, eventdata, handles)
check_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function butFFT_Callback(hObject, eventdata, handles)
c = {'off', 'on'};
stat = strcmp(get(handles.mCalc, 'Checked'), 'on');

if ~stat & get(handles.butFFT, 'Value')
  col = {[1 0.7 0.7], [0.7, 1, 0.7]};
  set(handles.mCalc, 'Checked', [c{~stat + 1}]);
  set(handles.butCalc, 'Value', ~stat,'BackgroundColor', [col{~stat + 1}]);
end

ePars_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function mFileScript_Callback(hObject, eventdata, handles)
if hObject==handles.mFileExecScript
    for k=1:size(handles.script,2),
        disp(handles.script{k})
        evalin('base',handles.script{k})
    end
else
    disp(' ')
    for k=1:size(handles.script,2),
        disp(handles.script{k})
    end
end

% --------------------------------------------------------------------
function mAbout_Callback(hObject, eventdata, handles)
msgbox({'FFTPLUGIN'; 'by Boris Epel and Alexey Silakov, 2003-07';''; ...
        'available parameters are:'; ...
        'awin - apodization window';...
        '    = bla/bar/con/cos/ham/han/wel/exp/gau/kai'; ...
        ''; ...
        'awidth, aalpha, ashift - apodization window pars';...
        ''; ...
        'zerofill - fill up to 2^(log2(size)+zerofill) points <0>';...
        ''; ...
        'rshift, lshift - shift of data (points) <0>';...
        ''; ...
        'phase0 - zero-order phase correction (deg) <0>';...
        ''; ...
        'xscale - x axis coefficient <1>';...
        ''; ...
        'opt = real/imag/qseq/cta/reim2/reim4 - type of FT <real>';...
        '   qseq - sequential acquisition';...
        '   cta  - cross-term averaging';...
        '   reim2 - 2D fft. FFT(FFT(real(y)))';...   
        '   reim4 - also 2D fft. FFT(FFT(real(y))) but gives all 4 quadrants';...           
        ''; ...
        'cta - cross-term averaging parameter <datasize/4>';...
        ''; ...
        'bsline = none/polyfit2D - automatic baseline';...
        '   n_pol (for ''polyfit2D'') - degree of the polynom to fit';...
        ''; ...
        'In plugin'},...
      'About', 'help')

% --------------------------------------------------------------------
function butFFTall_Callback(hObject, eventdata, handles)
ePars_Callback(hObject, eventdata, handles)

%%%%%%% AlSi 10.09.2004 %%%%%%%%%%%%%%%
function out = polyfit1(x, y, n, dir)
switch dir
   case 1
      str = '(:, ii)';
      tr = '';
      dd = 2;
   otherwise
      str = '(ii, :)';
      tr = '.''';
      dd = 1;
end
for ii = 1:size(y, dd)
   eval(['ty(:, 1) = real(y', str, tr, ');']);
   pp = polyfit(x, ty, n);
   eval(['out', str,' = polyval(pp, x',tr,');']);
end

%-----------------------------------------------------
function out = polyfit2(x, x1, y, n)
tout = y*0;
if size(y, 2) > 1,
   for ii = 1:size(y, 2)
      pp = polyfit(x, real(y(:, ii)), n);
      tout(:, ii) = polyval(pp, x);
   end   
end
out = tout;
if size(y, 1) >1
for jj = 1:size(y, 1)
   pp = polyfit(x1, (real(y(jj, :)) -tout(jj, :)).', n);
   out(jj, :) = polyval(pp, x1).' + tout(jj, :);
end
end

%-----------------------------------------------------
function AddIntBl_Callback(hObject, eventdata, handles)
switch hObject
    case handles.mBl1DPolyfit
        Scr = {'f = polyfit1(x, y, n, 1)'; 'n = 2'};
        titlestr = '1D polyfit';
    case handles.mBl2DPolyfit
        Scr = {'f = polyfit2(x, x1, y, n)'; 'n = 2'};        
        titlestr = '2D polyfit';
    otherwise
end
if ischar(get(handles.popSel, 'String')), 
  str = {get(handles.popSel, 'String')};
else
  str = get(handles.popSel, 'String');
end
val = length(str) + 1;
%val = get(handles.popSel, 'Value');

set(handles.list, 'String', Scr, 'Value', 1);
list_Callback(hObject, eventdata, handles);

str{val} = [num2str(val), ': ', titlestr];
set(handles.popSel, 'String', str, 'Value', val);


hand = guidata(handles.hh);
handles.sim{val}.Script = Scr;
handles.sim{val}.x = hand.src.ax.x;
handles.sim{val}.x1 = hand.src.ax.y;
handles.sim{val}.y = (hand.src.y)*0;
handles.sim{val}.stype = 'baseline';
handles.num = handles.num + 1;
guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);

guidata(handles.figure1, handles);
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
       
