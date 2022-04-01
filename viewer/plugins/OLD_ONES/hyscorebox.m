function varargout = hyscorebox(varargin)
% HYSCOREBOX Application M-file for hyscorebox.fig
%    Intended to use only with 'kazan' viewer

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 05-Oct-2004 22:54:47
% Alexey Silakov & Boris Epel 5-dec-2003, MPI
% boep 10-dec-2003, MPI

if nargin == 0  % LAUNCH GUI
    token = 'HYSCOREBOX_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,n,e] = fileparts(which('kazan'));
    inifilename = [fpath '\kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2; % how to behave if HYSCOREBOX window is already exist
    try opc = ini.KazanViewer.MultyKazan; end % try to get parameter from ini file
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end

	fig = openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
        load('ico.mat');
    handles.tbMain = uitoolbar(handles.figure1);    
    handles.ptChange = uipushtool(handles.tbMain, 'TooltipString', 'Change script by loading file', 'Tag', 'ptChange', ...
        'CData', ico.change, 'ClickedCallback', 'hyscorebox(''mLoad_Callback'',gcbo,[],guidata(gcbo))');    
    handles.ptOpen = uipushtool(handles.tbMain, 'TooltipString', 'Add script from file', 'Tag', 'ptOpen', ...
        'CData', ico.load, 'ClickedCallback', 'hyscorebox(''AddSim_Callback'',gcbo,[],guidata(gcbo))');
    handles.ptSave = uipushtool(handles.tbMain, 'TooltipString', 'Save current script', 'Tag', 'ptSave',...
        'CData', ico.save, 'ClickedCallback', 'hyscorebox(''mSave_Callback'',gcbo,[],guidata(gcbo))');
    handles.ptCopy = uipushtool(handles.tbMain, 'TooltipString', 'Duplicate Script', 'Tag', 'ptCopy', ...
        'CData', ico.copy, 'ClickedCallback', 'hyscorebox(''mFileCopy_Callback'',gcbo,[],guidata(gcbo))');
    handles.ptEdit = uipushtool(handles.tbMain, 'TooltipString', 'Change Script in editor', 'Tag', 'ptEdit', ...
        'CData', ico.edit, 'ClickedCallback', 'hyscorebox(''mFileOpenEditor_Callback'',gcbo,[],guidata(gcbo))');
    handles.ptClipboard = uipushtool(handles.tbMain, 'TooltipString', 'Copy to Clipboard', 'Tag', 'ptClipboard',...
        'CData', ico.clipboard,'ClickedCallback', 'hyscorebox(''mSave_Callback'',gcbo,[],guidata(gcbo))');    
    
    handles.hh = 0;
    handles.ploter = 0;
    handles.num = 1;
    handles.div = 0;
    handles.Spec = {};
    handles.sim{1}.Script = get(handles.list, 'String');
    handles.sim{1}.x1 = [];
    handles.sim{1}.x2 = [];    
    handles.sim{1}.y = [];
    handles.sim{1}.stype = 'sim';
    handles.CalcType = 2;
    handles.Quad = [1, 1, 1, 1];
    handles.Color =   [1, 1, 1;...
            1, 1, 1;...
            1, 1, 1;...
            1, 1, 1;...
            1, 1, 1;...
            1, 1, 1;...
            1, 1, 1];
    handles.sepstor = 1; % flag to store different HYSCORE sectra in different cels (default 0)
    set(handles.mLoad, 'Accelerator', 'L');
    set(handles.mSave, 'Accelerator', 'S');
    set(handles.mCalc, 'Accelerator', 'C');
    set(handles.AddSim, 'Accelerator', 'A');
    set(handles.mFileMoveUp, 'Accelerator', 'P');
    set(handles.mFileMoveDown, 'Accelerator', 'O');
    set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
    set(handles.ePars, 'UserData', ' ');
    set(handles.mSimAgr, 'Checked', 'on');
    handles.process = 1;
    Figure_ResizeFcn(handles.figure1, [], handles)
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        %     [varargout{1:nargout}] = 
        feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
        fig = openfig(mfilename,'reuse');
        handles = guihandles(fig);
        disperr(handles);
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

function FileUpdate(handles)

% function varargout = mFileSave(h, eventdata, handles, varargin)
% disp('I don''t know what this function must do');
% % handles = guidata(handles.figure1);
% % [fname,pname] = uiputfile(['\*.*'],'Save Script');
% % if fname == 0, return; end
% % 
% % y = real(handles.y);
% % x1 = handles.ax.x + handles.ax.dx;
% % x2 = handles.ax.y + handles.ax.dy;
% % p = [x,y];
% % save([pname, fname], 'p', '-ascii')
% --------------------------------------------------------------------
function varargout = Figure_ResizeFcn(h, eventdata, handles, varargin)
% 140 210
FigaPos = get(handles.figure1, 'Position');
if FigaPos(3)<130, FigaPos(3) = 130; set(handles.figure1, 'Position', FigaPos); end
if FigaPos(4)<210, FigaPos(4) = 210; set(handles.figure1, 'Position', FigaPos); end
brd = 5; % board
butCalcH = 24; % points
UpperFrameH = 25;
EditH = 15;
TopFigure = 20+brd;
CheckH = 15;
CalcPos = [brd, brd, FigaPos(3)-2*brd, butCalcH];
FramePos = [brd, brd*2 + butCalcH, FigaPos(3)-2*brd,...
             FigaPos(4) - UpperFrameH-TopFigure - 4*brd];
EditPos = [FramePos(1)+brd, FramePos(2)+brd,FramePos(3)-2*brd, EditH];
SlPos = [EditPos(1), EditPos(2) + EditPos(4)+brd-1, 25, 16];
PopPos = [SlPos(1) + SlPos(3)+brd, EditPos(2)+EditPos(4)+brd-3,...
        EditPos(3)-brd - SlPos(3),18];

ListPos = [SlPos(1), SlPos(2)+SlPos(4)+brd, EditPos(3),...
        FramePos(4)-EditPos(4)-SlPos(4)-CheckH- 4*brd];
CheckV = ceil((ListPos(3) - 3*brd)/4);
Check1Pos = [ListPos(1), ListPos(2)+ListPos(4), CheckV, CheckH];
Check2Pos = [ListPos(1)+CheckV+brd, ListPos(2)+ListPos(4), CheckV, CheckH];
Check3Pos = [ListPos(1)+CheckV*2+2*brd, ListPos(2)+ListPos(4), CheckV, CheckH];
Check4Pos = [ListPos(1)+CheckV*3+3*brd, ListPos(2)+ListPos(4), CheckV, CheckH];

TopFramePos = [FramePos(1), FramePos(2)+FramePos(4)+brd, FramePos(3), UpperFrameH];
PopSelPos = [TopFramePos(1)+brd, TopFramePos(2)+brd, ...
        3*(TopFramePos(3)-2*brd)/4-brd, TopFramePos(4)-2*brd];
ColPos = [TopFramePos(1)+ PopSelPos(3) + 2*brd, TopFramePos(2)+brd, ...
        (TopFramePos(3)-2*brd)/4, TopFramePos(4)-2*brd];

set(handles.butCalc,'Units', 'points',  'Position', CalcPos);
set(handles.frame1, 'Units', 'points', 'Position', FramePos);
set(handles.ePars, 'Units', 'points', 'Position', EditPos);
set(handles.slPars, 'Units', 'points', 'Position', SlPos);
set(handles.popup, 'Units', 'points', 'Position', PopPos);
set(handles.list, 'Units', 'points', 'Position', ListPos);
set(handles.chPlPl, 'Units', 'points', 'Position', Check1Pos);
set(handles.chPlMin, 'Units', 'points', 'Position', Check2Pos);
set(handles.chMinPl, 'Units', 'points', 'Position', Check3Pos);
set(handles.chMinMin, 'Units', 'points', 'Position', Check4Pos);
set(handles.frame3, 'Units', 'points', 'Position', TopFramePos);
set(handles.popSel, 'Units', 'points', 'Position', PopSelPos);
set(handles.butCol, 'Units', 'points', 'Position', ColPos);

% --------------------------------------------------------------------
function varargout = list_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');

[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(2:end)), ' '];

if strcmp(trim(name), 'f')
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'function');
    set(handles.slPars, 'Enable', 'inactive');  
    set(handles.popup, 'Enable', 'off');    
elseif (str2(1)=='''')
    primes = findstr(str2, '''');
    if size(primes, 2) < 2, primes(2) = size(str2, 2); end
    set(handles.ePars, 'String', str2(2:(primes(2)-1)));
    set(handles.ePars, 'UserData', 'string');
    set(handles.slPars, 'Enable', 'off');    
    set(handles.popup, 'Enable', 'off');    
elseif (~isempty(findstr(str2, 'lin')))|(str2(1)=='[')
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'range');    
    set(handles.slPars, 'Enable', 'off');
    set(handles.popup, 'Enable', 'off');
else
    dig = str2num(str2);
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'value');
    set(handles.slPars, 'Enable', 'on');
    set(handles.popup, 'Enable', 'on');    
end

% --------------------------------------------------------------------
function varargout = ePars_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ePars.
handles = guidata(handles.figure1);
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
    set(handles.ePars, 'String', newval);
    Res{val} = [name ' = ''', newval,''''];
case 'range'    
    set(handles.ePars, 'String', newval);
    Res{val} = [name ' = ', newval];
otherwise
    nums = str2num(newval);
    set(handles.ePars, 'String', trim(newval));
    Res{val} = [name ' = ' trim(newval)];
end
set(handles.list, 'String', Res);
list_Callback(h, eventdata, handles);
handles.sim{get(handles.popSel, 'Value')}.Script = Res;
guidata(handles.figure1, handles);

if isempty(strfind(name,'Shift.'))
    if strcmp(get(handles.mCalc, 'Checked'), 'on'), 
        nn = get(handles.popSel, 'Value');
        handles = hycalc(get(handles.popSel, 'Value'), handles);
        plothis(handles);
        if strcmp(get(handles.mSimAutoZoom, 'Checked'), 'on')
            ZZ_Callback(0, [], handles);
        end
    end
else
    plothis(handles);
end

% --------------------------------------------------------------------
function varargout = eResize_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
handles.div = str2num(get(handles.eResize, 'String'));
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function varargout = slPars_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;
shift = get(h, 'Value');
set(h, 'Value', 0.5);
CurVal = str2num(get(handles.ePars, 'String'));
num = get(handles.popup, 'Value');
str = get(handles.popup, 'String');
dim = str2num([str{num}]);
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.ePars, 'String', num2str(Val, 15));
ePars_Callback(h, eventdata, handles, varargin);

% --------------------------------------------------------------------
function varargout = slider2_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
shift = get(handles.slider2, 'Value');
set(handles.slider2, 'Value', 0.5);
num = get(handles.popup2, 'Value');
str = get(handles.popup2, 'String');
dim = str2num([str{num}]);

CurVal = str2num(get(handles.eResize, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.eResize, 'String', num2str(Val));
eResize_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function mSimShow_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
c = {'off', 'on'};
stat = strcmp(get(hObject, 'Checked'), 'on');
set(hObject, 'Checked', [c{~stat + 1}]);
if ~stat, 
    plothis(handles);
end 

% --------------------------------------------------------------------
function varargout = hycalc(num, handles)
handles = guidata(handles.figure1);
% load script 
Sys = []; Exp = []; Opt = []; Shift = [];
str = get(handles.list, 'String');
for k=1:size(str,1),
    try,eval([str{k} ';']);end
end
if Sys.S>=1
      disp('hycalc: only S = 1/2 available');
      return;
end
%%%%%%%%%% IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% in script:                                   %
% All angles must be in degree, NOT in radian  %
% Field in mT                                  %
% Frequency, will be in MHz                    %
% Exept mwFreq, this parameter has to be in GHz%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end

switch safeget(Opt, 'Sim', 'sketch')
    case 'sketch'
        handles.OverlayOptions{num} = [1, 1];
        handles.CalcType(num) = 1;
        handles.SimSrc = 1;
        guidata(handles.figure1, handles);
        sketch_only(Sys, Exp, Opt, Shift, handles);
    case 'sketchhfdec'
        handles.OverlayOptions{num} = [1, 1];
        handles.CalcType(num) = 1;
        handles.SimSrc = 1;
        guidata(handles.figure1, handles);
        sketchhfdec(Sys, Exp, Opt, Shift, handles);        
    case 'sketch2DeseB'
        handles.OverlayOptions{num} = [1, 1];        
        handles.CalcType(num) = 1;
        handles.SimSrc = 1;
        guidata(handles.figure1, handles);
        sketch2DeseB(Sys, Exp, Opt, Shift, handles);  
    case 'fd'
        set(handles.figure1, 'Name', 'HYSCORE: busy');
        cc = 1;
        if handles.sepstor, cc = num; end
        handles.sim{cc}.x1 = [];
        handles.sim{cc}.x2 = [];
        handles.sim{cc}.y  = [];
        switch Sys.I(1)
            case {1/2, 1, 3/2, 2, 5/2}
                [y, ax] = kv_hyscorefd(Sys, Exp, Opt, Shift);
        end
        handles.sim{cc}.x1 = ax.x*1e-6; %[Hz] -> [MHz]
        handles.sim{cc}.x2 = ax.y*1e-6; %[Hz] -> [MHz]
        handles.sim{cc}.y = y;
        if isfield(Shift, 'contour')
            handles.sim{cc}.contour = Shift.contour;
        end
        handles.OverlayOptions{num} = [2, 4];        
        if ~safeget(Shift, 'Overlay', 0)
            handles.CalcType(num) = 2;
            handles.SimSrc = 0; 
        else
            handles.CalcType(num) = 4;
            handles.SimSrc = 1;
        end
        guidata(handles.figure1, handles);
        set(handles.figure1, 'Name', 'HYSCORE');
        dispnorm(handles);
    case 'td'
        set(handles.figure1, 'Name', 'HYSCORE: busy');
        cc = 1;
        if handles.sepstor, cc = num; end
        handles.sim{cc}.x1 = [];
        handles.sim{cc}.x2 = [];
        handles.sim{cc}.y  = [];
        switch Sys.I(1)
            case {1/2, 1, 5/2}
                [y, ax] = kv_hyscoretd(Sys, Exp, Opt, Shift);
        end
        handles.sim{cc}.x1 = ax.x*1e-6; %[Hz] -> [MHz]
        handles.sim{cc}.x2 = ax.y*1e-6; %[Hz] -> [MHz]
        handles.sim{cc}.y = y;
        if isfield(Shift, 'contour')
            handles.sim{cc}.contour = Shift.contour;
        end
        handles.OverlayOptions{num} = [2, 4];        
        if ~safeget(Shift, 'Overlay', 0)
            handles.CalcType(num) = 2;
            handles.SimSrc = 0; 
        else
            handles.CalcType(num) = 4;
            handles.SimSrc = 1;
        end
        guidata(handles.figure1, handles);
        set(handles.figure1, 'Name', 'HYSCORE');
        dispnorm(handles);        
    case 'kvmr22d'
        set(handles.figure1, 'Name', 'HYSCORE: busy');
        cc = 1;
        if handles.sepstor, cc = num; end
        handles.sim{cc}.x1 = [];
        handles.sim{cc}.x2 = [];
        handles.sim{cc}.y  = [];
        axField = linspace(Exp.Field(1), Exp.Field(2), Exp.nFPoints);
        switch Sys.I(1)
            case {1/2, 1}
                ty = zeros(Exp.nPoints, Exp.nFPoints);
                for ii = 1:Exp.nFPoints
                    Exp.Field = axField(ii);
                    [x, y] = kvmr2(Sys, Exp, Opt);
                    ty(:, ii) = y(:);
                end
        end
        handles.sim{cc}.x1 = x;
        handles.sim{cc}.x2 = axField';
        handles.sim{cc}.y = ty;
        if isfield(Shift, 'contour')
            handles.sim{cc}.contour = Shift.contour;
        end
        handles.OverlayOptions{num} = [3, 4];
        if ~safeget(Shift, 'Overlay', 0)
            handles.CalcType(num) = 3;
            handles.SimSrc = 0; 
        else
            handles.CalcType(num) = 4;
            handles.SimSrc = 1;
        end
        guidata(handles.figure1, handles);
        set(handles.figure1, 'Name', 'HYSCORE');
        dispnorm(handles);    
    case 'saffron'
        % new version of EasySpin has the option to calculate HYSCORE
        %
        % BUT, tau is in [us] instead of [ns]
        % BUT, Apa/Qpa is in [rad], not in [deg]
        % BUT, there is no way it specifies gn/I, you have to use Sys.Nucs ( :/ )
        % BUT, no MaxFreq, one should use dt instead ( :/ )
        set(handles.figure1, 'Name', 'HYSCORE: busy');
        cc = 1;
        if handles.sepstor, cc = num; end
        handles.sim{cc}.x1 = [];
        handles.sim{cc}.x2 = [];
        handles.sim{cc}.y  = [];
        if ~safeget(Opt, 'Symmetry', 0)
            Opt.Symmetry = 'Ci';
        end
        if isfield(Exp, 'MaxFreq')
            Exp.dt = 1/2/Exp.MaxFreq;
        end

        Exp.tau = Exp.tau*1e-3; % ns to us
        if isfield(Sys, 'Apa')
            Sys.Apa = Sys.Apa*pi/180; % degree to rad
        end
        if isfield(Sys, 'Qpa')
            Sys.Qpa = Sys.Qpa*pi/180; % degree to rad        
        end
        if isfield(Exp, 'gField')
            Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
        end
        Exp.Sequence = 'HYSCORE';
        %%%%%%%%%% Stephan, you are an idiot !!!!
%         out=saffron(Sys, Exp, Opt);
        [x1, x2, S, out] = saffron(Sys, Exp, Opt);
            
        handles.sim{cc}.y = out.fd;
        handles.sim{cc}.x1 = out.f1(:);
        handles.sim{cc}.x2 = out.f2(:);
        
        if isfield(Shift, 'contour')
            handles.sim{cc}.contour = Shift.contour;
        end
        handles.OverlayOptions{num} = [3, 4];        
        if ~safeget(Shift, 'Overlay', 0)
            handles.CalcType(num) = 3;
            handles.SimSrc = 0; 
        else
            handles.CalcType(num) = 4;
            handles.SimSrc = 1;
        end
         guidata(handles.figure1, handles);
         
        set(handles.figure1, 'Name', 'HYSCORE');
        dispnorm(handles);   
    case 'mimstau'
        set(handles.figure1, 'Name', 'HYSCORE: busy');
        cc = 1;
        if handles.sepstor, cc = num; end
        handles.sim{cc}.x1 = [];
        handles.sim{cc}.x2 = [];
        handles.sim{cc}.y  = [];

        [y, ax] = mimstau(Sys, Exp, Opt);
        
        handles.sim{cc}.x1 = ax.x; %[MHz]
        handles.sim{cc}.x2 = ax.y; %[MHz]
        handles.sim{cc}.y = y;
        if isfield(Shift, 'contour')
            handles.sim{cc}.contour = Shift.contour;
        end
        handles.OverlayOptions{num} = [2, 4];        
        if ~safeget(Shift, 'Overlay', 0)
            handles.CalcType(num) = 2;
            handles.SimSrc = 0; 
        else
            handles.CalcType(num) = 4;
            handles.SimSrc = 1;
        end
        guidata(handles.figure1, handles);
        set(handles.figure1, 'Name', 'HYSCORE');
        dispnorm(handles);
    case 'gui'
        % alsi: sometimes i'm using this box as an interface for programs,
        % not related to KazanViewer
        handles.OverlayOptions{num} = [3, 3];        
        handles.CalcType(num) = 3;
        handles.SimSrc = 0;
        guidata(handles.figure1, handles);        
    otherwise
        disp('Wrong Opt.Sim parameter. Possible parameters are:');
        disp('       sketch, sketchhfdec, sketch2DeseB, fd, kvmr22d, saffron, gui');
end
if nargout > 0, 
    varargout{1} = handles;
end
% --------------------------------------------------------------------
function sketch_only(Sys, Exp, Opt, Shift ,handles)

tstart = cputime;
tau = safeget(Exp, 'tau', 0)*1e-9; %ns->s

%%% simulation for a few nuclei
Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I', 'nNucs'};
% lets find number of simmulaitons (parameters can be common, it means size(XX, 1)==1)
Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad
Sys.nNucs = safeget(Sys, 'nNucs', 1); % number of equivalent nuclei
sim_size = 1;
isfield_num = [];
for pk = 1:length(Parameters)
    if isfield(Sys, Parameters{pk})        
        tPar = getfield(Sys, Parameters{pk});
        isfield_num(end+1) = pk;        
        if size(tPar, 1)>1
            sim_size = size(tPar, 1);
        end
    end
end
if sim_size>1
    if size(Sys.I, 2) > 1,
        Sys.I = Sys.I.';
    end
    if size(Sys.gn, 2) > 1,
        Sys.gn = Sys.gn.';
    end
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk); % kind of optimisation
        tPar = getfield(Sys, Parameters{pk});
        if size(tPar, 1)==1
            Sys = setfield(Sys, Parameters{pk}, tPar(ones(sim_size, 1), :));
        elseif size(tPar, 1)~=sim_size
%             error({['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.']; ['Correct number is ', num2str(sim_size)]});
              error(['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.'], ['Correct number is ', num2str(sim_size)]);
        end
    end
end

if isfield(Sys, 'Aiso'),
    Sys.A = safeget(Sys, 'A', 0) + Sys.Aiso;
end
% 
%%% orientation selection

Opt.Treshold = safeget(Opt, 'Treshold', 1e-3);


%%%%%%%%%%%%%%%%%%%% new version 18.02.2006 %%%%%%%%%%%
Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq & ~isfield(Exp, 'Field'), ...
           error('ESEEM: there is no Exp.mwFreq to calculate B0');
   end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
else 
   Exp.Field = safeget(Exp, 'Field', 0);
end

[phi, theta, ak] = Local_orisel(Sys, Exp, Opt);

w1 = []; w2 = [];
num = get(handles.popSel, 'Value');
% Calculation
if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function

if ~isfield(Opt, 'Nuclei')
    Opt.Nuclei = 1:sim_size;
elseif max(Opt.Nuclei)>sim_size
    error('max(Opt.Nuclei) > number of nuclei');
end

Sy = Sys;
for ii = Opt.Nuclei
    Sy.A = Sys.A(ii, :);
    Sy.Apa = Sys.Apa(ii, :);

    Sy.I = Sys.I(ii);
    Sy.gn = Sys.gn(ii);
    
   switch Sy.I(1)
   case 0.5
      handles.sim{num}.x1 = [];
      handles.sim{num}.x2 = [];
      handles.sim{num}.y  = [];
      for cth =1:length(theta)
         if ak(cth)>Opt.Treshold,
            freq = kv_hscfreq(Sy, Exp, phi(cth), theta(cth));
            w1(end+1, :) = freq(1); % MHz
            w2(end+1, :) = freq(2); % MHz
         end
      end
         handles.sim{num}.x1 = [w1 w2];
         handles.sim{num}.x2 = [w2 w1];
         handles.sim{num}.y  = [1, 1];

   case {1,1.5}
      handles.sim{num}.x1 = [];
      handles.sim{num}.x2 = [];
      handles.sim{num}.y  = [];     
      if isfield(Sys, 'Q')
          Sy.Q = Sys.Q(ii, :);
          Sy.Qpa = Sys.Qpa(ii, :);
      else
          Sy.Q = [0, 0, 0];
          Sy.Qpa = [0, 0, 0];
      end
      for cth =1:length(theta)
            if ak(cth)>Opt.Treshold,       
               freq = kv_hscfreq(Sy, Exp, phi(cth), theta(cth));
               w1(end+1, :) = freq(:, 1)'; % MHz
               w2(end+1, :) = freq(:, 2)'; % MHz
            end   
      end
      
      if ~isfield(Opt, 'ShowCor'), Opt.ShowCor = [1, 1; 1, 2; 2, 1; 2, 2]; end
      
      if Opt.ShowCor == 0; 
          
          tRow = (1:(2*Sy.I+1));
          Opt.ShowCor = [reshape(tRow(ones((2*Sy.I+1), 1), :),(2*Sy.I+1)^2, 1), reshape(tRow(ones((2*Sy.I+1), 1), :).',(2*Sy.I+1)^2, 1)];
      end % Show all correlations
      
      for cc = 1:size(Opt.ShowCor, 1)   
         handles.sim{num}.x1 = [handles.sim{num}.x1 w1(:, Opt.ShowCor(cc, 1)) w2(:, Opt.ShowCor(cc, 2))];
         handles.sim{num}.x2 = [handles.sim{num}.x2 w2(:, Opt.ShowCor(cc, 2)) w1(:, Opt.ShowCor(cc, 1))];
         handles.sim{num}.y = [handles.sim{num}.y 1, 1];            
      end
      %          handles.sim{num}.x1 = [w1 w2];
      %          handles.sim{num}.x2 = [w2 w1];
      %          handles.sim{num}.y = [1, 1, 1, 1, 1, 1];
   otherwise
      disp('hycalc: only I = 1/2 and I = 1 available!!!!')
      return;
   end
end
handles.sim{num}.amp = abs(max(handles.sim{num}.y, [], 2)-min(handles.sim{num}.y, [], 2)).';
if isfield(Opt, 'LineStyle')
   handles.sim{num}.LineStyle = Opt.LineStyle;
end
if isfield(Opt, 'LineWidth')
   handles.sim{num}.LineWidth = Opt.LineWidth;
end
if isfield(Opt, 'Marker')
   handles.sim{num}.Marker = Opt.Marker;
end
% if max(handles.sim{num}.amp)==0
%     error('All spectra have zero amplitude');
% else
%     handles.sim{num}.y = handles.sim{num}.y./handles.sim{num}.amp(ones(size(handles.sim{num}.y))); 
% end

guidata(handles.figure1, handles);
dispnorm(handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate orientation selection and produce a set of angles
function [phi, theta, weights] = Local_orisel(Sys, Exp, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

if isfield(Exp, 'phi')& isfield(Exp, 'theta'),
    nphi = length(Exp.phi);
    ntheta = length(Exp.theta);

    theta = reshape(Exp.theta(ones(nphi, 1), :), nphi*ntheta, 1)*pi/180;
    phi = reshape(Exp.phi(ones(ntheta, 1), :).', nphi*ntheta, 1)*pi/180;
    ww = ones(nphi*ntheta, 1);
elseif isfield(Opt, 'Symmetry'),
    if ~Opt.Symmetry % if [0], use grid from MR2 (spiral [0:pi/2]-(0:2*pi))
        % taken from FORTRAN program of Ed Reijerse "MR2" 
        ncall = 0;
        nselect= 0;

        krid = safeget(Opt, 'nKnots', 20);

        theta = [];
        phi = [];
        last = 0;
        step = 0;
        while (~last)
            if(step==0 & ncall~=0), error(' PWDPEAL unitialized!'); end
            if(krid<0), error(' PWDPEAL illegal value of krid'); end
            if(ncall==0)
                itheta = krid;
                iphi   = 0;
                step = 90 / krid;
                nphi = 2*krid;
                thetaa = 90;
                dthe = -0.5 * step / nphi;
                dphi = 360 / nphi;
            end
            if(itheta<=0) error(' pdrpeal too many calls'); end
            if(iphi==nphi)
                iphi = 0;
                itheta = itheta -1;
                thetam = (itheta) * step;
                nphi = floor(sin(pi/180*thetam)*(4*krid));
                if(nphi<=0) error( 'BUG pdrpeal 1'); end
                % *              SMALLEST NPHI IS TYPICALLY 6
                thetaa = thetam + 0.5 * step;
                dthe = - step / (nphi);
                dphi = 360 / (nphi);
            end
            thetap = thetaa + (iphi) * dthe;
            phip   = (iphi) * dphi;
            iphi  = iphi + 1;
            ncall = ncall + 1;
            last = (iphi==nphi & itheta==1);
            if(last), step = 0; end

            theta(end+1, 1) = thetap*pi/180;
            phi(end+1, 1)   = phip*pi/180;
        end
        ww  = phi*0 + 1;
    else
        Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
        Opt.nKnots = safeget(Opt, 'nKnots', 20);
        [phi, theta, ww] = sphgrid(Opt.Symmetry, Opt.nKnots);
        phi = phi.'; theta = theta.'; ww = ww.';
    end

else
    t_theta = linspace(0, pi, safeget(Opt, 'nKnots', 20));
    Opt.nKnots = safeget(Opt, 'nKnots', 20);
    phi = [];
    theta = [];
    nnn = [];
    for cc = 1:Opt.nKnots
        num = floor(Opt.nKnots*abs(sin(t_theta(cc))));
        nnn(cc) = num;
        if num ~=0, 
            t_phi = linspace(0, 2*pi, num+1).';
            phi = [phi; t_phi(1:(end-1), 1)];
            theta = [theta; ones(num, 1)*t_theta(cc)];
        end
    end
    ww  = phi*0 + 1;

end

nangles = length(phi);
if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end % EasySpin function



% calculating orientation selection
ctheta = cos(theta);
stheta = sin(theta);
cphi = cos(phi);
sphi = sin(phi);
angle(:, 1)=stheta.*cphi;
angle(:, 2)=stheta.*sphi;
angle(:, 3) = ctheta;

Gx=angle(:, 1)*Sys.g(1);
Gy=angle(:, 2)*Sys.g(2);
Gz=angle(:, 3)*Sys.g(3);
geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

%%%%%%%%%%%%%%%%%%%% new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% include HF to the orientation selection 
useHFsel = safeget(Opt, 'useHFsel', 0);
if useHFsel
    projA = zeros(nangles, sum(Sys.nNucs));
    k = 0;
    for ii = 1:length(Sys.nNucs)
        for jj = 1:Sys.nNucs(ii)
            k = k+1;
            if sum(Sys.Apa(ii, :))
                R = erot(Sys.Apa(ii, :)*pi/180);
                At = R*diag(Sys.A(ii, :))*R';
            else
                At = diag(Sys.A(ii, :));
            end
            LA(:, 1) = angle(:, 1) * At(1, 1) + angle(:, 2) * At(2, 1) + angle(:, 3) * At(3, 1);
            LA(:, 2) = angle(:, 1) * At(1, 2) + angle(:, 2) * At(2, 2) + angle(:, 3) * At(3, 2);
            LA(:, 3) = angle(:, 1) * At(1, 3) + angle(:, 2) * At(2, 3) + angle(:, 3) * At(3, 3);
            projA(:, k) = sqrt(LA(:, 1).^2 + LA(:, 2).^2 + LA(:, 3).^2);
            mult(k) = 2*Sys.I(ii)+1;
            
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq & ~isfield(Exp, 'Field'), ...
           error('ESEEM: there is no Exp.mwFreq to calculate B0');
   end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
else 
    Exp.Field = safeget(Exp, 'Field', 0);
end

if (Exp.ExciteWidth > 0) & (Exp.mwFreq>0)
%     gcenter = fld2g(Exp.Field*1E-3, Exp.mwFreq*1E9); % killed by alsi 27.12.07

    feff = Exp.Field*1E-3*bmagn*geff/planck *1e-6; % Must be in MHZ
    if ~useHFsel
        
        ak = exp(-2*((feff-Exp.mwFreq*1e3)./Exp.ExciteWidth).^2); % all in MHz
    else
        % dE (epr) = nu_S + Sum(a*MIi)
        nTrans = prod(mult);
        if length(mult)==1,
            MIs(:, 1) = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
        else
            prod2 = prod(mult(2:end));
            tA = (-(mult(1)-1)/2 : 1 : (mult(1)-1)/2)';
            MIs(:, 1) = kron(tA, ones(prod2, 1));

            for ii = 2:(length(mult)-1)
                tA = (-(mult(ii)-1)/2 : 1 : (mult(ii)-1)/2)';
                prod1 = prod(mult(1:(ii-1)));
                ttA = kron(ones(prod1, 1), tA);
                prod2 = prod(mult((ii+1):end));
                MIs(:, ii) = kron(ttA, ones(prod2, 1));
            end
            tA = (-(mult(end)-1)/2 : 1 : (mult(end)-1)/2)';
            prod2 = prod(mult(1:(end-1)));
            MIs(:, length(mult)) = kron(ones(prod2, 1), tA);
       end
        ak = zeros(nangles, 1);
        for ii = 1:nTrans
            ffre = zeros(nangles, 1);
            for jj = 1:length(mult)
                ffre = ffre + projA(:, jj)*MIs(ii, jj);
            end
            ak = ak + exp(-2*((feff+ffre-Exp.mwFreq*1e3)./Exp.ExciteWidth).^2); % all in MHz
        end
        
    end

else
    ak = ones(nangles, 1);
end
%         Sy.Nucs = '57Fe';
%         Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I'};
%         Sy = Sys;
%         for kk = 1:length(isfield_num)
%             pk = isfield_num(kk);
%             tPar = getfield(Sys, Parameters{pk});
%             rmfield(Sy, Parameters{pk})
%             Sy = setfield(Sy, Parameters{pk}, tPar(1, :));
%         end
%         ak = orisel(Sy, Exp, [phi.', theta.']);    
if max(ak) < 1e-4, error('_orisel: no resonance'); end 
weights = ak.*ww;
weights = weights/max(weights);
% Weights = orisel(Sys,Par,Ori)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sketch plot for HF defence sequences
% --------------------------------------------------------------------
function sketchhfdec(Sys, Exp, Opt, Shift ,handles)
if isfield(Exp, 'phi')& isfield(Exp, 'theta'),
    nphi = length(Exp.phi);
    ntheta = length(Exp.theta);

    theta = reshape(Exp.theta(ones(nphi, 1), :), nphi*ntheta, 1)*pi/180;
    phi = reshape(Exp.phi(ones(ntheta, 1), :).', nphi*ntheta, 1)*pi/180;
else
    Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
    Opt.nKnots = safeget(Opt, 'nKnots', 20);
    [phi, theta] = sphgrid(Opt.Symmetry, Opt.nKnots);
    phi = phi.'; theta = theta.';
end
% calculating orientation selection
ctheta = cos(theta);
stheta = sin(theta);
cphi = cos(phi);
sphi = sin(phi);
angle(:, 1)=stheta.*cphi;
angle(:, 2)=stheta.*sphi;
angle(:, 3) = ctheta;

Gx=angle(:, 1)*Sys.g(1);
Gy=angle(:, 2)*Sys.g(2);
Gz=angle(:, 3)*Sys.g(3);
geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq, disp('hycalc: B0 = 0'); disperr(handles); return; end
   Exp.Field = (Exp.mwFreq*1e+9*planck)/(Exp.gField*bmagn*1e-3);  % mT
end
if (Exp.ExciteWidth > 0) & (Exp.mwFreq>0)
    gcenter = fld2g(Exp.Field*1E-3, Exp.mwFreq*1E9);
    gw      = gcenter*Exp.ExciteWidth*1e6*(planck)/(2*bmagn*1e-3)/Exp.Field /sqrt(2*log(2));
    ak = exp(-2*((geff-gcenter)./gw).^2);
else
    ak = ones(1,ntheta*nphi);
end
if max(ak) < 1e-3, disp('no resonanse'); disperr(handles); return; end 
% Weights = orisel(Sys,Par,Ori)

Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad

if isfield(Sys, 'Aiso'),
    Sys.A = safeget(Sys, 'A', 0) + Sys.Aiso;
end

if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end

Exp.tau = safeget(Exp, 'tau', 0);
Opt.Treshold = safeget(Opt, 'Treshold', 1e-3);
w1 = []; w2 = []; w1d = []; w2d = [];
num = get(handles.popSel, 'Value');
% Calculation
switch Sys.I
   case 0.5
      handles.sim{num}.x1 = [];
      handles.sim{num}.x2 = [];
      handles.sim{num}.y  = [];      
      Syst = Sys;
      Expt = Exp;
      Syst.A = Syst.A*safeget(Exp, 'Decoupling', 0);
      for cth =1:length(theta)
            if ak(cth)>Opt.Treshold,       
               freq = kv_hscfreq(Sys, Exp, phi(cth), theta(cth));
               w1(end+1, :) = freq(1); % MHz
               w2(end+1, :) = freq(2); % MHz
               
               freqdec = kv_hscfreq(Syst, Expt, phi(cth), theta(cth));
               w1d(end+1, :) = freqdec(1); % MHz
               w2d(end+1, :) = freqdec(2); % MHz               
            end   
      end
%          handles.sim{num}.x1 = [w1 w2]; % ESEEM one
%          handles.sim{num}.x2 = [w1d w2d]; % check ????????
%          handles.sim{num}.y  = [1, 1];

   case 1
      handles.sim{num}.x1 = [];
      handles.sim{num}.x2 = [];
      handles.sim{num}.y  = [];      
      Syst = Sys;
      Expt = Exp;
      Syst.A = Syst.A*0;
      for cth =1:length(theta)
            if ak(cth)>Opt.Treshold,       
               freq = kv_hscfreq(Sys, Exp, phi(cth), theta(cth));
               w1(end+1, :) = freq(:, 1)'; % MHz
               w2(end+1, :) = freq(:, 2)'; % MHz
               
               freqdec = kv_hscfreq(Syst, Expt, phi(cth), theta(cth));
               w1d(end+1, :) = freqdec(:, 1)'; % MHz
               w2d(end+1, :) = freqdec(:, 2)'; % MHz               
            end   
      end
    otherwise
      disp('hycalc: only I = 1/2 and I = 1 available!!!!')
      return;
   end     
      if ~isfield(Opt, 'ShowCor'), Opt.ShowCor = [1, 1; 1, 2; 2, 1; 2, 2]; end
      
      if Opt.ShowCor(1) < 2; % single freq 1nu
           handles.sim{num}.x1 = [w1, w2];
           handles.sim{num}.x2 = [w1d w2d];
           handles.sim{num}.y = [1, 1];
           if Opt.ShowCor(1)==0 % 1nu + 2nu
               handles.sim{num}.x1 = [handles.sim{num}.x1, w1, w2];
               handles.sim{num}.x2 = [handles.sim{num}.x2, w1d*2, w2d*2];
               handles.sim{num}.y = [1, 1, 1, 1];
           end
       else % only double frequencies (new HF decoupled defence) 
               handles.sim{num}.x1 = [w1, w2];
               handles.sim{num}.x2 = [w1d*2, w2d*2];
               handles.sim{num}.y = [1, 1];
       end

          
%           tRow = (1:(2*Sys.I+1));
%           Opt.ShowCor = [reshape(tRow(ones((2*Sys.I+1), 1), :),(2*Sys.I+1)^2, 1), reshape(tRow(ones((2*Sys.I+1), 1), :).',(2*Sys.I+1)^2, 1)];
%       end % Show all correlations
      
%       for cc = 1:size(Opt.ShowCor, 1)   
%          handles.sim{num}.x1 = [handles.sim{num}.x1 w1(:, Opt.ShowCor(cc, 1)) w2(:, Opt.ShowCor(cc, 2))];
%          handles.sim{num}.x2 = [handles.sim{num}.x2 w2(:, Opt.ShowCor(cc, 2)) w1(:, Opt.ShowCor(cc, 1))];
%          handles.sim{num}.y = [handles.sim{num}.y 1, 1];            
%       end
      
%                handles.sim{num}.x1 = [w1 w2];
%                handles.sim{num}.x2 = [(w2+w1)/2 (w2+w1)/2];
%                handles.sim{num}.y = [1, 1];

handles.sim{num}.amp = abs(max(handles.sim{num}.y, [], 2)-min(handles.sim{num}.y, [], 2)).';
if isfield(Opt, 'LineStyle')
   handles.sim{num}.LineStyle = Opt.LineStyle;
end
if isfield(Opt, 'LineWidth')
   handles.sim{num}.LineWidth = Opt.LineWidth;
end
if isfield(Opt, 'Marker')
   handles.sim{num}.Marker = Opt.Marker;
end
% if max(handles.sim{num}.amp)==0
%     error('All spectra have zero amplitude');
% else
%     handles.sim{num}.y = handles.sim{num}.y./handles.sim{num}.amp(ones(size(handles.sim{num}.y))); 
% end

guidata(handles.figure1, handles);
dispnorm(handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sketch plot for Field dependent ESE
% --------------------------------------------------------------------
function sketch2DeseB(Sys, Exp, Opt, Shift ,handles)
% if isfield(Exp, 'phi')& isfield(Exp, 'theta'),
%     nphi = length(Exp.phi);
%     ntheta = length(Exp.theta);
% 
%     theta = reshape(Exp.theta(ones(nphi, 1), :), nphi*ntheta, 1)*pi/180;
%     phi = reshape(Exp.phi(ones(ntheta, 1), :).', nphi*ntheta, 1)*pi/180;
%     nangles = nphi*ntheta
% else
%     Opt.Symmetry = safeget(Opt, 'Symmetry', 'Ci');
     Opt.nKnots = safeget(Opt, 'nKnots', 20);
%     [phi, theta] = sphgrid(Opt.Symmetry, Opt.nKnots);
%     phi = phi.'; theta = theta.';
%     nangles = length(phi);
% end
% % calculating orientation selection
% ctheta = cos(theta);
% stheta = sin(theta);
% cphi = cos(phi);
% sphi = sin(phi);
% angle(:, 1)=stheta.*cphi;
% angle(:, 2)=stheta.*sphi;
% angle(:, 3) = ctheta;
% 
% Gx=angle(:, 1)*Sys.g(1);
% Gy=angle(:, 2)*Sys.g(2);
% Gz=angle(:, 3)*Sys.g(3);
% geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0);
Exp.mwFreq = safeget(Exp, 'mwFreq', 0);
if ((length(safeget(Exp, 'gField', 0))<2)&(length(safeget(Exp, 'Field', 0))<2))
    warning('There is only one field point specified. Check "Exp.Field" or "Exp.gField"');
end
if isfield(Exp, 'gField'),
   if ~Exp.mwFreq, disp('hycalc: B0 = 0'); disperr(handles); return; end
   Exp.Field(1) = (Exp.mwFreq*1e+9*planck)/(Exp.gField(1)*bmagn*1e-3);  % mT
   Exp.Field(2) = (Exp.mwFreq*1e+9*planck)/(Exp.gField(2)*bmagn*1e-3);  % mT
end

if ~isfield(Exp, 'nFPoints')
    disp('no Exp.nFPoints have been found. Default number of field points [10] is taken');
end
nFPoints = safeget(Exp, 'nFPoints', 10);

axField = linspace(Exp.Field(1), Exp.Field(2), nFPoints);
Opt.Treshold = safeget(Opt, 'Treshold', 1e-3);
phi = zeros(nFPoints, Opt.nKnots*2);
theta = phi;
for ii = 1:nFPoints,
    gcenter = fld2g(axField(ii)*1E-3, Exp.mwFreq*1E9);
    gw      = gcenter*Exp.ExciteWidth*1e6*(planck)/(2*bmagn*1e-3)/axField(ii) /sqrt(2*log(2));
    [ph, th] = anorisel(Sys.g, gcenter, Opt.nKnots, 'Ci');
    phi(ii, :) = ph';
    theta(ii, :) = th';
    [phb1, thb1] = anorisel(Sys.g, gcenter+gw, Opt.nKnots, 'Ci');
    phib1(ii, :) = phb1';
    thetab1(ii, :) = thb1';
    [phb2, thb2] = anorisel(Sys.g, gcenter-gw, Opt.nKnots, 'Ci');
    phib2(ii, :) = phb2';
    thetab2(ii, :) = thb2';    
%     if (Exp.ExciteWidth > 0) & (Exp.mwFreq>0)
%         gcenter = fld2g(axField(ii)*1E-3, Exp.mwFreq*1E9);
%         gw      = gcenter*Exp.ExciteWidth*1e6*(planck)/(2*bmagn*1e-3)/axField(ii) /sqrt(2*log(2));
%         ak = exp(-2*((geff-gcenter)./gw).^2);
%         [val(ii),indmax(ii)] = min(abs(geff-gcenter));
%         [val1(ii),indbrd1(ii)] = min(abs(geff-gcenter-gw));
%         [val2(ii),indbrd2(ii)] = min(abs(geff-gcenter+gw));        
% %        indx{ii} = find(ak>Opt.Treshold);
%     else
%         indmax(ii) = floor(nangles/2);
%         indbrd1(ii) = 1;
%         indbrd2(ii) = nangles;
%     end
%     tstr = num2str(axField(ii));
%     if (val(ii)<Opt.Treshold), 
%         warning(['No resonanse at Bo = ',tstr, ' mT.']); 
%         disperr(handles); 
%     end 
    % Weights = orisel(Sys,Par,Ori)
end

Sys.Apa = safeget(Sys, 'Apa', [0, 0, 0])*pi/180; % degree -> rad
Sys.Qpa = safeget(Sys, 'Qpa', [0, 0, 0])*pi/180; % degree -> rad

if isfield(Sys, 'Aiso'),
    Sys.A = safeget(Sys, 'A', 0) + Sys.Aiso;
end

if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end

Exp.tau = safeget(Exp, 'tau', 0);

w1 = []; w2 = []; w1d = []; w2d = [];
num = get(handles.popSel, 'Value');
% Calculation
handles.sim{num}.x1 = [];
handles.sim{num}.x2 = [];
handles.sim{num}.y  = [];

switch Sys.I
    case 0.5
        Syst = Sys;
        Expt = Exp;
        for cth =1:nFPoints
            %if ak(cth)>Opt.Treshold,
            Expt.Field = axField(cth);
            for ii = 1:length(phi(1, :))
                freq= kv_hscfreq(Sys, Expt, phi(cth, ii), theta(cth, ii));
                freqb1 = kv_hscfreq(Sys, Expt, phib1(cth, ii), thetab1(cth, ii));
                freqb2 = kv_hscfreq(Sys, Expt, phi(cth, ii), theta(cth, ii));     
                w1(cth, ii) = freq(1);
                w1b1(cth, ii) = freqb1(1);
                w1b2(cth, ii) = freqb2(1); % MHz
                w2(cth, ii) = freq(2);
                w2b1(cth, ii) = freqb1(2);
                w2b2(cth, ii) = freqb2(2); % MHz
            end

            %end
        end
%         handles.sim{num}.x1 = [w1 w2]; % ESEEM one
%         handles.sim{num}.x2 = axFreq; % check ????????
%         handles.sim{num}.y  = [1, 1];
        
    case 1
        Syst = Sys;
        Expt = Exp;
        w1 = zeros(nFPoints, Opt.nKnots*6);
        w2 = w1; w1b1 = w1; w1b2 = w1; w2b1 = w1; w2b2 = w1;
        for cth =1:nFPoints
            Expt.Field = axField(cth);
            %if ak(cth)>Opt.Treshold,   
            tw1 = []; tw2 = []; tw1b1 = []; tw2b1 = []; tw1b2 = []; tw2b2 = [];
            for ii = 1:length(phi(1, :))
                freq = kv_hscfreq(Sys, Expt, phi(cth, ii), theta(cth, ii));
                freqb1 = kv_hscfreq(Sys, Expt, phib1(cth, ii), thetab1(cth, ii));
                freqb2 = kv_hscfreq(Sys, Expt, phib2(cth, ii), thetab2(cth, ii)); 
                tw1 = [tw1, freq(:, 1)'];
                tw2 = [tw2, freq(:, 2)'];
                tw1b1 = [tw1b1, freqb1(:, 1)'];
                tw2b1 = [tw2b1, freqb1(:, 2)'];
                tw1b2 = [tw1b2, freqb2(:, 1)'];
                tw2b2 = [tw2b2, freqb2(:, 2)'];                
           end
            w1(cth, :) = tw1;
            w1b1(cth, :) = tw1b1;
            w1b2(cth, :) = tw1b2; % MHz
            w2(cth, :) = tw2;
            w2b1(cth, :) = tw2b1;
            w2b2(cth, :) = tw2b2; % MHz            
%                 w2b1(cth, :) = freqb1(:, 2)';
%                 w2b2(cth, :) = freqb2(:, 2)'; % MHz
            %freqdec = kv_hscfreq(Syst, Expt, phi(cth), theta(cth));
            %w1d(cth, :) = freqdec(:, 1)'; % MHz
            %w2d(cth, :) = freqdec(:, 2)'; % MHz               
            %end   
            
        end
end      


Opt.ShowCor = safeget(Opt, 'ShowCor', 0);
    handles.sim{num}.x2 = axField;
    handles.sim{num}.y = [1, 1];
switch safeget(Opt, 'ShowCor', 0)
    case 0 % only center line to show
        handles.sim{num}.x1 = [w1, w2];
        handles.sim{num}.y = [1, 1];
    case 1 % only borther to show
        handles.sim{num}.x1 = [w1b1, w1b2, w2b1, w2b2];
        handles.sim{num}.y = [1, 1, 1, 1];
    otherwise % show everything
        handles.sim{num}.x1 = [w1b1, w1, w1b2, w2b1, w2, w2b2];
        handles.sim{num}.y = [1, 1, 1, 1, 1, 1];
end

handles.sim{num}.amp = abs(max(handles.sim{num}.y, [], 2)-min(handles.sim{num}.y, [], 2)).';
if isfield(Opt, 'LineStyle')
   handles.sim{num}.LineStyle = Opt.LineStyle;
end
if isfield(Opt, 'LineWidth')
   handles.sim{num}.LineWidth = Opt.LineWidth;
end
if isfield(Opt, 'Marker')
   handles.sim{num}.Marker = Opt.Marker;
end
% if max(handles.sim{num}.amp)==0
%     error('All spectra have zero amplitude');
% else
%     handles.sim{num}.y = handles.sim{num}.y./handles.sim{num}.amp(ones(size(handles.sim{num}.y))); 
% end

guidata(handles.figure1, handles);
dispnorm(handles)
return

% --------------------------------------------------------------------
function plothis(handles)
handles = guidata(handles.figure1);
colors = handles.Color;
% stat = strcmp(get(handles.mFixAxis, 'Checked'), 'on');

hand = guidata(handles.hh);

hand.out = {};
% handles.SimSrc = safeget(handles, 'SimSrc', 1);
cnum = get(handles.popSel, 'Value');

% handles.OverlayOptions 
str = handles.sim{cnum}.Script;
for k=1:size(str,1),
    if ~isempty(strfind(str{k},'Shift.'))
        try,eval([str{k} ';']);end
    end
end
if ~safeget(Shift, 'Overlay', 0)
    handles.CalcType(cnum) = handles.OverlayOptions{cnum}(1);
%      handles.SimSrc = 0;
else
    handles.CalcType(cnum) = handles.OverlayOptions{cnum}(2);
%      handles.SimSrc = 1;
end
 handles.SimSrc = sum( handles.CalcType(cnum)==[1, 4]);
 
if handles.SimSrc
   ax = hand.src.ax;
   y = hand.src.y;
   
   if handles.process, [y, ax] = processing(y, ax);
   else
      ax.diff = '';
      ax.filt = '';
   end
   
   c = {'off', 'on'};
   %stat = strcmp(get(handles.mShsumm, 'Checked'), 'on');
   %bdiff = strcmp(get(handles.mShdiff, 'Checked'), 'on');
   %basln = strcmp(get(handles.mSimSubBl, 'Checked'), 'on');
   %autobasln = strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on');
   hipx =0;  hipy = 0;
   errr = 0;
   
   hand.out{1}.ax = ax;
   hand.out{1}.y = y;
   hand.out{1}.ax.c = [0, 0, 1];
   hand.out{1}.ax.s = 0;
   hand.out{1}.xlabel = ax.xlabel;
   hand.out{1}.forceplot = 'image';
   hand.out{1}.title = 'Src';
end


if handles.CalcType(cnum) == 1 % sketch simulation
    for kk = 1:handles.num
        if isempty(handles.sim{kk}.x1) continue; end
        if handles.CalcType(kk)~=1  continue; end % skip "fd" simmulations
        Shift = struct('x',0,'y',0,'Scale',1,'wf',1);
        % NOTE: here parapeter Shift.Scale response for rescaling the axes.
        % If Shift.Scale = 1 -> [MHz]
        % load script 
        str = handles.sim{kk}.Script;
        for k=1:size(str,1),
            if ~isempty(strfind(str{k},'Shift.')) 
                try,eval([str{k} ';']);end
            end
        end
        if strcmp(handles.sim{kk}.stype, 'sim')
            data_mult = exp(handles.div/10);
        else
            data_mult = 1;
        end
        alpha = [1, -1, 1, -1];
        beta = [1, 1, -1, -1]; 
        if length(Shift.Scale) <2,  Shift.Scale(2) = Shift.Scale(1); end
        for as = 1:4
            if handles.Quad(as)
                x1Temp = handles.sim{kk}.x1*alpha(as)*Shift.Scale(1);
                x2Temp = handles.sim{kk}.x2*beta(as)*Shift.Scale(2);
                yTemp = handles.sim{kk}.y;
                %hand.out{end + 1}.ax = ax;
                hand.out{end+1}.ax.dx = 0;
                hand.out{end}.ax.x = x1Temp + Shift.x1;
                hand.out{end}.ax.y = x2Temp + Shift.x2;
                hand.out{end}.y = yTemp*data_mult;
                hand.out{end}.ax.Color = colors(kk, :);
                hand.out{end}.ax.LineStyle = safeget(handles.sim{kk}, 'LineStyle', '-');
                hand.out{end}.ax.LineWidth = safeget(handles.sim{kk}, 'LineWidth', 1);
                hand.out{end}.ax.Marker = safeget(handles.sim{kk}, 'Marker', 'none');
                hand.out{end}.ax.s = Shift.y*data_mult ; % Shift in relative units
                %hand.out{end}.forceplot = 'line';
                hand.out{end}.title = ['Out', num2str(kk, 3)];
            end
        end
    end
    %script = lscript(handles);        
    hand.PlotType = 5;
elseif handles.CalcType(cnum) == 2 % fd
    Shift = struct('x',0,'y',0,'Scale',1,'wf',1);
    % load script 
    str = handles.sim{cnum}.Script;
    cc = 1;
    if handles.sepstor, cc = cnum; end
    for k=1:size(str,1),
        if ~isempty(strfind(str{k},'Shift.')) 
            try,eval([str{k} ';']);end
        end
    end
    % What to draw VVVVVV
    if length(Shift.Scale) <2,  Shift.Scale(2) = Shift.Scale(1); end
    x1size = length(handles.sim{cc}.x1);
    x2size = length(handles.sim{cc}.x2);
    PPMM = [1, 0, 1, 0;...
            1, 1, 0, 0;...
            0, 1, 0, 1;...
            0, 0, 1, 1];    
    SP = [0; 0; 0; 0];
    for as = 1:4
        if handles.Quad(as)
               SP(:) = PPMM(:, as) + SP(:);
        end
    end
    XAxD = [];
    YAxD = [];
    if SP(3), XAxD = [XAxD, 1:floor(x1size/2+0.5);]; end
    if SP(4), YAxD = [YAxD, 1:floor(x2size/2+0.5);]; end
    if SP(1), XAxD = [XAxD, (floor(x1size/2+0.5)+1):x1size;]; end
    if SP(2), YAxD = [YAxD, (floor(x2size/2+0.5)+1):x2size;]; end    
    % AAAAAAAAAAAAAAAA
    %hand.out{end + 1}.ax = ax;
    hand.out{1}.ax.dx = 0;
    hand.out{1}.ax.x = handles.sim{cc}.x1(XAxD, 1)*Shift.Scale(1);
    hand.out{1}.ax.y = handles.sim{cc}.x2(YAxD, 1)*Shift.Scale(2);
    hand.out{1}.ax.xlabel = 'Frequency, MHz';
    hand.out{1}.ax.ylabel = 'Frequency, MHz';    
    hand.out{1}.y = handles.sim{cc}.y(XAxD, YAxD);
    hand.out{1}.ax.Color = colors(1, :);
    hand.out{1}.ax.LineStyle = safeget(handles.sim{cc}, 'LineStyle', '-');
    hand.out{1}.ax.LineWidth = safeget(handles.sim{cc}, 'LineWidth', 1);
    hand.out{1}.ax.Marker = safeget(handles.sim{cc}, 'Marker', 'none');
    hand.out{1}.ax.s = 0; % Shift in relative units
    hand.out{1}.Script = handles.sim{cnum}.Script;
    hand.out{1}.ax.title = 'HYSCORE simulation';

    if isfield(Shift, 'contour')
        hand.out{1}.ax.contour = Shift.contour;
    end
    %hand.out{1}.forceplot = 'line';
    hand.out{1}.title = ['Result'];    
    hand.PlotType = 4;
elseif handles.CalcType(cnum) == 3 % ESEEMvsBO fd
    Shift = struct('x',0,'y',0,'Scale',1,'wf',1);
    % load script 
    str = handles.sim{cnum}.Script;
    cc = 1;
    if handles.sepstor, cc = cnum; end
    for k=1:size(str,1),
        if ~isempty(strfind(str{k},'Shift.')) 
            try,eval([str{k} ';']);end
        end
    end
    % What to draw VVVVVV
    if length(Shift.Scale) <2,  Shift.Scale(2) = Shift.Scale(1); end
    x1size = length(handles.sim{cc}.x1);
    x2size = length(handles.sim{cc}.x2);
    
     % AAAAAAAAAAAAAAAA
    %hand.out{end + 1}.ax = ax;
    hand.out{1}.ax.dx = 0;
    hand.out{1}.ax.x = handles.sim{cc}.x1*Shift.Scale(1);
    hand.out{1}.ax.y = handles.sim{cc}.x2*Shift.Scale(2);
    hand.out{1}.ax.xlabel = 'Frequency, MHz';
    hand.out{1}.ax.ylabel = 'Frequency, MHz';    
    hand.out{1}.y = handles.sim{cc}.y;
    hand.out{1}.ax.Color = colors(1, :);
    hand.out{1}.ax.LineStyle = safeget(handles.sim{cc}, 'LineStyle', '-');
    hand.out{1}.ax.LineWidth = safeget(handles.sim{cc}, 'LineWidth', 1);
    hand.out{1}.ax.Marker = safeget(handles.sim{cc}, 'Marker', 'none');
    hand.out{1}.ax.s = 0; % Shift in relative units
    hand.out{1}.Script = handles.sim{cnum}.Script;
    hand.out{1}.ax.title = 'HYSCORE simulation';
    if isfield(Shift, 'contour')
        hand.out{1}.ax.contour = Shift.contour;
    end
    %hand.out{1}.forceplot = 'line';
    hand.out{1}.title = ['Result'];    
    hand.PlotType = 4;    
elseif handles.CalcType(cnum) == 4 % ESEEMvsBO fd image + contour
    for kk = 1:handles.num
        Shift = struct('x',0,'y',0,'Scale',1,'wf',1);
        % load script 
        str = handles.sim{kk}.Script;
        cc = kk;
        for k=1:size(str,1),
            if ~isempty(strfind(str{k},'Shift.')) 
                try,eval([str{k} ';']);end
            end
        end
        % What to draw VVVVVV
        if length(Shift.Scale) <2,  Shift.Scale(2) = Shift.Scale(1); end
        x1size = length(handles.sim{cc}.x1);
        x2size = length(handles.sim{cc}.x2);
        
        % AAAAAAAAAAAAAAAA
        %hand.out{end + 1}.ax = ax;
        hand.out{end+1}.ax.dx = 0;
        hand.out{end}.ax.x = handles.sim{cc}.x1*Shift.Scale(1);
        hand.out{end}.ax.y = handles.sim{cc}.x2*Shift.Scale(2);
        hand.out{end}.ax.xlabel = 'Frequency, MHz';
        hand.out{end}.ax.ylabel = 'Frequency, MHz';    
        hand.out{end}.y = handles.sim{cc}.y;
        hand.out{end}.ax.Color = colors(kk, :);
        hand.out{end}.ax.LineStyle = safeget(handles.sim{cc}, 'LineStyle', '-');
        hand.out{end}.ax.LineWidth = safeget(handles.sim{cc}, 'LineWidth', 1);
        hand.out{end}.ax.Marker = safeget(handles.sim{cc}, 'Marker', 'none');
        hand.out{end}.ax.s = 0; % Shift in relative units
        hand.out{end}.Script = handles.sim{cnum}.Script;
        hand.out{end}.ax.title = 'HYSCORE simulation';
        hand.out{end}.forceplot = 'contour';
        if isfield(Shift, 'contour'),
            hand.out{end}.ax.contour = Shift.contour;
        end
        %hand.out{1}.forceplot = 'line';
        hand.out{end}.title = ['Result ', num2str(kk)];    
    end
    %script = lscript(handles);        
    hand.PlotType = 5;
end
guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function handles = ZZ_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

switch hObject
    case handles.mSimAutoZoom
    c = {'off', 'on'};
    stat = strcmp(get(handles.mSimAutoZoom, 'Checked'), 'on');
    set(handles.mSimAutoZoom, 'Checked', [c{~stat + 1}]);
otherwise
    num = get(handles.popSel, 'Value');
    max1 = max(max(real(handles.sim{1}.y)));
    max2 = max(real(hand.src.y));
    if hObject==handles.ZZ
    else
        amp1 = abs(max1-min(min(real(handles.sim{1}.y))));
        amp2 = abs(max2-min(real(hand.src.y)));
    end
    a = max(10*log(amp2/amp1));
    if isempty(a), a=1; end
    set(handles.eResize, 'String', num2str(a));
    handles.div = a;
    guidata(handles.figure1, handles);
    plothis(handles);
end  

% --------------------------------------------------------------------
function mLoad_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);

full = get(handles.figure1, 'FileName');
[fpath,name,ext] = fileparts(full);

if ischar(get(handles.popSel, 'String')), 
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end
val = get(handles.popSel, 'Value');

switch hObject
case {handles.mLoad, handles.ptChange}
    [fname,pname] = uigetfile([fpath '\pepper_scr\*.m'],'Open Script');
    if fname == 0, return; end
    
    [name, script] = LoadScript([pname '\' fname], handles.sim{val}.stype);
    set(handles.list, 'String', script, 'Value', 1);
    
    str{val} = [num2str(val), ': ', name];
    set(handles.popSel, 'String', str);
    handles.sim{val}.Script = script;
    guidata(handles.figure1, handles);
    list_Callback(hObject, [], handles);
case handles.mLoadAll
    if exist('uigetdir')==5
        directory_name = [uigetdir(path,'Select directory'), '\'];
        
    else
        [a, directory_name] = uiputfile([fpath,'\select.dir'],'Select directory');
    end
    if ~directory_name, return; end
    
    files = dir([directory_name,'*.m']);
    scripts = {files.name};
    
    hand = guidata(handles.hh);
    stype = 'sim';
    for kk=1:size(scripts, 2)
        [name, script] = LoadScript(fullfile(directory_name, scripts{kk}), 'sim');
        str{handles.num + 1} = [num2str(handles.num + 1), ': ', name];
        handles.sim{handles.num + 1}.Script = script;
        handles.sim{handles.num + 1}.x = hand.src.ax.x;
        handles.sim{handles.num + 1}.y = zeros(size(hand.src.ax.x));
        handles.sim{handles.num + 1}.x1 = 0;
        handles.sim{handles.num + 1}.stype = stype;
        handles.Color(handles.num + 1, :) = handles.Color(mod(handles.num,7)+1, :);
        handles.num = handles.num + 1;
    end
    
    str = chngnum(str); % changing numbers in popSel list
    set(handles.list, 'String', script, 'Value', 1);
    set(handles.popSel, 'String', str, 'Value', handles.num);
    guidata(handles.figure1, handles);
    popSel_Callback(hObject, eventdata, handles);
    
end

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
Shift.a=[];
for k=1:size(str,2),
    try,eval([str{k} ';']);end
end
% result mofification
if ~isfield(Shift, 'x1'), str{end+1} = 'Shift.x1 = 0';end
if ~isfield(Shift, 'x2'), str{end+1} = 'Shift.x2 = 0';end
if nargout > 0 , 
    varargout = {fname , str.'};
end

% --------------------------------------------------------------------
function mSave_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

formfile = get(hand.MainFigure, 'FileName');
[path, name,e] = fileparts(formfile); 

switch hObject
case {handles.mSave, handles.ptSave}
    num = get(handles.popSel, 'Value');
    sims = get(handles.popSel, 'String');
    if ischar(sims), 
        string = sims;     
    else
        string = sims{num};
    end
    
    [fname,pname] = uiputfile([path, '\', 'pepper_scr', '\',string(3:end), '.m'],'Save Script');
    if fname == 0, return; end
    
    [fpath,fname,ext] = fileparts([pname, fname]); 
    
    if strcmp(ext ,'') , ext = '.m'; end
    fid = fopen([pname, fname, ext], 'w');
    str = get(handles.list, 'String');
    for kk = 1:size(str, 1)
        fprintf(fid, '%s\r\n', str{kk});
    end
    fclose(fid);
case  handles.ptClipboard
    num = get(handles.popSel, 'Value');
    sims = get(handles.popSel, 'String');
    if ischar(sims), string = sims;     
    else, string = sims{num};
    end
    Scr = get(handles.list, 'String');
    out = [];
    for kk = 1:length(Scr)
        out = sprintf( '%s%s\n', out, Scr{kk});
    end
    clipboard('copy',out);
    disp(['HYSCORE: script ''',sims{num},''' copied to the clipboard.']);    
case handles.mSaveAll
    if exist('uigetdir')==5
        directory_name = uigetdir(path,'Select directory');
    else
        [a, directory_name] = uiputfile([path,'\select.dir'],'Select directory');
    end
    if ~directory_name, return; end
    
    sims = get(handles.popSel, 'String');
    if ischar(sims), 
        cnum = 1;
    else
        cnum = 1:length(sims);
    end
    
    for num = cnum
        if ischar(sims), 
            [tmp, string] = strtokstr(sims, ': ');     
        else
            [tmp, string] = strtokstr(sims{num}, ': ');
        end
        [pa,name,ext] = fileparts(string);
        fid = fopen([directory_name,'\', name, '.m'], 'w');
        str = handles.sim{num}.Script;
        for kk = 1:length(str)
            fprintf(fid, '%s\r\n', str{kk});
        end
        fclose(fid); 
        disp(['Simplugin: simulation script ''', name,''' have been saved to ', directory_name]);
    end
   
end
% --------------------------------------------------------------------
% --- Executes on button press in check.
function check_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
c = {'off', 'on'};
col = {[1 0.7 0.7], [0.7, 1, 0.7]};
stat = strcmp(get(handles.mCalc, 'Checked'), 'on');
set(handles.mCalc, 'Checked', [c{~stat + 1}]);
set(handles.butCalc, 'Value', ~stat,'BackgroundColor', [col{~stat + 1}]);

if ~stat, 
    %   list_Callback(hObject, eventdata, handles)
    ePars_Callback(hObject, eventdata, handles)
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(handles.figure1);

% ------------------------------------------------------------
function popSel_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
num = get(handles.popSel, 'Value');
set(handles.list, 'String', handles.sim{num}.Script, 'Value', 2);
list_Callback(hObject, eventdata, handles);
%set(handles.eResize, 'String', num2str(handles.div));
set(handles.butCol, 'BackgroundColor', handles.Color(num, :));


% ------------------------------------------------------------
function butCol_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
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
handles = guidata(handles.figure1);
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
outstr = chngnum({str{ss, :}}');
set(handles.popSel, 'String', outstr, 'Value', 1);
guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function AddSim_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);

    stype = 'sim';
    title = 'Open Simulation Script';


full = get(handles.figure1, 'FileName');
[fpath,name,ext] = fileparts(full); 
[fname, pname] = uigetfile([fpath '\pepper_scr\*.m'],title);
if fname == 0, return; end

[name, script] = LoadScript([pname fname], stype);
set(handles.list, 'String', script, 'Value', 1);

if ischar(get(handles.popSel, 'String')),
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end

str = chngnum(str); % changing numbers in popSel list

str{handles.num + 1} = [num2str(handles.num + 1), ': ', name];
set(handles.popSel, 'String', str, 'Value', handles.num + 1);

hand = guidata(handles.hh);
handles.sim{handles.num + 1}.Script = script;
handles.sim{handles.num + 1}.x1 = 0;
handles.sim{handles.num + 1}.x2 = 0;
handles.sim{handles.num + 1}.y = 0;
handles.sim{handles.num + 1}.stype = stype;
handles.Color(handles.num + 1, :) = handles.Color(mod(handles.num,7)+1, :);
handles.num = handles.num + 1;

guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function varargout = mFileCopy_Callback(h, eventdata, handles, varargin)
num = get(handles.popSel, 'Value');

script = get(handles.list, 'String');
set(handles.list, 'Value', 2);

if ischar(get(handles.popSel, 'String')),
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end

str{handles.num + 1} = [str{num}, ' copy'];

str = chngnum(str);

set(handles.popSel, 'String', str, 'Value', handles.num + 1);

hand = guidata(handles.hh);
handles.sim{handles.num + 1}.Script = script;
handles.sim{handles.num + 1}.x1 = 0;
handles.sim{handles.num + 1}.x2 = 0;
handles.sim{handles.num + 1}.y = 0;
handles.sim{handles.num + 1}.stype = handles.sim{num}.stype;
if handles.num + 1 > 7
    handles.Color(handles.num + 1, :) = handles.Color(handles.num-7, :);
end
handles.num = handles.num + 1;

guidata(handles.figure1, handles);
popSel_Callback(h, eventdata, handles);

% --------------------------------------------------------------------
function mSimCol_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
num = get(handles.popSel, 'Value');
switch h
case handles.mColSum
    c = uisetcolor(handles.SumCol);
    if size(c, 2) > 1, handles.SumCol = c; 
    else, return;
    end
case handles.mColDiff
    c = uisetcolor(handles.DiffCol);
    if size(c, 2) > 1, handles.DiffCol = c;
    else, return;
    end        
end
guidata(handles.figure1, handles);
plothis(handles);


% --------------------------------------------------------------------
function mFixAxis_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
c = {'off', 'on'};
stat = strcmp(get(handles.mFixAxis, 'Checked'), 'on');
set(handles.mFixAxis, 'Checked', [c{~stat + 1}]);
if ~stat, 
    plothis(handles);
end

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

% --- Executes on button press in butCalc.
function butCalc_Callback(hObject, eventdata, handles)
check_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function script = lscript(handles)
disp('lscript: not yet relized');

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
    % Is token a new? /alsi
    if sum(strcmp(tokens, tok)), continue; end;
    % Is token a number ?
    if isempty(str2num(tok))
        tokens{end+1} = tok;
    end
end

% --------------------------------------------------------------------
function outstr = chngnum(str)
for k = 1:size(str, 1)
    tempstr = str{k};
    if length(tempstr)>2,
        [strstart, strend] = strtokstr(tempstr, ': ');
    else
        strend = tempstr;
    end
    if ~isempty(strend), str{k} = [num2str(k), ': ', strend];
    else
        str{k} = [num2str(k), ': ', strstart];
    end
end
outstr = str;

% --------------------------------------------------------------------
function mSimAgr_Callback(h, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimAgr, 'Checked'), 'on');
set(handles.mSimAgr, 'Checked', str{~stat + 1});
handles.process = ~stat;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function mFileAddPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');
answer = kv_constructfield;

if isempty(answer) , return; end

scrlen = length(scr);
script = {};
for count = 1:pos
    script{end+1} = scr{count};
end
script{end+1} = [answer];
if pos ~= scrlen,
    for count = pos+1:scrlen
        script{end+1} = scr{count};
    end
end

set(handles.list, 'String', script, 'Value', pos+1);
handles.sim{num}.Script = script;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);
% --------------------------------------------------------------------
function mFileChPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');

answer = kv_constructfield(scr{pos});
if isempty(answer) , return; end

scr{pos} = answer;

set(handles.list, 'String', scr, 'Value', pos);
handles.sim{num}.Script = scr;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);

% --------------------------------------------------------------------
function mFileDelPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');

button = questdlg('Do you want to continue?',...
    'Continue','Yes','No','Yes');
if strcmp(button,'No'), return; end

scrlen = length(scr);
script = {};
for count = 1:pos-1,
    script{end + 1} = scr{count};
end
if pos~=scrlen,
    for count = pos+1:scrlen
        script{end + 1} = scr{count};        
    end
end

set(handles.list, 'String', script, 'Value', pos);
handles.sim{num}.Script = script;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function mFileChName_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.popSel, 'String');
if ischar(scr), scr = {scr}; end

num = get(handles.popSel, 'Value');

prompt = {'Type name'};
dlg_title = 'Change simulation name';
num_lines= 1;
if length(scr{num}) > 2, 
    [en, st] = strtokstr(scr{num}, ': ');
else
    st = scr{num};
end

def     = {st};
answer  = inputdlg(prompt,dlg_title,num_lines,def);

scr{num} = [num2str(num), ': ',answer{1}];
set(handles.popSel, 'String', scr, 'Value', num);

% --------------------------------------------------------------------
function mAbout_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
    [fpath,n,e] = fileparts(which('kazan'));
    filename = [fpath '\help\hyscorebox.html'];
    %ini = hthelp(filename);
    stat = web(filename);
% --------------------------------------------------------------------
function varargout = mSimApplyFreq_Callback(h, eventdata, handles, varargin)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimApplyFreq, 'Checked'), 'on');
set(handles.mSimApplyFreq, 'Checked', str{~stat + 1});
handles.applyfreq = ~stat;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function UpdateScriptLine(handles, fld, newval)
hand = guidata(handles.hh);
handles = guidata(handles.figure1);
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');
for kk=1:length(Res)
    if ~isempty(strfind(Res{kk}, fld))
        Res{kk}=[fld,' = ',num2str(newval)];
        set(handles.list, 'String',Res);
        break;
    end
end

% --------------------------------------------------------------------
function res = lin(x1, x2, N)
% "linspace" like function
res = linspace(x1, x2, N);
return

% --------------------------------------------------------------------
function checkQ_Callback(h, eventdata, handles)
Val = get(h, 'Value');
switch h
    case handles.chPlPl
        handles.Quad(1) = Val;
    case handles.chPlMin
        handles.Quad(2) = Val;
    case handles.chMinPl
        handles.Quad(3) = Val;
    case handles.chMinMin
        handles.Quad(4) = Val;
end

guidata(handles.figure1, handles);
plothis(handles);
        
% --------------------------------------------------------------------
function mSimSrc_Callback(h, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimSrc, 'Checked'), 'on');
set(handles.mSimSrc, 'Checked', str{~stat + 1});
handles.process = ~stat;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function disperr(handles)
set(handles.figure1, 'Name', ['HYSCORE error']);
set(handles.butCalc, 'BackgroundColor',  [1, 0.7, 0.7]);
% --------------------------------------------------------------------
function dispnorm(handles)
set(handles.figure1, 'Name', ['HYSCORE']);
set(handles.butCalc, 'BackgroundColor',  [0.7, 1, 0.7]);

% --------------------------------------------------------------------
function mFileMove_Callback(hObject, eventdata, handles)
Str = get(handles.list, 'String');
Val = get(handles.list, 'Value');
if ischar(Str), return; end % ischar==1 means only 1 string in list.
maxval = length(Str);
Str1 = Str;
switch hObject
    case {handles.mFileMoveUp, handles.mMoveUp}
        if Val == 1, return; end
        Str1{Val-1} = Str{Val};
        Str1{Val} = Str{Val-1};        
        newval = Val-1;
    case {handles.mFileMoveDown, handles.mMoveDown}
        if Val == maxval, return; end
        Str1{Val+1} = Str{Val};
        Str1{Val} = Str{Val+1};                
        newval = Val+1;
end
set(handles.list, 'String', Str1, 'Value', newval);

% --------------------------------------------------------------------
function varargout = mFileOpenEditor_Callback(h, eventdata, handles, varargin)

num = get(handles.popSel, 'Value');
outscr = kv_scripteditor(get(handles.list,'String'));
% open(handles.sim{num}.ScriptName);
if isempty(outscr), return; end

set(handles.list, 'String', outscr, 'Value', 1);
handles.sim{get(handles.popSel, 'Value')}.Script = outscr;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);

% --------------------------------------------------------------------
function varargout = mSimSepSt_Callback(hObject, eventdata, handles, varargin)
handles = guidata(handles.figure1);
c = {'off', 'on'};
stat = strcmp(get(hObject, 'Checked'), 'on');
set(hObject, 'Checked', [c{~stat + 1}]);
handles.sepstor = c{~stat + 1};
guidata(handles.figure1, handles);