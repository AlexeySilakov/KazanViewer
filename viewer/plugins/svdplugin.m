function varargout = svdplugin(varargin)
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
    token = 'SVDPLUGIN_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,~,e] = fileparts(which('kazan'));
    inifilename = [fpath '\kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try opc = ini.KazanViewer.MultyKazan; end
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end

	fig = figure; %openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    set(fig, 'Units', 'pixels', 'Name', 'SVD plugin', 'NumberTitle', 'off');
    set(0, 'Units', 'pixels');
    vv = get( 0, 'ScreenSize');
    set(fig, 'Position', [319 100 250 min(800, vv(4)-200)], 'Tag', 'figure1', 'MenuBar', 'none'); %%% hide until everything is done
    set(fig, 'ResizeFcn', 'svdplugin(''figure1_ResizeFcn'',[],[],guidata(gcf))');
    
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  handles = loadGUI(handles);
  
  handles.hh = 0;
  handles.ploter = 0;
  handles.num = 1;
  handles.div = 0;
  handles.src = {};
%   handles.sim{1}.Script = get(handles.list, 'String');
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

%   set(handles.LoSc, 'Accelerator', 'L');
%   set(handles.SvSc, 'Accelerator', 'S');
%   set(handles.mCalc, 'Accelerator', 'C');
%   set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
  guidata(fig, handles);
%   list_Callback([], [], handles);
  figure1_ResizeFcn(fig, [], handles)
  
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FileUpdate(varargin);
return;

function handles = loadGUI(handles)
FName = get(0,'FixedWidthFontName'); FSize = 10;
handles.File = uimenu(handles.figure1, 'Tag', '', 'Label', 'File', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 

% handles.mFileDelete = uimenu(handles.File, 'Tag', 'mFileDelete', 'Label', 'Delete', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''butDel_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
% handles.LoSc = uimenu(handles.File, 'Tag', 'LoSc', 'Label', 'Load', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''LoSc_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'L' ); 
% handles.SvSc = uimenu(handles.File, 'Tag', 'SvSc', 'Label', 'Save', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''SvSc_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'S' ); 
% handles.mFileCopy = uimenu(handles.File, 'Tag', 'mFileCopy', 'Label', 'Copy', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mFileCopy_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
% handles.mFileScript = uimenu(handles.File, 'Tag', 'mFileScript', 'Label', 'Show script', 'Checked', 'off', 'Separator', 'on', 'Callback', 'fftplugin(''mFileScript_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
% handles.mFileExecScript = uimenu(handles.File, 'Tag', 'mFileExecScript', 'Label', 'Execute script', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mFileScript_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
% 
% handles.Options = uimenu(handles.figure1, 'Tag', '', 'Label', 'Options', 'Checked', 'off', 'Separator', 'off', 'Callback', '' ); 
% handles.mSimSubBl = uimenu(handles.Options, 'Tag', 'mSimSubBl', 'Label', 'Subtract baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mSimShow_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', '' ); 
% handles.mCalc = uimenu(handles.Options, 'Tag', 'mCalc', 'Label', 'Calculate', 'Checked', 'off', 'Separator', 'on', 'Callback', 'fftplugin(''check_Callback'',gcbo,[],guidata(gcbo))', 'Accelerator', 'C' ); 
% 
% handles.mAbout = uimenu(handles.figure1, 'Tag', 'mAbout', 'Label', 'About', 'Checked', 'off', 'Separator', 'off', 'Callback', 'fftplugin(''mAbout_Callback'',gcbo,[],guidata(gcbo))' ); 

handles.butLoad = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butLoad', 'String', 'LOAD DATA', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''butLoad_Callback'',gcbo,[],guidata(gcbo))'); 

handles.togSVD = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'togSVD', 'String', 'SVD', 'Value', 1, 'Units', 'pixels', 'Callback', 'svdplugin(''togMode_Callback'',gcbo,[],guidata(gcbo))'); 
handles.togGlo = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'togGlo', 'String', 'Glob.An.', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''togMode_Callback'',gcbo,[],guidata(gcbo))'); 
handles.togTar = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'togTar', 'String', 'Target', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''togMode_Callback'',gcbo,[],guidata(gcbo))'); 

handles.panelSVD =  uipanel(handles.figure1,  'Tag', 'panel1', 'Units', 'pixels',  'Visible', 'on', 'BackgroundColor', [0.95, 0.9, 0.9]); 
handles.panelGlo =  uipanel(handles.figure1,  'Tag', 'panel1',  'Units', 'pixels', 'Visible', 'off', 'BackgroundColor', [0.9, 1, 0.9]); 
handles.panelTar =  uipanel(handles.figure1,  'Tag', 'panel1', 'Units', 'pixels', 'Visible', 'off', 'BackgroundColor', [0.95, 0.95, 1]); 

handles.butSVD = uicontrol(handles.panelSVD, 'Style', 'pushbutton', 'Tag', 'butSVD', 'String', 'Calculate', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''butSVD_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 

handles.butTar = uicontrol(handles.panelTar, 'Style', 'pushbutton', 'Tag', 'butTar', 'String', 'Load Target', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''butTar_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 
handles.edTarURange = uicontrol(handles.panelTar, 'Style', 'edit', 'Tag', 'edTarURange', 'String', '1:2', 'Units', 'pixels', 'Callback', 'svdplugin(''edTarURange_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 

handles.popToShow = uicontrol(handles.panelSVD, 'Style', 'popupmenu', 'Tag', 'popToShow', 'String', ...
    {'U', ...
     'S', ...
     'V', ...
     'U*S', ...
     'S*V''',...
     'U(range)*S*V'''}, ...
     'UserData', {'U', 'S', 'V', 'U', 'V', 'U'},...
     'Value', 1, 'Units', 'pixels', 'Callback', 'svdplugin(''popToShow_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 
handles.edRange = uicontrol(handles.panelSVD, 'Style', 'edit', 'Tag', 'edRange', 'String', 'all', 'Units', 'pixels', 'Callback', 'svdplugin(''edRange_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 
handles.txtShowS = uicontrol(handles.panelSVD, 'Style', 'text', 'Tag', 'butShowS', 'String', 'S/max(S) table', 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''butShowS_Callback'',gcbo,[],guidata(gcbo))'); 
handles.lsShowS = uicontrol(handles.panelSVD, 'Style', 'listbox', 'Tag', 'lsShowS', 'String', {' .. press "Calculate" ..'}, ...
    'Value', 1, 'Units', 'pixels', 'Callback', '', ...
    'BackgroundColor', [1, 1, 1]*0.8, 'FontUnits', 'points', 'FontName', FName, 'FontSize', FSize, 'HorizontalAlignment', 'left');
 
handles.popModel = uicontrol(handles.panelGlo, 'Style', 'popupmenu', 'Tag', 'popModel', 'String', ...
    {'Kinetic', ...
     'E.Chem'}, ...
     'UserData', {'ode45', 'Echem'},...
     'Value', 1, 'Units', 'pixels', 'Callback', 'svdplugin(''popToShow_Callback'',[],[], guidata(gcf))', ...
    'BackgroundColor', [1, 1, 1]); 
handles.edNStates = uicontrol(handles.panelGlo, 'Style', 'edit', 'Tag', 'edNStates', 'String', '3', ...
    'Units', 'pixels', 'Callback', 'svdplugin(''edNStates'',[],[], guidata(gcf))', ...
    'BackgroundColor', [1, 1, 1], 'TooltipString', 'Number of States'); 


handles.butAddEq = uicontrol(handles.panelGlo, 'Style', 'pushbutton', 'Tag', 'butRunMod', 'String', 'Add', 'Value', 0, 'Units', 'pixels', ...
    'Callback', 'svdplugin(''butAddEq_Callback'',gcbo,[],guidata(gcbo))'); 

handles.edModel = uicontrol(handles.panelGlo, 'Style', 'edit', 'Tag', 'edStep', 'String', ... 
    {'y0 = [1, 0, 0]*k0', 'dy(1)=-k12*y(1)', 'dy(2)=+k12*y(1)-k23*y(2)', 'dy(3)=k23*y(2)'}, ...
     'Min', 0, 'Max', 4, 'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''edStep_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1], 'FontUnits', 'points', 'FontName', FName, 'FontSize', FSize, 'HorizontalAlignment', 'left'); 

handles.lsParams = uicontrol(handles.panelGlo, 'Style', 'listbox', 'Tag', 'lsParams', 'String', {' .. please update ..'}, ...
    'Value', 1, 'Units', 'pixels', 'Callback', 'svdplugin(''lsParams_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]);
handles.butRunMod = uicontrol(handles.panelGlo, 'Style', 'pushbutton', 'Tag', 'butRunMod', 'String', 'VVV Update VVV', 'Value', 0, 'Units', 'pixels', ...
    'Callback', 'svdplugin(''butRunMod_Callback'',gcbo,[],guidata(gcbo))'); 

handles.edStep = uicontrol(handles.panelGlo, 'Style', 'edit', 'Tag', 'edStep', 'String', '0.1', ...
     'Value', 0, 'Units', 'pixels', 'Callback', 'svdplugin(''edStep_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 

handles.ePars = spinbox(handles.panelGlo, 'Tag', 'ePars', 'Value', 0.5, 'Step', 0.1, ...
        'Units', 'pixels', 'Callback', 'svdplugin(''ePars_Callback'',[],[], guidata(gcf))');
handles.popGlobToShow = uicontrol(handles.panelGlo, 'Style', 'popupmenu', 'Tag', 'popToShow', 'String', ...
    {'Show: Spectra', ...
     'Show: Kinetics and Targets', ...
     'Show: Ut,Ue', ...
     'Show: Vt,Ve'}, ...
     'UserData', {'Sp', 'Ki', 'Us', 'Vs'},...
     'Value', 1, 'Units', 'pixels', 'Callback', 'svdplugin(''popGlobToShow_Callback'',gcbo,[],guidata(gcbo))', ...
    'BackgroundColor', [1, 1, 1]); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function varargout = togMode_Callback(h, ~, handles, varargin)

hButs = [handles.togSVD, handles.togGlo, handles.togTar];
hPanels = [handles.panelSVD, handles.panelGlo, handles.panelTar];
idx = find(hButs==h);
idx1 = find(hButs~=h);

set(hButs(idx), 'Value', 1); set(hButs(idx1), 'Value', 0);
set(hPanels(idx), 'Visible', 'on'); set(hPanels(idx1), 'Visible', 'off');

figure1_ResizeFcn(h, [], handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)

FigaPos = get(handles.figure1, 'Position');

if FigaPos(3)<130, FigaPos(3) = 130; set(handles.figure1, 'Position', FigaPos); end
if FigaPos(4)<210, FigaPos(4) = 210; set(handles.figure1, 'Position', FigaPos); end

brd = 5;


butH = 30;
butW = 60;
% Lets go from top to bottom

hhh  = FigaPos(4)-brd - butH;
pos = [brd, hhh, FigaPos(3)-2*brd, butH];
set(handles.butLoad, 'Position', pos);

hhh  = hhh-brd - butH;
bL = floor((FigaPos(3)-2*brd)/3);

    pos = [brd, hhh, bL, butH];
set(handles.togSVD, 'Position', pos);
    pos = [brd+bL, hhh, bL, butH];
set(handles.togGlo, 'Position', pos);
    pos = [brd+2*bL, hhh, bL, butH];
set(handles.togTar, 'Position', pos);

panelH = hhh - brd;
panPos = [brd, brd, FigaPos(3)-2*brd, panelH];

if get(handles.togSVD, 'Value');
    set(handles.panelSVD, 'Position', panPos);
        hhh  = hhh-4*brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.butSVD, 'Position', pos);
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.popToShow, 'Position', pos);
    
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.edRange, 'Position', pos);
    
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.txtShowS, 'Position', pos);
        
        listH = hhh-2*brd;
        hhh  = hhh-brd-listH;
        pos = [brd, hhh, panPos(3)-2*brd, listH];
    set(handles.lsShowS, 'Position', pos);
    
elseif get(handles.togGlo, 'Value');
    set(handles.panelGlo, 'Position', panPos);    
%     [           ][      ]
%     
%     [                   ]
%     [                   ]
%     [        TB         ]
%     [                   ]
%     [                   ]
%     
%     [    update         ]
% 
%     [                   ]
%     [                   ]
%     [        res        ]
%     [                   ]
%     [                   ]
    
        hhh  = hhh-4*brd - butH;
        pos = [brd, hhh, floor((panPos(3)-2*brd)*0.5), butH];
    set(handles.popModel, 'Position', pos);
        pos = [brd+pos(3), hhh, floor((panPos(3)-2*brd)/4), butH];
    set(handles.edNStates, 'Position', pos);
        pos = [pos(3)+pos(1), hhh, floor((panPos(3)-2*brd)/4), butH];
    set(handles.butAddEq, 'Position', pos);
    
     listH = floor((panelH-4*butH-9*brd)/2);
        hhh = hhh-listH-brd;
        pos = [brd, hhh, panPos(3)-2*brd, listH];
     set(handles.edModel, 'Position', pos);

        hhh  = hhh-brd - 1*butH;
        pos = [brd, hhh, (panPos(3)-2*brd), butH];
    set(handles.butRunMod, 'Position', pos);     
        
        hhh = hhh-listH-brd;
        pos = [brd, hhh, panPos(3)-2*brd, listH];
    set(handles.lsParams, 'Position', pos);    
    
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, floor((panPos(3)-2*brd)/3), butH];
    set(handles.edStep, 'Position', pos);
        pos(1) = pos(1)+pos(3);
        pos(3) = pos(3)*2;
    set(handles.ePars, 'Position', pos);
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.popGlobToShow, 'Position', pos);

    
else
    set(handles.panelTar, 'Position', panPos);    
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.edTarURange, 'Position', pos);
        hhh  = hhh-4*brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
    set(handles.butTar, 'Position', pos);
        hhh  = hhh-brd - butH;
        pos = [brd, hhh, panPos(3)-2*brd, butH];
        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = butLoad_Callback(h, ~, handles)
hand = guidata(handles.hh);
handles.src = hand.src;
[ppath,name,ext] = fileparts(safeget(handles.src, 'filename', 'noname.xxx'));
% n = size(hand.src.y, 2);
% str{1} = 'All';
% for ii=1:n, str{ii+1}=num2str(ii); end
% % set(handles.pmSlice, 'String', str, 'Value', 1);
% while length(handles.set) < n, handles.set(end+1) = handles.set(1); end    
set(h, 'BackgroundColor', [0.64 0.8 0.64]);
set(h, 'String', ['Data: ', name, ext])
a = size(hand.src.y);
handles.rangeU = [1:a(1)];
handles.rangeV = [1:a(2)];
handles.rangeS = [1:min(a)];
guidata(handles.figure1, handles);
Plott_Data(handles);

% pmSlice_Callback(h, eventdata, handles);
% Plott_SVD(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function getRange(handles)

num = get(handles.popToShow, 'Value')
str = get(handles.popToShow, 'UserData');

strG = str{num};

dim = getfield(handles, ['range', strG]);
errmsg = [];
if length(dim)>1
    if ~sum(diff(dim)-dim(2)+dim(1)) 
        if (dim(2)-dim(1) == 1)
            [stringPrint, errmsg] = sprintf('%d:%d', dim(1),dim(end));
        else
            [stringPrint, errmsg] = sprintf('%d:%d:%d', dim(1), dim(2)-dim(1), dim(end));
        end
    else
        [stringPrint, errmsg] = sprintf(' %d,', dim);
        stringPrint = stringPrint(2:(end-1));
    end
else
    if ~isreal(dim)
        stringPrint = 'all';
    else
        [stringPrint, errmsg] = sprintf(' %d,', dim);
        stringPrint = stringPrint(2:(end-1));
    end
end

if ~isempty(errmsg),
    stringPrint = 'all';
end

set(handles.edRange, 'String', stringPrint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = edTarURange_Callback(h, ~, handles)
str = get(handles.edTarURange, 'String');
ran = str2num(str);


Data = getfield(handles.SVD, 'U');
if min(ran)<1
    idx = find(ran>=1);
    ran = ran(idx);
end
if max(ran)>size(Data, 1)
    idx = find(ran<=size(Data, 1));
    ran = ran(idx);
end

% guidata(handles.figure1, handles);
butTar_Callback(h, 0, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = edRange_Callback(h, ~, handles)

str = get(handles.edRange, 'String');
ran = str2num(str);

num = get(handles.popToShow, 'Value');
str = get(handles.popToShow, 'UserData');

strG = str{num};

Data = getfield(handles.SVD, strG);
if min(ran)<1
    idx = find(ran>=1);
    ran = ran(idx);
end
if max(ran)>size(Data, 1)
    idx = find(ran<=size(Data, 1));
    ran = ran(idx);
end
handles = setfield(handles, ['range', strG], ran);
guidata(handles.figure1, handles);
Plott_SVD(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function butTar_Callback(h, ~, handles)

hand = guidata(handles.hh);
handles.target = hand.src;
[ppath,name,ext] = fileparts(safeget(hand.src, 'filename', 'noname.xxx'));
% n = size(hand.src.y, 2);
% str{1} = 'All';
% for ii=1:n, str{ii+1}=num2str(ii); end
% % set(handles.pmSlice, 'String', str, 'Value', 1);
% while length(handles.set) < n, handles.set(end+1) = handles.set(1); end    
set(handles.butTar, 'BackgroundColor', [0.64 0.8 0.64]);
set(handles.butTar, 'String', ['Data: ', name, ext])
a = size(handles.target.y);
handles.TarNsets = a(2);
handles.TarPoints = a(1);

% handles.rangeV = [1:a(2)];
% handles.rangeS = [1:min(a)];

guidata(handles.figure1, handles);

if ~isempty(handles.SVD.U)
    range = str2num(get(handles.edTarURange, 'String'));
    U = handles.SVD.U(:, range);
    P = U*U';
    
    for ii = 1:handles.TarNsets
        if sum(handles.src.ax.x-handles.src.ax.x)~=0
            TY = spline(handles.target.ax.x, handles.target.y(:, ii), handles.src.ax.x);
        else
            TY = handles.target.y(:, ii);
        end
        Res(:, ii) = P*TY;
        
    end
    
    
    handles.tarRes = handles.target;
    handles.tarRes.y = Res;
    handles.tarRes.ax = handles.src.ax;
    
    KVhand = guidata(handles.hh);

    KVhand.out = {};
    KVhand.out{1} = handles.target;
    KVhand.out{1}.title = 'Target';
    KVhand.out{1}.ax.Color = [1, 0, 0];
    
    KVhand.out{2} = handles.tarRes;
    KVhand.out{2}.title = 'Result';
    KVhand.out{2}.ax.Color = [0, 0, 1];
    
    guidata(handles.hh, KVhand);
    [ppath,name,ext] = fileparts(get(handles.hh, 'FileName'));
    eval([name '(''SetDataSource'', 2, KVhand)']);
    
else
    warndlg('No SVD calculation results found');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plott_Data(handles)
h = guidata(handles.figure1);
KVhand = guidata(h.hh);

KVhand.out = {};
KVhand.out{1} = handles.src;

guidata(h.hh, KVhand);
[ppath,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, KVhand)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plott_SVD(handles)

ha = guidata(handles.figure1);
KVhand = guidata(ha.hh);

KVhand.out = {};

if isempty(ha.src),
    warning('No reference data loaded');
    set(ha.butLoad, 'BackgroundColor', [1 0.64 0.64]);
    return;
end

KVhand.out{1}.ax = ha.src.ax;
switch get(ha.popToShow, 'Value')
    case 1 % U
        KVhand.out{1}.y = ha.SVD.U(:, ha.rangeU);
        KVhand.out{1}.ax.x = ha.src.ax.x;
        KVhand.out{1}.ax.xlabel = ha.src.ax.xlabel;
        KVhand.out{1}.ax.y = ha.rangeU;
        KVhand.out{1}.ax.ylabel = 'pnts';
        col = getLineColors(length(ha.rangeU));
        KVhand.out{1}.title = 'U';
        KVhand.out{1}.ax.Color = col;
    case 2 % S
        dia = diag(ha.SVD.S);
        KVhand.out{1}.y = dia(ha.rangeS, 1);  
        KVhand.out{1}.ax.x = [ha.rangeS];
        KVhand.out{1}.ax.xlabel = 'pnts';
        KVhand.out{1}.ax.y = 1;
        KVhand.out{1}.ax.ylabel = 'pnts';
        KVhand.out{1}.ax.Marker = '*';
        KVhand.out{1}.title = 'S';
        
    case 3 % V
        KVhand.out{1}.y = ha.SVD.V(:, ha.rangeV);
        KVhand.out{1}.ax.x = ha.src.ax.y;
        KVhand.out{1}.ax.xlabel = ha.src.ax.ylabel;
        KVhand.out{1}.ax.y = [ha.rangeV];
        KVhand.out{1}.ax.ylabel = 'pnts';
        KVhand.out{1}.title = 'V';
    case 4 %U*S
        KVhand.out{1}.y = ha.SVD.U(:, ha.rangeU)*ha.SVD.S(ha.rangeU, : );
        KVhand.out{1}.ax.x = ha.src.ax.x;
        KVhand.out{1}.ax.xlabel = ha.src.ax.xlabel;
        KVhand.out{1}.ax.y = ha.rangeU;
        KVhand.out{1}.ax.ylabel = 'pnts';
        col = getLineColors(length(ha.rangeU));
        KVhand.out{1}.ax.Color = col;
        KVhand.out{1}.title = 'U*S';
    case 5 % S*V'
        KVhand.out{1}.y = ha.SVD.S(:, ha.rangeV)*ha.SVD.V(:, ha.rangeV).';
        KVhand.out{1}.ax.x = ha.src.ax.x;
        KVhand.out{1}.ax.xlabel = ha.src.ax.xlabel;
        KVhand.out{1}.ax.y = [1:size(ha.SVD.U, 1)];
        KVhand.out{1}.ax.ylabel = 'pnts';
        col = getLineColors(length(ha.rangeV));
        KVhand.out{1}.ax.Color = col;    
        KVhand.out{1}.title = 'S*V''';
    case 6 % U*S*V'
        KVhand.out{1}.ax = ha.src.ax;
        KVhand.out{1}.y = ha.src.y;
        KVhand.out{1}.title = 'SRC';
        KVhand.out{2}.ax = ha.src.ax;
        KVhand.out{2}.y = ha.SVD.U(:, ha.rangeU)*ha.SVD.S(ha.rangeU, :)*ha.SVD.V.';
        KVhand.out{2}.ax.x = ha.src.ax.x;
        KVhand.out{2}.ax.xlabel = ha.src.ax.xlabel;
        KVhand.out{2}.ax.y = [1:size(ha.SVD.U, 1)];
        KVhand.out{2}.ax.ylabel = 'pnts';
        KVhand.out{2}.ax.Color = 'r';       
        KVhand.out{2}.title = 'U*S*V''';
end

guidata(ha.hh, KVhand);
[ppath,name,ext] = fileparts(get(ha.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, KVhand)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plott_GlobAn(handles)

ha = guidata(handles.figure1);
KVhand = guidata(ha.hh);

KVhand.out = {};

if isempty(ha.src),
    warning('No reference data loaded');
    set(ha.butLoad, 'BackgroundColor', [1 0.64 0.64]);
    return;
end

Mode = get(handles.popGlobToShow, 'Value');

switch Mode
    case 1
        KVhand.out{1}.title = 'Spec.Comp.';
        KVhand.out{1}.ax.x = handles.Glob.ResX;
        KVhand.out{1}.ax.y = [1:size(handles.Glob.ResSpec, 2)];
        KVhand.out{1}.y = handles.Glob.ResSpec;
        KVhand.out{1}.ax.xlabel = ha.src.ax.xlabel;
        KVhand.out{1}.ax.ylabel = 'Spec #';
        KVhand.out{1}.ax.Color = 'b';
    case 2
        KVhand.out{1}.title = 'Kinetics';
        KVhand.out{1}.ax.x = handles.Glob.ResTime;
        KVhand.out{1}.ax.y = [1:size(handles.Glob.ResKin, 2)];
        KVhand.out{1}.y = handles.Glob.ResKin;
        KVhand.out{1}.ax.xlabel = ha.src.ax.ylabel;
        KVhand.out{1}.ax.xlabel = 'Spec #';
        KVhand.out{1}.ax.Color = 'b';

        KVhand.out{2}.title = 'Kinetics';
        KVhand.out{2}.ax.x = handles.Glob.ResTime;
        KVhand.out{2}.ax.y = [1:size(handles.Glob.ResKinTarget, 2)];
        KVhand.out{2}.y = handles.Glob.ResKinTarget;
        KVhand.out{2}.ax.xlabel = ha.src.ax.ylabel;
        KVhand.out{2}.ax.xlabel = 'Spec #';
        KVhand.out{2}.ax.Color = 'r';
        
end
        
guidata(ha.hh, KVhand);
[ppath,name,ext] = fileparts(get(ha.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, KVhand)']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = butSVD_Callback(h, ~, handles)

if min(size(handles.src.y))==1,
    set(h, 'BackgroundColor', [1, 0.8, 0.8]);
    hx = errordlg('Loaded dataset has only one dimention. \n SVD is useless.', 'Wrond data !!!');
    return;
end

[U, S, V] = svd(handles.src.y);
handles.SVD.U = U;
handles.SVD.S = S;
handles.SVD.V = V;
guidata(handles.figure1, handles);
set(h, 'BackgroundColor', [0.8, 1, 0.8]);
Update_ListS(handles);
Plott_SVD(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Update_ListS(handles)
Sval = diag(handles.SVD.S);
idx = find(Sval> max(Sval)*1e-2);
idx = max(idx);
str = cell(length(Sval), 1);
for ii = 1:length(Sval)
    str{ii} = sprintf('%4d: %f', ii, Sval(ii)/max(Sval));
end
set(handles.lsShowS, 'String', str, 'Value', idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = popToShow_Callback(h, ~, handles)
getRange(handles);
Plott_SVD(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function col = getLineColors(N)

a = {[0         0    1.0000], ...
    [     0    0.5000         0], ...
    [1.0000         0         0], ...
    [     0    0.7500    0.7500], ...
    [0.7500         0    0.7500], ...
    [0.7500    0.7500         0], ...
    [0.2500    0.2500    0.2500], ...
    [0, 0, 0]};

b = repmat(a ,1, floor(N/8)+1);

col = b(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = butRunMod_Callback(h, ~, handles)
handles = guidata(handles.figure1);
Mode = get(handles.popModel, 'Value');

if Mode  ==1
    prepKineticAn(handles);
elseif Mode ==2
    prepEChemAn(handles);
end
handles = guidata(handles.figure1);
Update_ListPars(handles);
set(h, 'BackgroundColor', [0.8, 1, 0.8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepEChemAn(handles)
handles = guidata(handles.figure1);
handles.Glob.Mode = 2;
string = get(handles.edModel, 'String');
[variables, stats]=ParseString(string);

handles.Glob.variables = variables;
handles.Glob.stats = stats;

if isfield(handles.Glob, 'K')
    if length(handles.Glob.K)>stats.countR % cut
        handles.Glob.K = handles.Glob.K(1:stats.countR);
    elseif length(handles.Glob.K)<stats.countR % extend
        a = handles.Glob.K;
        len = length(handles.Glob.K);
        handles.Glob.K = ones(1, stats.countR);
        handles.Glob.K(1:len) = a;
    end
else % create
    handles.Glob.K = ones(1, stats.countR);
end
guidata(handles.figure1, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepKineticAn(handles)
handles = guidata(handles.figure1);
handles.Glob.Mode = 1;
string = get(handles.edModel, 'String');
[variables, stats]=ParseString(string);
if stats.countDC <stats.countC
    set(h, 'BackgroundColor', [1, 0.8, 0.8]);
    hx = errordlg({'Number of "dy" is less than "y". ', 'check the model.'}, 'N(dy)<N(y) !!!');
    return;
end

if stats.countDC <stats.countC
    set(h, 'BackgroundColor', [1, 0.8, 0.8]);
    hx = errordlg({'Number of "dy" is less than "y". ', 'check the model.'}, 'N(dy)<N(y) !!!');
    return;
end

% test ....
y0 = []; 
dy = ones(stats.countDC, 1);
y = ones(stats.countDC, 1);

for ii = 1:stats.countR
    eval([variables{ii}, ' = 1;']);
end
for ii = 1:length(string)
    try
        eval([string{ii}, ';']);
    catch
        set(h, 'BackgroundColor', [1, 0.8, 0.8]);
        hx = errordlg({string{ii}, lasterr}, 'matlab expression error');
        return;
    end
end

if isempty(y0)
    set(h, 'BackgroundColor', [1, 0.8, 0.8]);
    hx = errordlg({'Could not find y0. ', 'please specify y0 = [1, 0 ... etc ] ', 'in the beginning of the script.'}, 'no initial concentration !!!');
    return;
end
if size(y0)~=stats.countDC
    set(h, 'BackgroundColor', [1, 0.8, 0.8]);
    hx = errordlg({'Number of elements in "y0" is less than "dy". ',' check the model.'}, 'N(y0)<N(dy) !!!');
    return;
end

handles.Glob.variables = variables;
handles.Glob.stats = stats;

if isfield(handles.Glob, 'K')
    if length(handles.Glob.K)>stats.countR % cut
        handles.Glob.K = handles.Glob.K(1:stats.countR);
    elseif length(handles.Glob.K)<stats.countR % extend
        a = handles.Glob.K;
        len = length(handles.Glob.K);
        handles.Glob.K = ones(1, stats.countR);
        handles.Glob.K(1:len) = a;
    end
else % create
    handles.Glob.K = ones(1, stats.countR);
end

script_sort = cell(handles.Glob.stats.countDC, 1);
scrrr = {};
for ii = 1:length(string)
    if strfind(string{ii}, 'y0')
        handles.Glob.y0 = string{ii};
    else
       scrrr{end+1} = string{ii}; % clumzy, but I am lazy to do it better rightnow
%         [a, b] = strtok(string{ii}, '=');
%         num = sscanf(a, 'dy(%d)');
%         if ~isempty(num)
%             script_sort{num} = b(1:end);
%         end
    end
end

% handles.Glob.script = '[';
% for ii = 1:handles.Glob.stats.countDC
%     handles.Glob.script = [handles.Glob.script, script_sort{ii}, ';'];
% end
% handles.Glob.script = [handles.Glob.script(1:(end-1)),  ']'];

[pa, na, ex] = fileparts(which('svdplugin'));
fid = fopen([pa, '\svd_tmpfunc.m'], 'w');
fprintf(fid, 'function dy = svd_tmpfunc(t, y) \n');
    str = 'global'
    for ii = 1:stats.countR
        str = [str, ' ', variables{ii}]; 
    end
    str = [str, '; \n'];
fprintf(fid, str);
fprintf(fid, 'dy = zeros(%d,1); \n', stats.countDC);
for ii = 1:length(scrrr)
    fprintf(fid, '%s ; \n', scrrr{ii});
end
fclose(fid);
% handles.Glob.fHandle = inline(handles.Glob.script, 't', 'y');
guidata(handles.figure1, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variables,  stats, infunction]=ParseString(str)
rates = 'k';
amps = 'a';
concentr = 'y';
difconc = 'dy';

nernst = 'nernst';

toSpace = {'y0', './', '.*', '.^', '+', '-', '=', '/', '*', '^'};
toDelete = {'(', ')'};

strRead = '';
for ii = 1:length(str)
    strRead = [strRead, str{ii}, ' '];
end

for ii = 1:length(toSpace)
    idx = strfind(strRead, toSpace{ii});
    strRead(idx) = ' ';
end

ss = ones(length(strRead), 1);
for ii = 1:length(toDelete)
    idx = strfind(strRead, toDelete{ii});
    ss(idx) = 0;
end
idx = find(ss);
strRead = strRead(idx);

countR = 0;
countA = 0;
countC = 0;
countDC = 0;
countN = 0;
countNR = [];
Srates = '';
Samps = '';
Sconc = '';
SDconc = '';
Snernst = '';
while ~isempty(strRead)
    [a, strRead] = strtok(strRead);
    if strfind(a, rates)
        if isempty(strfind(Srates, [a, ' ']))
            countR = countR+1;
            allRates{countR} = a;
            Srates = [Srates, a, ' '];
        end
    elseif strfind(a, amps)
        if isempty(strfind(Samps, [a, ' ']))
            countA = countA+1;
            allAmps{countA} = a;
            Samps = [Samps, a, ' '];
        end
    elseif strfind(a, difconc)   
        if isempty(strfind(SDconc, [a, ' ']))
            countDC = countDC+1;
            allDConc{countDC} = a;
            SDconc = [SDconc, a, ' '];
        end
    elseif strfind(a, concentr)
        if isempty(strfind(Sconc, [a, ' ']))
            countC = countC+1;
            allConc{countC} = a;
            Sconc = [Sconc, a, ' '];
        end
    elseif strfind(a, nernst)
        if isempty(strfind(Snernst, [a, ' ']))
            countN = countN+1;
            allNernst{countN} = a;
            Snernst = [Snernst, a, ' '];
            Nnernst = str2num(a(7:end));

            if Nnernst>0
                countNR(countN) = Nnernst;
                for ii = 1:Nnernst
                    countR = countR+1;
                    allRates{countR} = sprintf('E%d_%d', countN, ii);
                end
                countR = countR+1;
                allRates{countR} = sprintf('Ne_%d', countN);
            end
        end
    end
    
end

variables = allRates;
% conc = allConc;
% Dconc = allDConc;
stats.countR = countR;
stats.countA = countA;
stats.countC = countC;
stats.countDC = countDC;
stats.countN = countN;
stats.countNR = countNR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Update_ListPars(handles)
if ~isfield(handles, 'Glob'), return; end
if ~isfield(handles.Glob, 'stats'), return; end

str = cell(handles.Glob.stats.countR, 1);
for ii = 1:handles.Glob.stats.countR
    str{ii} = sprintf('%s = %f', handles.Glob.variables{ii}, handles.Glob.K(ii));
end

set(handles.lsParams, 'String', str, 'Value', 1);
lsParams_Callback(handles.figure1, [], handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = lsParams_Callback(h, ~, handles)

Res = get(handles.lsParams, 'String');
val = get(handles.lsParams, 'Value');
if size(Res)<2, 
   val = 1; 
end

num= sscanf(Res{val}, '%*s = %f');

set(handles.ePars, 'String', num2str(num));
set(handles.ePars, 'Value', num);
set(handles.ePars, 'UserData', 'value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = edStep_Callback(h, ~, handles)
ss = get(h, 'String');
num = str2num(ss);
set(handles.ePars, 'Step', num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ePars_Callback(h, ~, handles)
Res = get(handles.lsParams, 'String');
val = get(handles.lsParams, 'Value');
newval = get(handles.ePars, 'Value');

name = strtok(Res{val});

Res{val} = sprintf('%s = %f', name, newval);

set(handles.lsParams, 'String', Res);

for ii = 1:length(Res)
    num = sscanf(Res{ii}, '%*s = %f');
    handles.Glob.K(ii) = num;
end

guidata(handles.figure1, handles);
if handles.Glob.Mode == 1 % kinetics
    [T, Y] = calc_eq(handles.Glob.y0, handles.Glob.variables, handles.Glob.K, handles.src.ax.y);
elseif handles.Glob.Mode == 2 % echem
%     Y = [];
%     
%     for ii = 1:handles.Glob.nEchem
        [T, Y] = calc_EChem(handles.Glob.stats.countNR, handles.Glob.variables, handles.Glob.K, handles.src.ax.y);
%         [T, tY] = calc_EChem(handles.Glob.E0s{ii}, handles.Glob.Nel(ii), handles.src.ax.y);
%         Y = [Y, tY];
%     end
end
[Spec, yR, tarY] = globalAn(handles.SVD.U, handles.SVD.S, handles.SVD.V, T, Y);
handles.Glob.ResSpec = Spec;
handles.Glob.ResKin = Y;
handles.Glob.ResTime = T;
handles.Glob.ResSU = yR;
handles.Glob.ResKinTarget = tarY;
handles.Glob.ResX =  handles.src.ax.x;

guidata(handles.figure1, handles);

Plott_GlobAn(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, Y] = calc_eq( y0str, vars, vals, t)
str = 'global ';
for ii = 1:length(vars)
    str = [str, ' ', vars{ii}];
end
str = [str, ';'];
eval(str);
for ii = 1:length(vars)
    eval([vars{ii},' = ',sprintf('%f', vals(ii))]);
end

ODEoptions  = [];
% ODEoptions = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
eval(y0str);
[T,Y] = ode45('svd_tmpfunc',t, y0, ODEoptions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vs, Y] = calc_EChem(nVars, vars, vals, Vs)
count = 0;
Y = [];
Con = 50;
for kk = 1:numel(nVars) %%%%%%%% number of independent EChem processes
    
    E0s = sort(vals( [1:nVars(kk)]+count ) , 'ascend');
    Ne = vals(nVars(kk)+1);
    count = nVars(kk)+1;
    %%%%%%%%%%% can be optimised further for speed!
    nP = numel(E0s)+1;
    Mat = zeros(nP);
    Mat(1, :) = 1;
    for ii = 2:nP
        Mat(ii, ii-1)= 1;
    end
    V0 = zeros(nP, 1); V0(1, 1)=1;
    
    for vv = 1:numel(Vs)
        for ii = 2:nP
            Mat(ii, ii)= -exp( (E0s(ii-1)-Vs(vv))/Con/Ne);
        end
        dM = det(Mat);
        for ii = 1:nP
            tM = Mat;
            tM(:, ii) = V0;
            tY(vv, ii) = det(tM)/dM;
        end
    end
    Y = [Y, tY];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Spec, yR, tarY] = globalAn(U, S, V, T, Y)
nP = size(T);
nSpec = size(Y, 2);

RR = S(1:nSpec, 1:nSpec)*V(:, 1:nSpec).'/Y';
Spec = U(:, 1:nSpec)*RR;
yR = (RR*Y').';

%%%%%%%% target testing kinetics
tV = V(:, 1:nSpec);
P = tV*tV';
for ii = 1:nSpec
    tarY(:, ii) = P*Y(:, ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popGlobToShow_Callback(h, ~, handles)
guidata(handles.figure1, handles);
Plott_GlobAn(handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function butAddEq_Callback(h, ~, handles)
Mode = get(handles.popModel, 'Value');
str = get(handles.edModel, 'String');
[variables,  stats] = ParseString(str);

nPrev = stats.countR;
nNew = str2num(get(handles.edNStates, 'String'));

if Mode == 1

nPrev = stats.countR;
nNew = str2num(get(handles.edNStates, 'String'));

for ii = 1:nNew
    tS = sprintf('dy(%d)=-y(%d)*(', nPrev+ii, nPrev+ii);
    for jj = 1:nNew
        if ii==jj, continue; end
        tS = sprintf('%s+k%d%d', tS, nPrev+ii, nPrev+jj);
    end
    tS = [tS, ')'];
    for jj = 1:nNew
        if ii==jj, continue; end
        tS = sprintf('%s+k%d%d*y(%d)', tS, nPrev+jj, nPrev+ii, nPrev+jj);
    end
    str{end+1} = tS;
end


end

if Mode == 2
    if stats.countN == 0
        str = {};
    end
    
    str{end+1} = sprintf('Y%d=a%d*nernst(%d)', stats.countN+1, stats.countN+1, nNew-1);

end 
% disp(str)
    
set(handles.edModel, 'String', str);







