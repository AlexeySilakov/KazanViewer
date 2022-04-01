function varargout = CompBox(varargin)
%@plugin = 'CompBox';
%@helpfile = 'none';
% COMPBOX Application M-file for CompBox.fig
%    FIG = COMPBOX launch CompBox GUI.
%    COMPBOX('callback_name', ...) invoke the named callback.

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.0 02-Jan-2005 11:41:30
% Alexey Silakov & Boris Epel 5-dec-2003, MPI
% boep 21-dec-2003, MPI

if nargin == 0  % LAUNCH GUI
    token = 'COMPBOX_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,n,e] = fileparts(which('kazan'));
    inifilename = [fpath , filesep, 'kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try opc = ini.KazanViewer.MultyKazan; end
    if (opc == 1| opc == 2) & ~isempty(oldfig), set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; end

	fig = figure; %openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    set(fig, 'Units', 'pixels', 'Name', 'CompBOX', 'NumberTitle', 'off');
    set(fig, 'Position', [319 600 200 400], 'Tag', 'figure1', 'Visible', 'off', 'MenuBar', 'none'); %%% hide until everything is done
    set(fig, 'ResizeFcn', 'CompBox(''figure1_ResizeFcn'',[],[],guidata(gcf))');
    
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.figure1 = fig;
    %%%%%%%%%%%%%% selector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.popfiles = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'popfiles', ...
        'Units','pixels', 'String', '--nothing loaded--', 'Callback', 'CompBox(''popfiles_Callback'',gcbo,[],guidata(gcbo))');
    %%%%%%%%%%% color/title
    handles.text6 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text6', ...
        'Units','pixels', 'String', '');
    handles.CollCur = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'CollCur', ...
        'Units','pixels', 'String', 'Color', 'Callback', 'CompBox(''CollCur_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbTitle = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'pbTitle', ...
        'Units','pixels', 'String', 'Title', 'Callback', 'CompBox(''pbTitle_Callback'',gcbo,[],guidata(gcbo))');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shifts %%%%%%%%%%%%%%%%%%%%%%
    handles.frame2 = uicontrol(handles.figure1,'Style', 'frame', 'Tag', 'frame2', ...
        'Units','pixels', 'String', '', 'BackgroundColor', [0.5, 0.5, 0.5]);
    
    handles.text1 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text1', ...
        'Units','pixels', 'String', 'Shift X', 'BackgroundColor', [0.5, 0.5, 0.5], 'FontSize', 9);
    handles.text4 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text4', ...
        'Units','pixels', 'String', 'Shift Y', 'BackgroundColor', [0.5, 0.5, 0.5], 'FontSize', 9);
    handles.SXdim = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'SXdim', ...
        'Units','pixels', 'Callback', 'CompBox(''SXdim_Callback'',gcbo,[],guidata(gcbo))');
    handles.SYdim = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'SYdim', ...
        'Units','pixels', 'Callback', 'CompBox(''SYdim_Callback'',gcbo,[],guidata(gcbo))');
    
    set([handles.SXdim, handles.SYdim], 'String', {'1e-7','1e-6','1e-5', '1e-4', '0.001','0.01','0.1','1','10','100','1000','1e4','1e5','1e6','1e7','1e8'});
    set([handles.SXdim, handles.SYdim], 'Value', 7);
    
    handles.ShiftX = spinbox(handles.figure1, 'String', '0', 'ButtonStyle', 'right', 'Direction', 'horizontal', ...
        'EditBackground', [1, 1, 1], 'Step', 0.1, 'Callback', 'CompBox(''ShiftX_Callback'',[],[],guidata(gcf))');
    handles.ShiftY = spinbox(handles.figure1, 'String', '0', 'ButtonStyle', 'right', 'Direction', 'vertical', ...
        'EditBackground', [1, 1, 1], 'Step', 0.1, 'Callback', 'CompBox(''ShiftY_Callback'',[],[],guidata(gcf))');
    
    handles.ClearX = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'ClearX', ...
        'Units','pixels', 'String', 'Set Zero', 'Callback', 'CompBox(''ClearX_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbNorm = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'pbNorm', ...
        'Units','pixels', 'String', 'Norm', 'Callback', 'CompBox(''pbNorm_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbToFirst = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'pbToFirst', ...
        'Units','pixels', 'String', '->1st', 'Callback', 'CompBox(''pbToFirst_Callback'',gcbo,[],guidata(gcbo))');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Zooms %%%%%%%%%%%%%%%%%%%%%%
    handles.frame1 = uicontrol(handles.figure1,'Style', 'frame', 'Tag', 'frame1', ...
        'Units','pixels', 'String', 'ZOOM', 'BackgroundColor', [0.5, 0.5, 0.5]);
    
    handles.text2 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text2', ...
        'Units','pixels', 'String', 'Zoom X', 'BackgroundColor', [0.5, 0.5, 0.5], 'FontSize', 9);
    handles.text3 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text3', ...
        'Units','pixels', 'String', 'Zoom Y', 'BackgroundColor', [0.5, 0.5, 0.5], 'FontSize', 9);
    
    handles.ZXdim = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'ZXdim', ...
        'Units','pixels', 'Callback', 'CompBox(''ZXdim_Callback'',gcbo,[],guidata(gcbo))');
    handles.ZYdim = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'ZYdim', ...
        'Units','pixels', 'Callback', 'CompBox(''ZYdim_Callback'',gcbo,[],guidata(gcbo))');
    
    set([handles.ZXdim, handles.ZYdim], 'String', {'1e-7','1e-6','1e-5', '1e-4', '0.001','0.01','0.1','1','10','100','1000','1e4','1e5','1e6','1e7','1e8'});
    set([handles.ZXdim, handles.ZYdim], 'Value', 7);
    
    handles.ZoomX = spinbox(handles.figure1, 'String', '1', 'ButtonStyle', 'right', 'Direction', 'horizontal', ...
        'Step', 0.1, 'Callback', 'CompBox(''ZoomX_Callback'',[],[],guidata(gcf))');
    handles.ZoomY = spinbox(handles.figure1, 'String', '1', 'ButtonStyle', 'right', 'Direction', 'vertical', ...
        'Step', 0.1, 'Callback', 'CompBox(''ZoomY_Callback'',[],[],guidata(gcf))');
    
    handles.ClearY = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'ClearY', ...
        'Units','pixels', 'String', 'Set Zero', 'Callback', 'CompBox(''ClearY_Callback'',gcbo,[],guidata(gcbo))');
    handles.text9 = uicontrol(handles.figure1,'Style', 'text', 'Tag', 'text9', ...
        'Units','pixels', 'String', 'Y scale');
    handles.popUnits = uicontrol(handles.figure1,'Style', 'popupmenu', 'Tag', 'popUnits', ...
        'Units','pixels', 'String', {'Linear', 'dB'}, 'Callback', 'CompBox(''popUnits_Callback'',gcbo,[],guidata(gcbo))');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.eSecondDim = uicontrol(handles.figure1,'Style', 'edit', 'Tag', 'eSecondDim', ...
        'Units','pixels', 'String', '', 'Callback', 'CompBox(''editctl_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.AddButt = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'AddButt', ...
        'Units','pixels', 'String', 'ADD NEW DATA', 'Callback', 'CompBox(''AddButt_Callback'',gcbo,[],guidata(gcbo))');
    handles.pushbutton8 = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'pushbutton8', ...
        'Units','pixels', 'String', 'Inverse Y', 'Callback', 'CompBox(''pbInvY_Callback'',gcbo,[],guidata(gcbo))');
    handles.delCurr = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'delCurr', ...
        'Units','pixels', 'String', 'Delete Current', 'Callback', 'CompBox(''DelCurr_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbChange = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'lDirlist', ...
        'Units','pixels', 'String', 'Replace', 'Callback', 'CompBox(''pbChange_Callback'',gcbo,[],guidata(gcbo))');
    handles.Delete = uicontrol(handles.figure1,'Style', 'pushbutton', 'Tag', 'lDirlist', ...
        'Units','pixels', 'String', 'delete ALL', 'Callback', 'CompBox(''Dell_Callback'',gcbo,[],guidata(gcbo))');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MENU
    handles.mProces = uimenu(handles.figure1, 'Label','Processing',              'Tag', 'mProces');
    handles.mFittFst = uimenu(handles.mProces, 'Label','Fit Current to 1st',        'Tag', 'mFittFst',      'Callback', 'CompBox(''pbToFirst_Callback'',gcbo,[],guidata(gcbo))');
    handles.mFittAllFst = uimenu(handles.mProces, 'Label','Fit All to 1st',         'Tag', 'mFittAllFst',   'Callback', 'CompBox(''mFittAllFst_Callback'',gcbo,[],guidata(gcbo))');
    handles.mNormCurrent = uimenu(handles.mProces, 'Label','Normalise Current',     'Tag', 'mNormCurrent',  'Callback', 'CompBox(''pbNorm_Callback'',gcbo,[],guidata(gcbo))');
    handles.mNormSlices = uimenu(handles.mProces, 'Label','Normalise Current 2D',   'Tag', 'mNormSlices',   'Callback', 'CompBox(''mNormSlices_Callback'',gcbo,[],guidata(gcbo))');
    handles.mNormalize = uimenu(handles.mProces, 'Label','Normalise All Traces',    'Tag', 'mNormalize',    'Callback', 'CompBox(''mNormalize_Callback'',gcbo,[],guidata(gcbo))');
    handles.mPreserveZero = uimenu(handles.mProces, 'Label','mPreserveZero',        'Tag', 'mPreserveZero', 'Callback', 'CompBox(''mPreserveZero_Callback'',gcbo,[],guidata(gcbo))');
    handles.mWaterfall = uimenu(handles.mProces, 'Label','Waterfall',               'Tag', 'mWaterfall',    'Callback', 'CompBox(''mWaterfall_Callback'',gcbo,[],guidata(gcbo))');
    handles.mWaterfall2D = uimenu(handles.mProces, 'Label','Waterfall 2D',          'Tag', 'mWaterfall2D',  'Callback', 'CompBox(''mWaterfall_Callback'',gcbo,[],guidata(gcbo))');
    handles.mFittAddToOut = uimenu(handles.mProces,'Label','Add Current To Output', 'Tag', 'mFittAddToOut', 'Callback', 'CompBox(''mFittAddToOut_Callback'',gcbo,[],guidata(gcbo))');
    handles.m1Dto2D = uimenu(handles.mProces, 'Label','Join Traces',                'Tag', 'm1Dto2D',       'Callback', 'CompBox(''m1Dto2D_Callback'',gcbo,[],guidata(gcbo))');
    handles.m2Dto1D = uimenu(handles.mProces, 'Label','Ungroup Data',               'Tag', 'm2Dto1D',       'Callback', 'CompBox(''m2Dto1D_Callback'',gcbo,[],guidata(gcbo))');
    handles.mSortIndex = uimenu(handles.mProces, 'Label','Sort By Y Value',         'Tag', 'mSortIndex',    'Callback', 'CompBox(''mSortIndex_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.LO = uimenu(handles.figure1, 'Label','Load Options');
    handles.mPrAgr = uimenu(handles.LO, 'Label','Sort By Y Value',         'Tag', 'mPrAgr',    'Callback', 'CompBox(''mPrAgr_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.x = 0;
    handles.y = 0;
    handles.trace = [];
    handles.Color = [0, 0, 1;...
                     1, 0, 0;...
                     0, 0, 1;...
                     0, 0, 0;...
                     1, 0, 1;...
                     0, 1, 1;...
                     1, 1, 0];
    handles.wfshift = 0;
    handles.hh = 0;
    handles.ZoomUnits = 'Norm';
    set(handles.mPrAgr, 'Checked', 'on');
    handles.process = 1;
    handles.sumall = 0;
    handles.change2D = 0;
    
    set(fig, 'Visible', 'on');
    guidata(fig, handles);
    

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
function varargout = FileUpdate(varargin)

% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
t = get(handles.figure1, 'position');
figW = t(3);
figH = t(4);
sp = 5;
butH = ((figH-(sp*18))/12);
if butH<12
    butH = 12;
    figH = 12*butH+sp*18;
    set(handles.figure1, 'Position', [t(1:3), figH]);
end

if (figW-3*sp)<150 
    figW = 150+3*sp; 
    set(handles.figure1, 'Position', [t(1:2), figW, figH]);
end

hh = sp;

butW = floor((figW-3*sp)/2);
set(handles.delCurr, 'Position', [sp, hh, butW, butH]);
set(handles.Delete, 'Position', [2*sp+butW, hh, butW, butH]);
    hh = hh+butH+sp;
set(handles.pushbutton8, 'Position', [sp, hh, butW, butH]);
set(handles.pbChange, 'Position', [2*sp+butW, hh, butW, butH]);
    hh = hh+butH+sp;
set(handles.AddButt, 'Position', [sp, hh, floor(figW-2*sp), butH]);
    hh = hh+butH+2*sp;
    
set(handles.eSecondDim, 'Position', [sp, hh, floor(figW-2*sp), butH]);
    hh = hh+butH+sp;

    %%%%%%%%%%%%%%%%%%%%%% ZOOM panel %%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.frame1, 'Position', [sp/2, hh, floor(figW-sp), 3*butH+4*sp]);
    hh = hh+sp;
    
txW = 50;

butW = floor((figW-4*sp-txW)/2);
set(handles.ClearY, 'Position', [sp, hh, butW, butH]); 
set(handles.text9, 'Position', [2*sp+butW, hh, txW, butH]); 
set(handles.popUnits, 'Position', [3*sp+butW+txW, hh, butW, butH]); 
    hh = hh+butH+sp;
ppW = 70;
set(handles.text3, 'Position', [sp, hh, txW, butH]);
set(handles.ZYdim, 'Position', [2*sp+txW, hh, ppW, butH]);
set(handles.ZoomY, 'Position', [2*sp+txW+ppW, hh, figW-txW-ppW-4*sp, butH]);
    hh = hh+butH+sp;
set(handles.text2, 'Position', [sp, hh, txW, butH]);
set(handles.ZXdim, 'Position', [2*sp+txW, hh, ppW, butH]);
set(handles.ZoomX, 'Position', [2*sp+txW+ppW, hh, figW-txW-ppW-4*sp, butH]);
    hh = hh+butH+2*sp;
    
%%%%%%%%%%%%%%%%%%%%%% SHIFT panel %%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.frame2, 'Position', [sp/2, hh, floor(figW-sp), 3*butH+4*sp]);
    hh = hh+sp;    
butW = floor((figW-4*sp)/3);    
set(handles.ClearX, 'Position', [sp, hh, butW, butH]); 
set(handles.pbNorm, 'Position', [2*sp+butW, hh, butW, butH]); 
set(handles.pbToFirst, 'Position', [3*sp+2*butW, hh, butW, butH]); 
    hh = hh+butH+sp;
ppW = 70;
set(handles.text4, 'Position', [sp, hh, txW, butH]);
set(handles.SYdim, 'Position', [2*sp+txW, hh, ppW, butH]);
set(handles.ShiftY, 'Position', [2*sp+txW+ppW, hh, figW-txW-ppW-4*sp, butH]);
    hh = hh+butH+sp;
set(handles.text1, 'Position', [sp, hh, txW, butH]);
set(handles.SXdim, 'Position', [2*sp+txW, hh, ppW, butH]);
set(handles.ShiftX, 'Position', [2*sp+txW+ppW, hh, figW-txW-ppW-4*sp, butH]);
   hh = hh+butH+2*sp; 
   
butW = floor((figW-4*sp)/3);    
sph = 10;
sz = floor((butH-sph)/2);
set(handles.text6, 'Position', [sp, hh+sz, butW, sph]); 
set(handles.CollCur, 'Position', [2*sp+butW, hh, butW, butH]); 
set(handles.pbTitle, 'Position', [3*sp+2*butW, hh, butW, butH]);    

hh = hh+butH+2*sp; 
set(handles.popfiles, 'Position', [sp, hh, floor(figW-2*sp), butH]);


% --------------------------------------------------------------------
function varargout = figure1_DeleteFcn(h, eventdata, handles, varargin)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.figure1, handles.hh)']);


% --------------------------------------------------------------------
function varargout = popfiles_Callback(h, eventdata, handles, varargin)
UpdateInterface(handles);

% --------------------------------------------------------------------
function varargout = ShiftX_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

% shift =get(handles.Shift_X, 'Value');
% set(handles.Shift_X, 'Value', 0.5);
% num = get(handles.SXdim, 'Value');
% str = get(handles.SXdim, 'String');
% dim = str2num([str{num}]);
try
    CurVal = get(handles.ShiftX, 'Value');
    if ~isempty(CurVal)
        %     Val = CurVal + dim*(shift - 0.5)*2;
        %     set(handles.edit1, 'String', num2str(Val));
        num = get(handles.popfiles, 'Value');
        handles.trace{num}.ShiftX = CurVal;
        guidata(handles.figure1, handles);
        Plott(handles);
    end
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ZoomX_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

% shift =get(handles.Zoom_X, 'Value');
% set(handles.Zoom_X, 'Value', 0.5);
% num = get(handles.ZXdim, 'Value');
% str = get(handles.ZXdim, 'String');
% dim = str2num([str{num}]);
try
    CurVal = get(handles.ZoomX, 'Value');
%     Val = CurVal + dim*(shift - 0.5)*2;
%     set(handles.edit3, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ZoomX = CurVal;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ShiftY_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

% shift =get(handles.Shift_Y, 'Value');
% set(handles.Shift_Y, 'Value', 0.5);
% num = get(handles.SYdim, 'Value');
% str = get(handles.SYdim, 'String');
% dim = str2num([str{num}]);
try
    CurVal = get(handles.ShiftY, 'Value');
%     Val = CurVal + dim*(shift - 0.5)*2;
%     set(handles.edit2, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ShiftY = CurVal;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ZoomY_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

% shift =get(handles.Zoom_Y, 'Value');
% set(handles.Zoom_Y, 'Value', 0.5);
% num = get(handles.ZYdim, 'Value');
% str = get(handles.ZYdim, 'String');
% dim = str2num([str{num}]);
try
    CurVal = get(handles.ZoomY, 'Value');
%     Val = CurVal + dim*(shift - 0.5)*2;
%     set(handles.edit4, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ZoomY = CurVal;
    guidata(handles.figure1, handles);
    Plott(handles);    
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = DelCurr_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

num = get(handles.popfiles, 'Value');
endnum = length(handles.trace);

if endnum == 1, 
    Dell_Callback(h, eventdata, handles, varargin); 
    return; 
end

switch num,
case endnum
    coll = [1:endnum-1];
case 1
    coll = [2:endnum];
otherwise
    coll = [1:num-1, num + 1:endnum];
end

handles.trace = {handles.trace{coll}};
guidata(handles.figure1, handles);

str = get(handles.popfiles, 'String');
set(handles.popfiles, 'String', str(coll), 'Value', 1);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = Dell_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

handles.trace = {};
guidata(handles.figure1, handles);

set(handles.popfiles, 'String', '----------', 'Value', 1);
set(handles.text6, 'BackgroundColor', [0.75 0.75 0.75]);

Plott(handles);

% --------------------------------------------------------------------
function varargout = CollCur_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

num = get(handles.popfiles, 'Value');
Col = uisetcolor(handles.trace{num}.Data.ax.Color);
if length(Col)== 1, return; end

handles.trace{num}.Data.ax.Color = Col;
guidata(handles.figure1, handles);
set(handles.text6, 'BackgroundColor', Col);
Plott(handles);
   

% --------------------------------------------------------------------
function varargout = ClearX_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

for kk=1:length(handles.trace)
    handles.trace{kk}.ShiftX = 0;
    handles.trace{kk}.ShiftY = 0;
end

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = ClearY_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

if get(handles.popUnits, 'Value') == 1, zoom = 1;
else, zoom = 0;
end

for kk=1:length(handles.trace)
    handles.trace{kk}.ZoomX = 1;
    handles.trace{kk}.ZoomY = zoom;
end
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = editctl_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

num = get(handles.popfiles, 'Value');
handles.trace{num}.ShiftX = str2num(get(handles.ShiftX, 'String'));
handles.trace{num}.ShiftY = str2num(get(handles.ShiftY, 'String'));
handles.trace{num}.ZoomX = str2num(get(handles.ZoomX, 'String'));
handles.trace{num}.ZoomY = str2num(get(handles.ZoomY, 'String'));
handles.trace{num}.Data.ax.y = str2num(get(handles.eSecondDim, 'String'))';
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = AddButt_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
if isempty(hand),  
    msgbox('Can not find main box handles', 'Error');
    return;
end

[path,name,ext] = fileparts(hand.fname);
file_str = name;
String = {};
src = hand.src;
% if handles.process, [src.y, src.ax] = processing(src.y, src.ax);
% else
%     src.ax.filt = '';
%     src.ax.diff = '';
% end

counter = length(handles.trace)+1;
trace.Name = file_str;
trace.Data = src;
trace.Data.ax.Color = handles.Color(mod(counter, 7)+1,:);
if length(src.ax.y)==1 & src.ax.y(1)==0
    trace.Data.ax.y=counter;
end
trace.ShiftX = 0;
trace.ShiftY = 0;
trace.ZoomX = 1;
trace.ZoomY = 1;
trace.YSign = 1;
handles.trace{end+1} = trace;

% for the security reasons recreate the Str 
Str = {};
for ii=1:counter
    Str{ii} = [num2str(ii),':',handles.trace{ii}.Name];
end
set(handles.popfiles, 'String', Str, 'Value', counter);

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function Plott(ha)
h = guidata(ha.figure1);
hand = guidata(h.hh);
Num = length(h.trace);
hand.out = {};

for kk = 1:Num
    trace = h.trace{kk};
    if h.process, [trace.Data.y, trace.Data.ax] = processing(trace.Data.y, trace.Data.ax);
    else
        trace.Data.ax.filt = '';
        trace.Data.ax.diff = '';
    end
    
    hand.out{kk} = trace.Data;
    hand.out{kk}.ax.x = trace.Data.ax.x*trace.ZoomX + trace.ShiftX; % <-alsi 04.10.2005
    switch h.ZoomUnits
    case 'Norm'
        zoom = trace.ZoomY;
    case 'dB'
        zoom = 10^(trace.ZoomY/10);        
    end
    if safeget(trace, 'change2D', 0)
        sh = safeget(trace, 'ShiftY2D', 0)+trace.ShiftY;
        zm = safeget(trace, 'ZoomY2D', 1);
        hand.out{kk}.y = (trace.Data.y*zoom*trace.YSign - sh).*zm;
    else
        hand.out{kk}.y = trace.Data.y*zoom*trace.YSign +trace.ShiftY; % <-alsi 04.10.2005
    end
    
    hand.out{kk}.ax.s = 0;% trace.ShiftY <-alsi 04.10.2005
    hand.out{kk}.ax.dx= 0;% trace.ShiftX <-alsi 04.10.2005
    hand.out{kk}.title = ['Trace', num2str(kk)];
end

guidata(h.hh, hand);
[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function UpdateInterface(handles)
if length(handles.trace) == 0
    set(handles.ShiftX, 'Value', 0);
    set(handles.ShiftY, 'Value', 0);    
    set(handles.ZoomX, 'Value', 0);    
    set(handles.ZoomY, 'Value', 0);    
    set(handles.text6, 'BackgroundColor', [0.75,0.75,0.75]);
    set(handles.eSecondDim, 'String', '-');    
else
    num = get(handles.popfiles, 'Value');
    set(handles.ShiftX, 'Value', handles.trace{num}.ShiftX);
    set(handles.ShiftY, 'Value', handles.trace{num}.ShiftY);
    set(handles.ZoomX, 'Value', handles.trace{num}.ZoomX);
    set(handles.ZoomY, 'Value', handles.trace{num}.ZoomY);
    set(handles.text6, 'BackgroundColor', handles.trace{num}.Data.ax.Color);
    str = num2str(handles.trace{num}.Data.ax.y');
    set(handles.eSecondDim, 'String', str);    
end    

% --------------------------------------------------------------------
function varargout = figure1_CloseRequestFcn(h, eventdata, handles, varargin)
% Stub for CloseRequestFcn of the figure handles.figure1.
% try
%     hand = guidata(handles.hh);
%     axeshand = get(hand.figure1, 'CurrentAxes');
%     set(axeshand, 'NextPlot', 'replace');
%     plot(hand.ax.x, hand.y, 'Parent', axeshand);
%     axis(axeshand, 'tight');
%     guidata(hand.figure1, hand);
% end
closereq;

% --- Executes on selection change in SXdim.
function SXdim_Callback(hObject, eventdata, handles)
num = get(handles.SXdim, 'Value');
str = get(handles.SXdim, 'String');
dim = str2num([str{num}]);
set(handles.ShiftX, 'Step', dim);

% --- Executes on selection change in SYdim.
function SYdim_Callback(hObject, eventdata, handles)
num = get(handles.SYdim, 'Value');
str = get(handles.SYdim, 'String');
dim = str2num([str{num}]);
set(handles.ShiftY, 'Step', dim);

% --- Executes on selection change in ZXdim.
function ZXdim_Callback(hObject, eventdata, handles)
num = get(handles.ZXdim, 'Value');
str = get(handles.ZXdim, 'String');
dim = str2num([str{num}]);
set(handles.ZoomX, 'Step', dim);

% --- Executes on selection change in ZYdim.
function ZYdim_Callback(hObject, eventdata, handles)
num = get(handles.ZYdim, 'Value');
str = get(handles.ZYdim, 'String');
dim = str2num([str{num}]);
set(handles.ZoomY, 'Step', dim);

% -----------------------------------------------------
function pbToFirst_Callback(h, eventdata, handles)
if length(handles.trace) < 2, return; end

num = get(handles.popfiles, 'Value');
if num==1, return; end

handles.trace{num}.ZoomX = 1;
handles.trace{num}.ShiftX = 0;
y  = handles.trace{1}.Data.y;
y1 = handles.trace{num}.Data.y;
zmy = (max(real(y))-min(real(y)))/(max(real(y1))-min(real(y1)));
yprime = y1*zmy;
handles.trace{num}.ShiftY = min(real(y))-min(real(yprime));

switch handles.ZoomUnits
case 'Norm'
    handles.trace{num}.ZoomY = zmy;
case 'dB'
    handles.trace{num}.ZoomY = log10(zmy)*10;        
end

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function mFittAllFst_Callback(hObject, eventdata, handles)
if length(handles.trace) < 2, return; end

y  = handles.trace{1}.Data.y;

for num=2:length(handles.trace)
    handles.trace{num}.ZoomX = 1;
    handles.trace{num}.ShiftX = 0;
    y1 = handles.trace{num}.Data.y;
    if strcmp(get(handles.mPreserveZero, 'Checked'),'on')
        zmy = (max(real(y))-min(real(y)))/(max(real(y1))-min(real(y1)));
        handles.trace{num}.ShiftY = 0;
    else
        zmy = max(real(y))/max(real(y1));
        yprime = y1*zmy;
        handles.trace{num}.ShiftY = min(real(y))-min(real(yprime));
    end
    
    switch handles.ZoomUnits
    case 'Norm'
        handles.trace{num}.ZoomY = zmy;
    case 'dB'
        handles.trace{num}.ZoomY = log10(zmy)*10;        
    end
end

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function mSumAll_Callback(hObject, eventdata, handles)

str = {'off', 'on'};
stat = strcmp(get(handles.mSumAll, 'Checked'), 'on');
set(handles.mSumAll, 'Checked', str{~stat + 1});
handles.sumall = ~stat;

guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function mShowll_Callback(hObject, eventdata, handles)

str = {'off', 'on'};
stat = strcmp(get(handles.mSumAll, 'Checked'), 'on');
set(handles.mSumAll, 'Checked', str{~stat + 1});
handles.sumall = ~stat;

guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function mWaterfall_Callback(hObject, eventdata, ha)
if length(ha.trace) == 0, return; end

wfc = inputdlg('wf','Enter the Waterfall parameter',1,{'1'});
wfc = str2num(wfc{1});
h = guidata(ha.figure1);
hand = guidata(h.hh);
Num = length(h.trace);

maxy = -1E19; miny = 1E19;

for kk = 1:Num
    maxy = max([maxy,max(max(real(h.trace{kk}.Data.y.*h.trace{kk}.ZoomY+h.trace{kk}.ShiftY)))]);
    miny = min([miny,min(min(real(h.trace{kk}.Data.y.*h.trace{kk}.ZoomY+h.trace{kk}.ShiftY)))]);
end

wf = (maxy-miny)*wfc;

switch hObject
case ha.mWaterfall,
    for kk = 1:Num
        h.trace{kk}.ShiftY = h.trace{kk}.ShiftY - (kk-1)*wf;
    end
case ha.mWaterfall2D,
    for kk = 1:Num
        [sz1,sz2] = size(h.trace{kk}.Data.y);
        sh = [1:sz2]-1;
        %h.trace{kk}.Data.y = h.trace{kk}.Data.y - sh(ones(sz1,1),:)*wf;
        h.trace{kk}.ShiftY2D = sh(ones(sz1,1),:)*wf; 
        h.trace{kk}.change2D = 1;
    end
    
end
guidata(ha.figure1, h)
UpdateInterface(h);
Plott(h);

% --------------------------------------------------------------------
function pbInvY_Callback(hObject, eventdata, handles)
if length(handles.trace) == 0, return; end

num = get(handles.popfiles, 'Value');
handles.trace{num}.YSign =  -handles.trace{num}.YSign;

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function mPrAgr_Callback(hObject, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(handles.mPrAgr, 'Checked'), 'on');
set(handles.mPrAgr, 'Checked', str{~stat + 1});
handles.process = ~stat;
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function popUnits_Callback(h, eventdata, handles)
str = {'Norm', 'dB'};
num = get(handles.popUnits, 'Value');
ZoomUnits = str{num};

if strcmp(ZoomUnits, handles.ZoomUnits), return; end;
handles.ZoomUnits = ZoomUnits;

for num=1:length(handles.trace)
    switch handles.ZoomUnits
    case 'Norm'
        handles.trace{num}.ZoomY = 10^(handles.trace{num}.ZoomY/10);
    case 'dB'
        handles.trace{num}.ZoomY = log10(handles.trace{num}.ZoomY)*10;        
    end
end

guidata(handles.figure1, handles);
UpdateInterface(handles);

% --------------------------------------------------------------------
function m1Dto2D_Callback(hObject, eventdata, handles)
Num = length(handles.trace);
if Num == 0, return; end

Out = [];
dim2 = [];
xx = handles.trace{1}.Data.ax.x*handles.trace{1}.ZoomX + handles.trace{1}.ShiftX;
minax = xx(1);
maxax = xx(end);
npax = length(xx);
axdx  = (maxax - minax)/length(handles.trace{1}.Data.ax.x);
notsameax = 0;
for kk = 2:Num
    xx = handles.trace{kk}.Data.ax.x*handles.trace{kk}.ZoomX + handles.trace{kk}.ShiftX;
    notsameax = notsameax + minax-xx(1) + maxax-xx(end) + npax - length(xx);
    if xx(1)   > minax, minax = xx(1); end
    if xx(end) < maxax, maxax = xx(end); end    
end
npoints = floor((maxax - minax)/axdx + 0.5);

if notsameax
    newax = linspace(minax, maxax, npoints)';
else
    newax = xx(:);
end
for kk = 1:Num
    if handles.process, [handles.trace{kk}.Data.y, handles.trace{kk}.Data.ax] = ...
            processing(handles.trace{kk}.Data.y, handles.trace{kk}.Data.ax);
    else
        handles.trace{kk}.Data.ax.filt = '';
        handles.trace{kk}.Data.ax.diff = '';
    end
  x = handles.trace{kk}.Data.ax.x*handles.trace{kk}.ZoomX + handles.trace{kk}.ShiftX;
  y = handles.trace{kk}.Data.y*handles.trace{kk}.ZoomY*handles.trace{kk}.YSign + handles.trace{kk}.ShiftY;
%   if notsameax
      Out(:, kk) = spline(x, y, newax);
%   else
%       Out(:, kk) = y;      
%   end
    dim2(end+1,1) = handles.trace{kk}.Data.ax.y;
end

ax = handles.trace{1}.Data.ax;
if isfield(handles.trace{1}.Data, 'dsc')
    dsc = handles.trace{1}.Data.dsc;
else
    dsc = [];
end
trace.Data = {};
trace.Data.ax = ax;
trace.Data.ax.x = newax;
trace.Data.ax.y = dim2;
trace.Data.ax.ylabel = 'trace,';
trace.Data.ax.c = 'b';
trace.Data.ax.s = 0;
trace.Data.ax.dx = 0;
trace.Data.ax.filt = '';
trace.Data.ax.diff = '';
trace.Data.ax.Color = handles.Color(1,:);


trace.Data.dsc = dsc;

trace.Data.y = Out;
trace.Data.title = 'Out';
trace.ShiftX = 0;
trace.ShiftY = 0;
trace.ZoomX = 1;
trace.ZoomY = 1;
trace.YSign = 1;
trace.Name = safeget(handles.trace{1}, 'Name', '');
handles.trace = {trace};

set(handles.popfiles, 'String', 'Group', 'Value', 1);
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function m2Dto1D_Callback(hObject, eventdata, handles)
if length(handles.trace) < 1, return; end

num = get(handles.popfiles, 'Value');
trace = handles.trace{num};
DelCurr_Callback(hObject, eventdata, handles);
handles = guidata(handles.figure1);

str = get(handles.popfiles, 'String');
endnum = length(handles.trace);
if ~endnum, str = {}; end;

counter  = size(trace.Data.y,2);
if isfield(trace.Data, 'dsc')
    dsc = trace.Data.dsc;
else
    dsc = [];
end

for ii=1:counter
   handles.trace{end+1}.Data.ax = trace.Data.ax;
   handles.trace{end}.Data.dsc = dsc;
   handles.trace{end}.Data.y  = trace.Data.y(:,ii);
   handles.trace{end}.Data.ax.y = trace.Data.ax.y(ii);
   handles.trace{end}.Data.title = [trace.Data.ax.title,'::',num2str(ii)];
   handles.trace{end}.Data.ax.Color = handles.Color(mod(ii,7)+1,:);
   str{end+1} = trace.Data.ax.title;
   handles.trace{end}.ShiftX = trace.ShiftX;
   handles.trace{end}.ShiftY = trace.ShiftY;
   handles.trace{end}.ZoomX  = trace.ZoomX;
   handles.trace{end}.ZoomY  = trace.ZoomY;
   handles.trace{end}.YSign  = trace.YSign;
   handles.trace{end}.Name   = [num2str(endnum+ii),':',trace.Data.ax.title,'::',num2str(ii)];
end

set(handles.popfiles, 'String', str, 'Value', 1);
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);
return

% --------------------------------------------------------------------
function varargout = pbChange_Callback(h, eventdata, handles, varargin)
Num = length(handles.trace);
if Num < 1, AddButt_Callback(h, eventdata, handles, varargin);
    return
end

hand = guidata(handles.hh);
if isempty(hand), msgbox('Can not find main box handles', 'Error');
    return;
end

Num = get(handles.popfiles, 'Value');
str = get(handles.popfiles, 'String');
trace = handles.trace{Num};

[path,name,ext] = fileparts(hand.fname);
file_str = name;

src = hand.src;
handles.trace{Num}.Name = file_str;
handles.trace{Num}.Data = src;
handles.trace{Num}.Data.ax.Color = trace.Data.ax.Color;

if iscell(str), str{Num} = [num2str(Num),':',handles.trace{Num}.Name];
else, str = [num2str(Num),':',handles.trace{Num}.Name];       
end

set(handles.popfiles, 'String', str);
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = edit7_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = mSortIndex_Callback(h, eventdata, handles, varargin)

if length(handles.trace) < 2, return; end

yy = [];
Str = {};
for kk=1:length(handles.trace), yy(end+1)=handles.trace{kk}.Data.ax.y;end;
[yyy,idx] = sort(yy);

trace = handles.trace;
for kk=1:length(handles.trace), 
    handles.trace{kk}=trace{idx(kk)};
    Str{kk} = [num2str(kk),':',handles.trace{kk}.Name];
end;

set(handles.popfiles, 'String', Str);
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = mPreserveZero_Callback(h, eventdata, handles, varargin)
str = {'off', 'on'};
stat = strcmp(get(handles.mPreserveZero, 'Checked'), 'on');
set(handles.mPreserveZero, 'Checked', str{~stat + 1});

% --------------------------------------------------------------------
function varargout = mNormSlices_Callback(h, eventdata, handles, varargin)
num = get(handles.popfiles, 'Value');
parentmenu = get(h, 'Parent');
if isempty(safeget(handles, 'mDisChandes', [])),
    handles.mDisChandes = uimenu(parentmenu, 'Label', 'Discard ..2D', 'Tag', 'mDisChandes', ...
        'callback', 'compbox(''mDisChanges_Callback'', guidata(gcbo))');
else
    set(handles.mDisChandes, 'Visible', 'on');
end
    handles.trace{num}.ZoomY2D = handles.trace{num}.Data.y*0;
    handles.trace{num}.ShiftY2D = handles.trace{num}.Data.y*0;
for slice=1:size(handles.trace{num}.Data.y,2)
    maxy = max(real(handles.trace{num}.Data.y(:,slice)));
    miny = min(real(handles.trace{num}.Data.y(:,slice)));
    
    if strcmp(get(handles.mPreserveZero, 'Checked'),'on')
        if maxy > 0 & miny > 0, miny=0;
        elseif maxy < 0 & miny < 0, maxy=0;
        end
    end
    
    switch handles.ZoomUnits
    case 'Norm'
        handles.trace{num}.ZoomY = 1;
    case 'dB'
        handles.trace{num}.ZoomY = 0;        
    end
    handles.trace{num}.ShiftY = 0;
    
    zmy = 1/(maxy-miny);
    %handles.trace{num}.Data.y(:,slice) = (handles.trace{num}.Data.y(:,slice) - miny) * zmy;
    handles.trace{num}.ZoomY2D(:, slice) = zmy(ones(size(handles.trace{num}.Data.y, 1), 1), 1);
    handles.trace{num}.ShiftY2D(:, slice) = miny(ones(size(handles.trace{num}.Data.y, 1), 1), 1);  
    handles.trace{num}.change2D = 1;
end

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);
% --------------------------------------------------------------------
function varargout = pbNorm_Callback(h, eventdata, handles, varargin)
num = get(handles.popfiles, 'Value');

maxy = max(max(real(handles.trace{num}.Data.y)));
miny = min(min(real(handles.trace{num}.Data.y)));

zmy = 1/(maxy-miny);
if strcmp(get(handles.mPreserveZero, 'Checked'),'on')
    handles.trace{num}.ShiftY = 0;
else
    handles.trace{num}.ShiftY = -miny*zmy;
end 

switch handles.ZoomUnits
case 'Norm'
    handles.trace{num}.ZoomY = zmy;
case 'dB'
    handles.trace{num}.ZoomY = log10(zmy)*10;        
end

guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);


% --------------------------------------------------------------------
function varargout = mNormalize_Callback(h, eventdata, handles, varargin)

for num=1:length(handles.trace)
    maxy = max(max(real(handles.trace{num}.Data.y)));
    miny = min(min(real(handles.trace{num}.Data.y)));
    
    zmy = 1/(maxy-miny);
    if strcmp(get(handles.mPreserveZero, 'Checked'),'on')
        handles.trace{num}.ShiftY = 0;
    else
        handles.trace{num}.ShiftY = -miny*zmy;
    end 
    
    switch handles.ZoomUnits
    case 'Norm'
        handles.trace{num}.ZoomY = zmy;
    case 'dB'
        handles.trace{num}.ZoomY = log10(zmy)*10;        
    end
end
guidata(handles.figure1, handles);
UpdateInterface(handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = mFittAddToOut_Callback(h, eventdata, handles, varargin)
h = guidata(handles.figure1);
hand = guidata(h.hh);

kk = get(handles.popfiles, 'Value');
trace = h.trace{kk};
if h.process, [trace.Data.y, trace.Data.ax] = processing(trace.Data.y, trace.Data.ax);
else
    trace.Data.ax.filt = '';
    trace.Data.ax.diff = '';
end
hand.out{end+1}.ax = trace.Data.ax;
hand.out{end}.ax.x = trace.Data.ax.x*trace.ZoomX;
switch h.ZoomUnits
case 'Norm'
    zoom = trace.ZoomY;
case 'dB'
    zoom = 10^(trace.ZoomY/10);        
end
hand.out{end}.y = trace.Data.y*zoom*trace.YSign;
hand.out{end}.ax.s = trace.ShiftY;
hand.out{end}.ax.dx= trace.ShiftX;
hand.out{end}.title = ['Trace', num2str(kk)];

guidata(h.hh, hand);
[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);


% --------------------------------------------------------------------
function varargout = pbTitle_Callback(h, eventdata, handles, varargin)
Str = get(handles.popfiles, 'String');
num = get(handles.popfiles, 'Value');

newtitle = inputdlg('Change Data Title', 'Data Title', 1, {handles.trace{num}.Name});
Str = SetTitle(handles.popfiles, newtitle{1}, num);

handles.trace{num}.Name = newtitle{1};
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function Str = SetTitle(combo, newtitle, num)
Str = get(combo, 'String');
Str{num} = [num2str(num),':',newtitle];
set(combo, 'String', Str);

% --------------------------------------------------------------------
function mDisChanges_Callback(handles)
Num = length(handles.trace);
for kk = 1:Num
    handles.trace{kk}.ShiftY2D = 0;
    handles.trace{kk}.ZoomY2D = 1;    
    handles.trace{kk}.ShiftX2D = 0;
    handles.trace{kk}.ZoomX2D = 1;    
end
guidata(handles.figure1, handles);
if ~isempty(safeget(handles, 'mDisChandes', [])),
    set(handles.mDisChandes, 'Visible', 'off');
end
Plott(handles);
