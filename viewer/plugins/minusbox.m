function varargout = minusbox(varargin)
% MINUSBOX Application M-file for minusbox.fig
%    FIG = MINUSBOX launch minusbox GUI.
%    MINUSBOX('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 08-Nov-2013 10:18:16

if nargin == 0  % LAUNCH GUI
    token = 'MINUSBOX_open';
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
    handles.ShiftX = 0;
    handles.set(1).ShiftY = 0;
    handles.ZoomX = 1;
    handles.set(1).ZoomY = 1;    
    hadnles.hh = [];
    handles.YSign = 1;
    handles.ShRes = 0;
    handles.valShRes = 0;
    handles.process = 0;
    handles.operator = '-';
    handles.Ref = [];
    handles.view = [1, 1, 1];
    handles.fileupdate = 0;
    handles.SubMean = 0;
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

% --- Executes during object deletion, before destroying properties.
function varargout = figure1_DeleteFcn(h, eventdata, handles, varargin)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.figure1, handles.hh)']);


% --------------------------------------------------------------------
function varargout = slider_Callback(h, eventdata, handles, varargin)

slice = getslicestring(handles);

switch h
case handles.slShiftX
    handedit = handles.edShiftX;
    handpm = handles.pmShiftX;
    varname = 'ShiftX';
    
case handles.slShiftY
    handedit = handles.edShiftY;
    handpm = handles.pmShiftY;
    varname = ['set(ii).ShiftY'];
    
case handles.slZoomX
    handedit = handles.edZoomX;
    handpm = handles.pmZoomX;
    varname = 'ZoomX';
    
case handles.slZoomY
    handedit = handles.edZoomY;
    handpm = handles.pmZoomY;
    varname = ['set(ii).ZoomY'];
end

shift =get(h, 'Value');
set(h, 'Value', 0.5);
num = get(handpm, 'Value');
str = get(handpm, 'String');
dim = str2num([str{num}]);
try
    CurVal = str2num(get(handedit, 'String'));
    Val = CurVal + dim*(shift - 0.5)*2;
    set(handedit, 'String', num2str(Val));
    for ii=slice
        eval(['handles.',varname,' = Val;']);
    end
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = edit_Callback(h, eventdata, handles, varargin)

slice = getslicestring(handles);

switch h
case handles.edShiftX
    varname = 'ShiftX';
case handles.edShiftY
    varname = ['set(ii).ShiftY'];
case handles.edZoomX
    varname = 'ZoomX';
case handles.edZoomY
    varname = ['set(ii).ZoomY'];
end

for ii=slice
    ss = strfind(lower(get(h, 'String')), 'auto');
    if isempty(ss)
        eval(['handles.',varname,' = str2num(get(h, ''String''));']);
    else
        eval(['handles.',varname,' = i']);
    end
end
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = pbRef_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
handles.Ref = hand.src;
n = size(hand.src.y, 2);
str{1} = 'All';
for ii=1:n, str{ii+1}=num2str(ii); end
set(handles.pmSlice, 'String', str, 'Value', 1);
while length(handles.set) < n, handles.set(end+1) = handles.set(1); end    
set(h, 'BackgroundColor', [0.64 0.8 0.64]);
guidata(handles.figure1, handles);
pmSlice_Callback(h, eventdata, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = pbShZero_Callback(h, eventdata, handles, varargin)
handles.ShiftX = 0;
handles.set(1).ShiftY = 0;
set(handles.edShiftX, 'String', '0');
set(handles.edShiftY, 'String', '0');
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = pbZmZero_Callback(h, eventdata, handles, varargin)
handles.ZoomX = 1;
handles.set(1).ZoomY = 1;
set(handles.edZoomX, 'String', '1');
set(handles.edZoomY, 'String', '1');
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = checkbox_Callback(h, eventdata, handles, varargin)
hands = [handles.chkPluse, handles.chkMinus, handles.chkDiff, handles.chkMult, handles.chkRef];
set(hands, 'Value', 0);
set(h, 'Value', 1);
handles.operator = get(h, 'String');
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = popupmenu_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function Plott(handles)
h = guidata(handles.figure1);

isref = ~get(h.chkRef, 'Value');
if isref & isempty(h.Ref), error('Please set the reference.'); end

hand = guidata(h.hh);
hand.out = {};

tmpSec = hand.src;

if ~isref, 
    tmpRef = tmpSec; 
else
    tmpRef = h.Ref;
end

if h.process, 
    [ tmpRef.y,  tmpRef.ax] = processing( h.Ref.y,  h.Ref.ax);
    [ tmpSec.y,  tmpSec.ax] = processing( hand.src.y,  hand.src.ax);
else
     tmpRef.ax.filt = '';
     tmpRef.ax.diff = '';
     tmpSec.ax.filt = '';
     tmpSec.ax.diff = '';     
end
if h.SubMean
    tmpRef.y = tmpRef.y-mean(tmpRef.y);
    tmpSec.y = tmpSec.y-mean(tmpSec.y);
end

tmpRef.ax.c = [0, 0, 1];
tmpRef.title = ['Ref'];

tmpSec.ax.x =tmpSec.ax.x*h.ZoomX + h.ShiftX;
if imag(h.set(1).ShiftY) % automatic shift
    for ii=1:size(tmpSec.y,2);
        tmpSec.y(:,ii) = tmpSec.y(:,ii)*h.set(ii).ZoomY + tmpRef.y(1)-tmpSec.y(1,ii)*h.set(ii).ZoomY ;
    end
else
    for ii=1:size(tmpSec.y,2);
        tmpSec.y(:,ii) = tmpSec.y(:,ii)*h.set(ii).ZoomY + h.set(1).ShiftY;
    end
end
tmpSec.ax.Color = [1, 0, 0];
tmpSec.ax.s = 0;
tmpSec.ax.dx= 0;
tmpSec.title = ['Op 2'];


if isref
    br = 0;
    if (length(tmpRef.ax.x)==length(tmpSec.ax.x))
        if sum(tmpRef.ax.x - tmpSec.ax.x)
            br = 1;
        end
    else
        br=1;
    end
    if br
        x1 = max([tmpRef.ax.x(1), tmpSec.ax.x(1)]);
        x2 = min([tmpRef.ax.x(end), tmpSec.ax.x(end)]);
        npoints  = size(tmpRef.ax.x, 1);
        xax = linspace(x1, x2, npoints).';
        for ii=1:size(tmpRef.y, 2);
            tmp1(:,ii) = spline(tmpRef.ax.x, tmpRef.y(:,ii), xax);
            tmp2(:,ii) = spline(tmpSec.ax.x, tmpSec.y(:,ii), xax);
        end
    else
        tmp1 = tmpRef.y;
        tmp2 = tmpSec.y;
        xax = tmpRef.ax.x;
    end
    eval(['y = tmp1', h.operator, 'tmp2;']);
else
    y = tmpSec.y;
    xax = tmpSec.ax.x;
end

tmpRes = tmpRef;
tmpRes.ax.x =xax;
tmpRes.y = y*h.YSign;
tmpRes.ax.Color = [0, 1, 0]*0.8;
tmpRes.ax.s = handles.ShRes*handles.valShRes;
tmpRes.ax.dx= 0;
tmpRes.title = ['Result'];
if handles.view(1), hand.out{end+1} = tmpRef; end
if handles.view(2) & ~get(h.chkRef, 'Value'), hand.out{end+1} = tmpSec; end
if handles.view(3), hand.out{end+1} = tmpRes; end

guidata(h.hh, hand);
[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);


% --------------------------------------------------------------------
function varargout = chkRevRes_Callback(h, eventdata, handles, varargin)
handles.YSign = handles.YSign*(-1);
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = pbToRef_Callback(h, eventdata, handles, varargin)
if isempty(handles.Ref), return; end

hand = guidata(handles.hh);
y  = handles.Ref.y;
y1 = hand.src.y;
kk = (max(real(y))-min(real(y)))/(max(real(y1))-min(real(y1)));
zmx = 1;
zmy = kk;
yprime = y1*kk;
shx = 0;
shy = min(real(y))-min(real(yprime));

set(handles.edShiftX, 'String', num2str(shx));
set(handles.edShiftY, 'String', num2str(shy));
set(handles.edZoomX, 'String', num2str(zmx));
set(handles.edZoomY, 'String', num2str(zmy));    

handles.ShiftX = shx;
handles.set(1).ShiftY = shy;
handles.ZoomX = zmx;
handles.set(1).ZoomY = zmy;
guidata(handles.figure1, handles);

Plott(handles);

% --------------------------------------------------------------------
function varargout = chkSftRes_Callback(h, eventdata, handles, varargin)
handles.ShRes = ~handles.ShRes;
guidata(handles.figure1, handles);
Plott(handles);

% --------------------------------------------------------------------
function varargout = edSftRes_Callback(h, eventdata, handles, varargin)
try
    handles.valShRes = str2num(get(h, 'String'));
    marker = 1;
catch
    set(h, 'String', handles.valShRes);
    marker = 0;
end
if marker
    guidata(handles.figure1, handles);
    Plott(handles);
end

% --------------------------------------------------------------------
function mView_Callback(h, eventdata, handles, varargin)
c = {'off', 'on'};
stat = strcmp(get(h, 'Checked'), 'on');
set(h, 'Checked', [c{~stat + 1}]);

switch h
case handles.mShRef
    handles.view(1) = ~handles.view(1);
case handles.mShTrace
    handles.view(2) = ~handles.view(2);    
case handles.mShRes   
    handles.view(3) = ~handles.view(3);    
end
guidata(handles.figure1, handles);
Plott(handles);

function mSub_Callback(h, eventdata, handles, varargin)
c = {'off', 'on'};
stat = strcmp(get(h, 'Checked'), 'on');
set(h, 'Checked', [c{~stat + 1}]);
handles.SubMean = ~stat;
guidata(handles.figure1, handles);
Plott(handles);
% --------------------------------------------------------------------
function mAbout_Callback(h, eventdata, handles, varargin)
msgbox({'MINUSBOX'; 'by Boris Epel and Alexey Silakov, 2003-04'; ''; ...
        'Result = Reference +-*/ CurData(x*zoomx + shiftx)*zoomy+shifty';...
        'or Result = CurData(x*zoomx + shiftx)*zoomy+shifty';
        'If abscissas of the Reference and Current Data are not identical,'; ' ''spline'' procedure is used to interpolate the data';...
        'The Result can be inverted and shifted along Y axis';
        ''; ...
        'IMPORTANT: The Ref. and Cur. Data must have same number of columns'}, ...
    'About', 'help');

% --------------------------------------------------------------------
function FileUpdate(handles)
if handles.fileupdate
    Plott(handles);
end

% --------------------------------------------------------------------
function mFileUp_Callback(h, eventdata, handles, varargin)
c = {'off', 'on'};
stat = strcmp(get(h, 'Checked'), 'on');
set(h, 'Checked', [c{~stat + 1}]);

handles.fileupdate = ~handles.fileupdate;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function pmSlice_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --------------------------------------------------------------------
function pmSlice_Callback(h, eventdata, handles)
n = get(handles.pmSlice, 'Value');
if n > 1, n = n - 1; end;

set(handles.edShiftY, 'String', handles.set(n).ShiftY);
set(handles.edZoomY, 'String', handles.set(n).ZoomY);    

% --------------------------------------------------------------------
function sl = getslicestring(handles)
n = get(handles.pmSlice, 'Value');
if n==1, sl=1:length(get(handles.pmSlice, 'String'));
else sl = n-1;
end
