function varargout = mathplugin(varargin)
% MATHPLUGIN Application M-file for mathplugin.fig
%    FIG = MATHPLUGIN launch mathplugin GUI.
%    MATHPLUGIN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 21-Nov-2005 15:10:59
% Alexey Silakov & Boris Epel 22-Mar-2004, MPI
% alsi 24-Mar-2004, MPI

if nargin == 0  % LAUNCH GUI
    token = 'MATHPLUGIN_open';
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
    handles.FstD.st = 1;
    handles.FstD.end = 1;
    handles.SndD.st = 1;
    handles.SndD.end = 1;
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
% --------------------------------------------------------------------
function FileUpdate(handles);
handles = guidata(handles.figure1);

if ~isfield(handles, 'hh'), disp('No handles from MainFigure'); return; end
hand = guidata(handles.hh);

handles.FstD.st = 1;
set(handles.ed1Dst, 'String', 1);

handles.SndD.st = 1;
set(handles.ed2Dst, 'String', 1);

handles.FstD.end = size(hand.src.y, 1);
set(handles.ed1Dend, 'String', size(hand.src.y, 1));

handles.SndD.end = size(hand.src.y, 2);
set(handles.ed2Dend, 'String', size(hand.src.y, 2));
return;
% --------------------------------------------------------------------
function varargout = Exec_Callback(h, eventdata, handles, varargin)
Plott(handles);

% ----------------------------------------------------
function sl1Dst_Callback(hObject, eventdata, handles)
shift =get(handles.sl1Dst, 'Value');
set(handles.sl1Dst, 'Value', 0.5);
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
try
    CurVal = str2num(get(handles.ed1Dst, 'String'));
    Val = min(floor(CurVal + (shift - 0.5)*2), size(hand.src.y, 1));
    if Val <1, Val = 1; end
    set(handles.ed1Dst, 'String', num2str(Val));
    handles.FstD.st = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% ----------------------------------------------------
function sl2Dst_Callback(hObject, eventdata, handles)
shift =get(handles.sl2Dst, 'Value');
set(handles.sl2Dst, 'Value', 0.5);
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
try
    CurVal = str2num(get(handles.ed2Dst, 'String'));
    Val = min(floor(CurVal + (shift - 0.5)*2), size(hand.src.y, 2));
    if Val <1, Val = 1; end
    set(handles.ed2Dst, 'String', num2str(Val));
    handles.SndD.st = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% ----------------------------------------------------
function sl1Dend_Callback(hObject, eventdata, handles)
shift =get(handles.sl1Dend, 'Value');
set(handles.sl1Dend, 'Value', 0.5);
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

try
    CurVal = str2num(get(handles.ed1Dend, 'String'));
    Val = min(floor(CurVal + (shift - 0.5)*2), size(hand.src.y, 1));
    if Val <1, Val = 1; end    
    set(handles.ed1Dend, 'String', num2str(Val));
    handles.FstD.end = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% ----------------------------------------------------
function sl2Dend_Callback(hObject, eventdata, handles)
set(handles.sl2Dend, 'Value', 0.5);
shift =get(handles.sl2Dend, 'Value');
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

try
    CurVal = str2num(get(handles.ed2Dend, 'String'));
    Val = min(floor(CurVal + (shift - 0.5)*2), size(hand.src.y, 2));
    if Val <1, Val = 1; end    
    set(handles.ed2Dend, 'String', num2str(Val));
    handles.SndD.end = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% ----------------------------------------------------
function ed1Dst_Callback(hObject, eventdata, handles)
num = str2num(get(handles.ed1Dst, 'String'));
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
if isempty(num)
    set(handles.ed1Dst, 'String', handles.FstD.st);
    return;
else
    if floor(num)<1, num = 1; end
    handles.FstD.st = min(floor(num), size(hand.src.y, 1));
end
set(handles.ed1Dst, 'String', handles.FstD.st);
guidata(handles.figure1, handles);
if get(handles.tgCutData, 'Value')
    Plott(handles);
end

% ----------------------------------------------------
function ed2Dst_Callback(hObject, eventdata, handles)
num = str2num(get(handles.ed2Dst, 'String'));
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
if isempty(num)
    set(handles.ed2Dst, 'String', handles.SndD.st);
    return;
else
    if floor(num)<1, num = 1; end
    handles.SndD.st = min(floor(num), size(hand.src.y, 2));
end
set(handles.ed2Dst, 'String', handles.SndD.st);
guidata(handles.figure1, handles);
if get(handles.tgCutData, 'Value')
    Plott(handles);
end

% ----------------------------------------------------
function ed1Dend_Callback(hObject, eventdata, handles)
num = str2num(get(handles.ed1Dend, 'String'));
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

if isempty(num)
    set(handles.ed1Dend, 'String', handles.FstD.end);
    return;
else
    if floor(num)<1, num = 1; end
    handles.FstD.end = min(floor(num), size(hand.src.y, 1));
end
set(handles.ed1Dend, 'String', handles.FstD.end);
guidata(handles.figure1, handles);
if get(handles.tgCutData, 'Value')
    Plott(handles);
end

% ----------------------------------------------------
function ed2Dend_Callback(hObject, eventdata, handles)
num = str2num(get(handles.ed2Dend, 'String'));
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

if isempty(num)
    set(handles.ed2Dend, 'String', handles.SndD.end);
    return;
else
    if floor(num)<1, num = 1; end
    handles.SndD.end = min(floor(num), size(hand.src.y, 2));
end
set(handles.ed2Dend, 'String', handles.SndD.end);
guidata(handles.figure1, handles);
if get(handles.tgCutData, 'Value')
    Plott(handles);
end

% ----------------------------------------------------
function tgCutData_Callback(hObject, eventdata, handles)
guidata(handles.figure1, handles);
Plott(handles);

% ----------------------------------------------------
function Plott(handles)
handles = guidata(handles.figure1);

hand = guidata(handles.hh);
hand.out = {};
tmp = hand.src;

if get(handles.tgCutData, 'Value')
    str1 = ['min(',num2str(handles.FstD.st), ', end) : min(', num2str(handles.FstD.end),', end)'];
    if size(hand.src.y, 2)>1
        if get(handles.cbExpression2Dim, 'Value')
            str2 = get(handles.eCut2Dim, 'String');
        else
            str2 = ['min(',num2str(handles.SndD.st), ', end) : min(', num2str(handles.SndD.end),', end)'];    
        end
    else
        str2 = ':';
    end
    eval(['tmp.y = tmp.y(',str1,',',str2,');']);
    eval(['tmp.ax.x = tmp.ax.x(',str1,', 1);']);
    eval(['tmp.ax.y = tmp.ax.y(',str2,', 1);']);
end

if get(handles.pbTranspose, 'Value')
    tmp.y = tmp.y.';
    ax = tmp.ax;
    tmp.ax.x = ax.y;
    if isfield(ax, 'ylabel')
        tmp.ax.xlabel = ax.ylabel;
    end
    tmp.ax.y = ax.x;
    if isfield(ax, 'xlabel')
        tmp.ax.ylabel = ax.xlabel;
    end
end

if get(handles.tgAbs, 'Value')
    tmp.y = abs(tmp.y);
end

if get(handles.tbSum, 'Value')
%     tmp.y = sum(tmp.y,2);
    tmp.y = mean(tmp.y,2);
    tmp.ax.y = 1;
end

[sz1,sz2]=size(tmp.y);
if get(handles.pbPhase, 'Value')
    phase = str2num(get(handles.ePhase, 'String'));
    phase1 = str2num(get(handles.ePhase1, 'String'));
    phase1piv = str2num(get(handles.ePhase1pivot, 'String'));
    phase1sh = str2num(get(handles.ePhase1sh, 'String'));
    for jj=1:sz2
        tmp.y(:,jj) = tmp.y(:,jj).*exp(-i*...
            (phase/180*pi+(phase1+phase1sh*(jj-1))*pi/180*([1:sz1]'-phase1piv)));
    end
    ax = tmp.ax;
end

hand.out{1} = tmp;
guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function varargout = phase_Callback(h, eventdata, handles, varargin)
guidata(handles.figure1, handles);
switch h
case handles.slPhase
    shift =get(handles.slPhase, 'Value');
    set(handles.slPhase, 'Value', 0.5);
    CurVal=str2num(get(handles.ePhase, 'String'));
    cc = safeget(handles, 'PhaseFactor',1.);
    Val = max(min(CurVal + (shift - 0.5)*2*cc, 180), -180);
    set(handles.ePhase, 'String', num2str(Val))
end
Plott(handles);

% --------------------------------------------------------------------
function varargout = slPhase_Callback(h, eventdata, handles, varargin)
guidata(handles.figure1, handles);
shift =get(handles.slPhase, 'Value');
set(handles.slPhase, 'Value', 0.5);
CurVal=str2num(get(handles.ePhase1, 'String'));
cc = safeget(handles, 'PhaseFactor',1.);
Val = max(min(CurVal + (shift - 0.5)*2*cc, 180), -180);
set(handles.ePhase1, 'String', num2str(Val))
Plott(handles);

% --------------------------------------------------------------------
function varargout = ePhase1_Callback(h, eventdata, handles, varargin)
Plott(handles);

% --------------------------------------------------------------------
function varargout = slPhase1_Callback(h, eventdata, handles, varargin)
guidata(handles.figure1, handles);
shift =get(handles.slPhase1, 'Value');
set(handles.slPhase1, 'Value', 0.5);
CurVal=str2num(get(handles.ePhase1, 'String'));
cc = safeget(handles, 'PhaseFactor',1.);
Val = max(min(CurVal + (shift - 0.5)*2*cc, 180), -180);
set(handles.ePhase1, 'String', num2str(Val))
Plott(handles);

% --------------------------------------------------------------------
function varargout = ePhase1pivot_Callback(h, eventdata, handles, varargin)
Plott(handles);

% --------------------------------------------------------------------
function varargout = ePhase1sh_Callback(h, eventdata, handles, varargin)
Plott(handles);

% --------------------------------------------------------------------
function varargout = slPhase2_Callback(h, eventdata, handles, varargin)
shift =get(handles.slPhase2, 'Value');
set(handles.slPhase2, 'Value', 0.5);
CurVal=str2num(get(handles.ePhase1sh, 'String'));
cc = safeget(handles, 'PhaseFactor',1.);
Val = max(min(CurVal + (shift - 0.5)*.2*cc, 180), -180);
set(handles.ePhase1sh, 'String', num2str(Val))
Plott(handles);


% --------------------------------------------------------------------
function varargout = pb01_Callback(h, eventdata, handles, varargin)
if get(handles.pb01, 'Value'), handles.PhaseFactor = 0.1; set(handles.pb10, 'Value',0); 
else handles.PhaseFactor = 1; end
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = pb10_Callback(h, eventdata, handles, varargin)
if get(handles.pb10, 'Value'), handles.PhaseFactor = 10; set(handles.pb01, 'Value',0); 
else handles.PhaseFactor = 1; end
guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function eCut2Dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eCut2Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function eCut2Dim_Callback(hObject, eventdata, handles)
% hObject    handle to eCut2Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eCut2Dim as text
%        str2double(get(hObject,'String')) returns contents of eCut2Dim as a double


% --- Executes on button press in cbExpression2Dim.
function cbExpression2Dim_Callback(hObject, eventdata, handles)
% hObject    handle to cbExpression2Dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbExpression2Dim


