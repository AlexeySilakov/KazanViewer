function varargout = FTIRBox(varargin)
% COMPBOX Application M-file for CompBox.fig
%    FIG = COMPBOX launch CompBox GUI.
%    COMPBOX('callback_name', ...) invoke the named callback.

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 29-Oct-2008 16:41:15
% Alexey Silakov & Boris Epel 5-dec-2003, MPI
% boep 21-dec-2003, MPI

if nargin == 0  % LAUNCH GUI
    token = 'FTIRBox_open';
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
    handles.x = 0;
    handles.y = 0;
    handles.trace = [];
    handles.Color = [0, 0, 1;...
                     1, 0, 0;...
                     0, 1, 0;...
                     0, 0, 0;...
                     1, 0, 1;...
                     0, 1, 1;...
                     1, 1, 0];
    handles.wfshift = 0;
    hahdles.hh = 0;
    handles.ZoomUnits = 'Norm';
    set(handles.mPrAgr, 'Checked', 'on');
    handles.process = 1;
    handles.sumall = 0;
    handles.showall = 1;
    handles.change2D = 0;
    
    handles.PTSidx = [];
    handles.PTSx = [];
    handles.PTSdata = [];
    handles.currdata = 1;
    handles.Ref = [];
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
function varargout = FileUpdate(varargin)

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

shift =get(handles.Shift_X, 'Value');
set(handles.Shift_X, 'Value', 0.5);
num = get(handles.SXdim, 'Value');
str = get(handles.SXdim, 'String');
dim = str2num([str{num}]);
try
    CurVal = str2num(get(handles.edit1, 'String'));
    Val = CurVal + dim*(shift - 0.5)*2;
    set(handles.edit1, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ShiftX = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ZoomX_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

shift =get(handles.Zoom_X, 'Value');
set(handles.Zoom_X, 'Value', 0.5);
num = get(handles.ZXdim, 'Value');
str = get(handles.ZXdim, 'String');
dim = str2num([str{num}]);
try
    CurVal = str2num(get(handles.edit3, 'String'));
    Val = CurVal + dim*(shift - 0.5)*2;
    set(handles.edit3, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ZoomX = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ShiftY_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

shift =get(handles.Shift_Y, 'Value');
set(handles.Shift_Y, 'Value', 0.5);
num = get(handles.SYdim, 'Value');
str = get(handles.SYdim, 'String');
dim = str2num([str{num}]);
try
    CurVal = str2num(get(handles.edit2, 'String'));
    Val = CurVal + dim*(shift - 0.5)*2;
    set(handles.edit2, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ShiftY = Val;
    guidata(handles.figure1, handles);
    Plott(handles);
catch
    msgbox('wrong data');
end

% --------------------------------------------------------------------
function varargout = ZoomY_Callback(h, eventdata, handles, varargin)
if length(handles.trace) == 0, return; end

shift =get(handles.Zoom_Y, 'Value');
set(handles.Zoom_Y, 'Value', 0.5);
num = get(handles.ZYdim, 'Value');
str = get(handles.ZYdim, 'String');
dim = str2num([str{num}]);
try
    CurVal = str2num(get(handles.edit4, 'String'));
    Val = CurVal + dim*(shift - 0.5)*2;
    set(handles.edit4, 'String', num2str(Val));
    num = get(handles.popfiles, 'Value');
    handles.trace{num}.ZoomY = Val;
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
handles.trace{num}.ShiftX = str2num(get(handles.edit1, 'String'));
handles.trace{num}.ShiftY = str2num(get(handles.edit2, 'String'));
handles.trace{num}.ZoomX = str2num(get(handles.edit3, 'String'));
handles.trace{num}.ZoomY = str2num(get(handles.edit4, 'String'));
%handles.trace{num}.Data.ax.y = str2num(get(handles.eSecondDim, 'String'))';
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
file_str = [name, ext];
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

if isempty(ha.Ref),
    warning('No reference data loaded');
    set(ha.pbSetRef, 'BackgroundColor', [1 0.64 0.64]);
    return;
end
hand.out{1} = ha.Ref;
hand.out{1}.title = 'Ref';
hand.out{1}.ax.Color = [0, 0, 1];

yones = ones(size(ha.Ref.y, 2), 1);

if get(ha.chkUseBl, 'Value')
    if ~isfield(ha, 'BL'), warning('No baseline is loaded'); 
    elseif isempty(ha.BL), warning('No baseline is loaded');
    else
        diff1 = 0;
        if sum(size(ha.BL.ax.x)~= size(ha.Ref.ax.x)) 
            diff1 = 1;
        else
            if sum(ha.BL.ax.x - ha.Ref.ax.x), diff1 = 1; end
        end
        
        if diff1,
            BL = ha.BL;
            BL.y = spline(ha.BL.ax.x, ha.BL.y, ha.Ref.ax.x);
            BL.ax.x = ha.Ref.ax.x;
        else
            BL = ha.BL;
        end
        
        hand.out{end+1} = BL;
        hand.out{end}.title = 'BL';
        hand.out{end}.ax.Color = [1, 0, 0];
        BLsc = str2num(get(ha.edBLscale, 'String'));
        hand.out{end}.y = BL.y*BLsc;
        
        hand.out{end+1} = hand.out{1};
        hand.out{end}.y = -1*log10((ha.Ref.y)./(BL.y(:, yones)*BLsc));
        ha.Ref.y = hand.out{end}.y;
    end
end

if ha.showall
    for kk = 1:Num
        trace = h.trace{kk};
        if h.process, [trace.Data.y, trace.Data.ax] = processing(trace.Data.y, trace.Data.ax);
        else
            trace.Data.ax.filt = '';
            trace.Data.ax.diff = '';
        end
        
        hand.out{end+1}.ax = trace.Data.ax;
        hand.out{end}.ax.x = trace.Data.ax.x*trace.ZoomX + trace.ShiftX; % <-alsi 04.10.2005
        switch h.ZoomUnits
            case 'Norm'
                zoom = trace.ZoomY;
            case 'dB'
                zoom = 10^(trace.ZoomY/10);        
        end
        if safeget(trace, 'change2D', 0)
            sh = safeget(trace, 'ShiftY2D', 0)+trace.ShiftY;
            zm = safeget(trace, 'ZoomY2D', 1);
            hand.out{end}.y = (trace.Data.y*zoom*trace.YSign - sh).*zm;
        else
            hand.out{end}.y = trace.Data.y*zoom*trace.YSign +trace.ShiftY; % <-alsi 04.10.2005
        end
        
        hand.out{end}.ax.s = 0; % trace.ShiftY <-alsi 04.10.2005
        hand.out{end}.ax.dx= 0; % trace.ShiftX <-alsi 04.10.2005
        hand.out{end}.title = ['Trace', num2str(kk)];
    end
end

tmpy=[];
if get(ha.pbRefSub, 'Value')
    tmpy = ha.Ref.y;
    pol = str2num(get(ha.edSavgolPol, 'String'));
    tc = str2num(get(ha.edSavgolTC, 'String'));
    for ii = 1:size(tmpy, 2)
         ttt(:, ii) = cumsum(smooth([0; diff(tmpy(:, ii), 1, 1)],tc,'savgol',pol), 1);
    end
    ttt = [ttt-ttt(end, 1)]+ tmpy(end, 1);%
    tmpy = tmpy-ttt;
    
    hand.out{end+1} =  ha.Ref;
    hand.out{end}.y = ttt;
    hand.out{end}.ax.Color = [0, 1, 0];
    hand.out{end}.title = 'smooth ref';

end
    
if Num>0
    if isempty(tmpy)
        tmpy = ha.Ref.y;
    end
    hand.out{end+1} =  ha.Ref;
    for kk = 1:Num
        trace = h.trace{kk};
        if h.process,
            [trace.Data.y, trace.Data.ax] = processing(trace.Data.y, trace.Data.ax);
        else
            trace.Data.ax.filt = '';
            trace.Data.ax.diff = '';
        end        
        
%         hand.out{end+1}.ax = trace.Data.ax;
%         hand.out{end}.ax.x = trace.Data.ax.x*trace.ZoomX + trace.ShiftX; % <-alsi 04.10.2005
        switch h.ZoomUnits
            case 'Norm'
                zoom = trace.ZoomY;
            case 'dB'
                zoom = 10^(trace.ZoomY/10);        
        end
        diff1 = 0;
        if sum(size(trace.Data.ax.x)~= size(ha.Ref.ax.x)) 
            diff1 = 1;
        else
            if sum(trace.Data.ax.x - ha.Ref.ax.x), diff1 = 1; end
        end
        
        if diff1,
            trace.Data.y = spline(trace.Data.ax.x+trace.ShiftX, trace.Data.y, ha.Ref.ax.x);
            trace.Data.ax.x = ha.Ref.ax.x;
        end
            
        if safeget(trace, 'change2D', 0)
            sh = safeget(trace, 'ShiftY2D', 0)+trace.ShiftY;
            zm = safeget(trace, 'ZoomY2D', 1);
            tty = (trace.Data.y*zoom*trace.YSign - sh).*zm;
        else
            tty = trace.Data.y*zoom*trace.YSign +trace.ShiftY; % <-alsi 04.10.2005
        end
        
        tmpy = tmpy - tty(:, yones);
        
%         hand.out{end}.ax.s = 0;% trace.ShiftY <-alsi 04.10.2005
%         hand.out{end}.ax.dx= 0;% trace.ShiftX <-alsi 04.10.2005
%         hand.out{end}.title = ['Trace', num2str(kk)];
    end
%     hand.out{end}.y = tmpy;
%     hand.out{end}.ax.Color = [0, 0, 0];
%     hand.out{end}.title = 'Difference';
end
if ~isempty(tmpy);
    hand.out{end+1} =  ha.Ref;
    hand.out{end}.y = tmpy;
    hand.out{end}.ax.Color = [0, 0, 0];
    hand.out{end}.title = 'Difference';
end
if get(ha.chkFilterHF, 'Value')
    if isempty(tmpy);
        tmpy = ha.Ref.y;
    end
    ffty = (ifft(tmpy));
    fttt = ffty(:, 1)*0;
    lentty = floor(length(ffty)/2);
    x = [1:lentty]';
    x01 = str2num(get(ha.edFilterBeg, 'String'));
    x01 = max(x01, 1);
    f(:, 1) = -1./(1+exp((x-x01)/0.02));
    x02 = str2num(get(ha.edFilterEnd, 'String'));
    x02 = min(lentty, x02);
    f(:, 1) = 1-f(:, 1) - (1./(1+exp((x-x02)/0.02)));
    
    
%     yy = spline([x01; x02], real(ffty([x01, x02], 1)), [x01:x02]');
%     yy = yy + i*spline([x01; x02], imag(ffty([x01, x02], 1)), [x01:x02]');
    
    fttt(1:lentty) = f;
    fttt(end - [1:lentty]) = f(end:-1:1);

%     tmpy = fft(ffty.*f(:, yones));
    tmpy = real(fft(ffty.*fttt(:, yones)));

    hand.out{end+1} =  ha.Ref;
    hand.out{end}.y = tmpy;
    hand.out{end}.ax.Color = [0, 0, 1];
    hand.out{end}.title = 'Filtered';
end
tt =  hand.out{end};
if get(ha.chkDiff, 'Value')
    order = str2num(get(ha.edDiffOrder, 'String'));
    amp = str2num(get(ha.edDiffAmp, 'String'));

    for ii = 1:length(amp)
        diffs = Local_makediff(tt.y, order, amp(ii));

        tt.y = tt.y - diffs;  
    end
    hand.out{end+1} = hand.out{1};
    hand.out{end}.y = diffs;
    hand.out{end}.ax.Color = [0.5, 0.5, 0.5];
    hand.out{end}.title = 'dy';
    
    hand.out{end+1} = hand.out{1};
    hand.out{end}.y = tt.y;
    

end
tt =  hand.out{end};
if length(ha.PTSidx) > 0
    if get(ha.cbPointBased, 'Value')
        tty = ha.PTSdata(:, 2);
        ttx = ha.PTSdata(:, 1);
        avery = tty;
    else
        tty = tt.y(ha.PTSidx,:);
        ttx = tt.ax.x(ha.PTSidx);
        avery = tt.y(ha.PTSidx,:);
    end
    spl.y    = tty;
    spl.ax.x = ttx;
    spl.ax.Marker    = '*';
    spl.ax.Color     = 'r';
    spl.ax.LineStyle = 'none';
    spl.title = 'Marker';
    hand.out{end+1} = spl;
    
    hand.out{end+1} = spl;
    hand.out{end}.y = spl.y(ha.currdata);
    hand.out{end}.ax.x = spl.ax.x(ha.currdata);
    hand.out{end}.ax.Marker = 'o';
    hand.out{end}.ax.Color     = 'k';
%     avery = (tt.y(ha.PTSidx,:) + tt.y(ha.PTSidx-1,:) + tt.y(ha.PTSidx+1,:))/3;
    
    for ii=1:size(tty,2)
        warning off;
        spl1.y(:,ii) = tt.y(:, ii)*0;
        idxx = min(ha.PTSidx):max(ha.PTSidx);
        switch get(ha.pmMethod, 'Value')
        case 1
            Polnum = str2num(get(ha.ednPol, 'String'));
            [p]=polyfit(ttx, avery(:,ii),min(Polnum, length(ha.PTSidx)-1));
            spl1.y(idxx,ii) = polyval(p, tt.ax.x(idxx));
        case 2
            spl1.y(idxx,ii) = spline(ttx, ...
                avery(:,ii), tt.ax.x(idxx));
        end
        
    end
    spl1.ax.x = tt.ax.x;
    spl1.ax.y = tt.ax.y;
    spl1.ax.Color   = 'r';
    spl1.title       = 'Fit';
    hand.out{end+1} = spl1;
    
    spldif.ax.x = tt.ax.x;
    spldif.ax.y = tt.ax.y;
    spldif.ax.Color   = 'g';
    spldif.title       = 'FitDiff';
    spldif.ax.xlabel =     safeget(tt.ax, 'xlabel', '?, ?');
    spldif.ax.ylabel =     safeget(tt.ax, 'ylabel', '?, ?'); 
    if get(ha.cbPointBased, 'Value')
        spldif.y = tt.y - spl1.y(:, yones);
    else
        spldif.y = tt.y - spl1.y;
    end
    makeitcomplex = 1; % for a convenience 
    if makeitcomplex
        spldif.y = real(spldif.y) + i*real(spldif.y);
    end
    hand.out{end+1} = spldif;
%     spldif.y = cumsum(spldif.y);
%     spldif.ax.Color   = 'k';
%     spldif.title       = 'FitDiff_int';
%     
%     hand.out{end+1} = spldif;

    
end

guidata(h.hh, hand);
[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function UpdateInterface(handles)
if length(handles.trace) == 0
    set(handles.edit1, 'String', '0');
    set(handles.edit2, 'String', '0');    
    set(handles.edit3, 'String', '1');    
    set(handles.edit4, 'String', '1');    
    set(handles.text6, 'BackgroundColor', [0.75,0.75,0.75]);
    %set(handles.eSecondDim, 'String', '-');    
else
    num = get(handles.popfiles, 'Value');
    set(handles.edit1, 'String', num2str(handles.trace{num}.ShiftX));
    set(handles.edit2, 'String', num2str(handles.trace{num}.ShiftY));
    set(handles.edit3, 'String', num2str(handles.trace{num}.ZoomX));
    set(handles.edit4, 'String', num2str(handles.trace{num}.ZoomY));
    set(handles.text6, 'BackgroundColor', handles.trace{num}.Data.ax.Color);
    str = num2str(handles.trace{num}.Data.ax.y');
   % set(handles.eSecondDim, 'String', str);    
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

% --- Executes on selection change in SYdim.
function SYdim_Callback(hObject, eventdata, handles)

% --- Executes on selection change in ZXdim.
function ZXdim_Callback(hObject, eventdata, handles)

% --- Executes on selection change in ZYdim.
function ZYdim_Callback(hObject, eventdata, handles)

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
handles.showall = ~stat;

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
for ii=1:counter
   handles.trace{end+1}.Data.ax = trace.Data.ax;
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
        'callback', 'ftirbox(''mDisChanges_Callback'', guidata(gcbo))');
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


% --- Executes on button press in pbSetRef.
function pbSetRef_Callback(h, eventdata, handles)
% From MinuxBox
hand = guidata(handles.hh);
handles.Ref = hand.src;
[path,name,ext] = fileparts(safeget(handles.Ref, 'filename', 'noname.xxx'));
% n = size(hand.src.y, 2);
% str{1} = 'All';
% for ii=1:n, str{ii+1}=num2str(ii); end
% % set(handles.pmSlice, 'String', str, 'Value', 1);
% while length(handles.set) < n, handles.set(end+1) = handles.set(1); end    
set(h, 'BackgroundColor', [0.64 0.8 0.64]);
set(h, 'String', ['Ref: ', name, ext])
guidata(handles.figure1, handles);
% pmSlice_Callback(h, eventdata, handles);
Plott(handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, ha)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 88888)']);
if isempty(x); return; end

if ~isempty(ha.PTSdata)
    ha.PTSdata = [ha.PTSdata; [x, y]];
else
    ha.PTSdata = [x, y];
end
[ha.PTSdata(:, 1), idx] = sort(ha.PTSdata(:, 1));
ha.PTSdata(:, 2) = ha.PTSdata(idx, 2);

hand = guidata(ha.hh);
axx = ha.Ref.ax.x;

[mm, idx] = min(abs(axx(:, ones(length(x), 1)) - x(:, ones(length(axx), 1))'), [], 1);
ha.PTSidx = [ha.PTSidx, idx];
ha.PTSx = [ha.PTSx; axx(idx)];
% for ii=1:length(x)
%     [mm,idx] = min(abs(axx-x(ii)));
%     ha.PTSidx(end+1) =idx;
%     ha.PTSx(end+1)   =axx(idx);
% end

[ha.PTSx, idx] = sort(ha.PTSx);
ha.PTSidx      = ha.PTSidx(idx);

guidata(ha.figure1, ha);

PointsList1(ha);
Plott(ha);

% --- Executes on button press in tbAvSubtracrt.
function tbAvSubtracrt_Callback(hObject, eventdata, handles)
% hObject    handle to tbAvSubtracrt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plott(handles);

function eAverage_Callback(hObject, eventdata, handles)
% hObject    handle to eAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plott(handles)

% --- Executes on selection change in lbPoints.
function lbPoints_Callback(hObject, eventdata, handles)
num = get(hObject, 'Value');
    if get(handles.cbPointBased, 'Value')
        tty = handles.PTSdata(num, 2);
        ttx = handles.PTSdata(num, 1);
        str = sprintf('%10.6f', tty);    
        set(handles.edPtsYaxis, 'String', str);
    else
        ttx = handles.PTSidx(num);
    end
str = sprintf('%10.6f', ttx);    
set(handles.edPtsXaxis, 'String', ttx);
handles.currdata = num;
guidata(handles.figure1, handles);
Plott(handles)
% --- Executes on button press in pbClear.
function pbClear_Callback(hObject, eventdata, handles)
% hObject    handle to pbClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PTSidx = [];
handles.PTSx = [];
handles.PTSdata = [];
handles.currdata = 1;
guidata(handles.figure1, handles);
PointsList1(handles)
Plott(handles)

% --- Executes on button press in cbPointBased.
function cbPointBased_Callback(hObject, eventdata, handles)
Plott(handles)
% hObject    handle to cbPointBased (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on selection change in pmMethod.
function pmMethod_Callback(hObject, eventdata, ha)
Plott(ha)
% hObject    handle to pmMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pmMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmMethod


% --- Executes on button press in pbAutoPeaks.
function pbAutoPeaks_Callback(hObject, eventdata, ha)
% hObject    handle to pbAutoPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hand = guidata(ha.hh);
axx = ha.Ref.ax.x;

tt = ha.Ref;

% if get(ha.tbAvSubtracrt, 'Value')
%     tt.ax.filt     = 'savgol';
%     tt.ax.filtpar1 = '1D';
%     tt.ax.pol      = 2;
%     tt.ax.tc       = str2num(get(ha.eAverage, 'String'));
%     [tt1.y, tt.ax] = processing(tt.y, tt.ax);
%     tt.y = tt.y - tt1.y;
% end

tt.ax.filt     = 'savgol';
tt.ax.filtpar1 = '1D';
tt.ax.pol      = 2;
tt.ax.tc       = str2num(get(ha.eTolerance, 'String'));
[tt.y, tt.ax] = processing(tt.y, tt.ax);

idx = findpeaks(-tt.y);
for ii=1:length(idx)
    ha.PTSidx(end+1) =idx(ii);
    ha.PTSx(end+1)   =axx(idx(ii));
end

% [yyy, idx] = sort(ha.PTSx);
% idx = idx( idx > 3 & idx < length(axx)-2);
% 
ha.PTSidx      = ha.PTSidx(2:end-1);
ha.PTSx        = ha.PTSx(2:end-1);

guidata(ha.figure1, ha);

PointsList1(ha);
Plott(ha)


function eTolerance_Callback(hObject, eventdata, handles)
% hObject    handle to eTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eTolerance as text
%        str2double(get(hObject,'String')) returns contents of eTolerance as a double


% --- Executes on button press in pbChange1.
function pbChange1_Callback(hObject, eventdata, ha)
% hObject    handle to pbChange1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 2)']);

if length(x)~=2, return; end;

hand = guidata(ha.hh);
axx = ha.Ref.ax.x;

% [mm,pidx] = min(abs(ha.PTSx-x(1)));
[mm,pidx] = min(abs(ha.PTSdata(:, 1)-x(1)));
pidx = pidx(1);

[mm,idx] = min(abs(axx-x(2)));
ha.PTSidx(pidx) =idx;
ha.PTSx(pidx)   =axx(idx);
ha.PTSdata(pidx, :) = [x(2), y(2)];

[ha.PTSdata(:, 1), idx] = sort(ha.PTSdata(:, 1));
ha.PTSdata(:, 2) = ha.PTSdata(idx, 2);

[ha.PTSx, idx] = sort(ha.PTSx);
ha.PTSidx      = ha.PTSidx(idx);

guidata(ha.figure1, ha);

PointsList1(ha);
Plott(ha)


function PointsList1(handles)
Str = {};

for ii=1:length(handles.PTSidx)
    Str{end+1} = sprintf('[%5d] %5.3f %5.3f',handles.PTSidx(ii), handles.PTSdata(ii, 1), handles.PTSdata(ii, 2));
end

set(handles.lbPoints,'String',Str, 'Visible', 'on');

function [ind,peaks] = findpeaks(y)
% FINDPEAKS  Find peaks in real vector.
%  ind = findpeaks(y) finds the indices (ind) which are
%  local maxima in the sequence y.  
%
%  [ind,peaks] = findpeaks(y) returns the value of the peaks at 
%  these locations, i.e. peaks=y(ind);

y = y(:)';

switch length(y)
case 0
    ind = [];
case 1
    ind = 1;
otherwise
    dy = diff(y);
    not_plateau_ind = find(dy~=0);
    ind = find( ([dy(not_plateau_ind) 0]<0) & ([0 dy(not_plateau_ind)]>0) );
    ind = not_plateau_ind(ind);
    if y(1)>y(2)
        ind = [1 ind];
    end
    if y(end-1)<y(end)
        ind = [ind length(y)];
    end
end

if nargout > 1
    peaks = y(ind);
end


% --- Executes on button press in pbRefSub.
function pbRefSub_Callback(hObject, eventdata, handles)
Plott(handles);


function edSavgolTC_Callback(hObject, eventdata, handles)
Plott(handles);

function edSavgolPol_Callback(hObject, eventdata, handles)
Plott(handles);

function pbDelete1_Callback(hObject, eventdata, ha)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 1)']);

if length(x)~=1, return; end;

hand = guidata(ha.hh);
axx = ha.Ref.ax.x;

% [mm,pidx] = min(abs(ha.PTSx-x(1)));
[mm,pidx] = min(abs(ha.PTSdata(:, 1)-x(1)));
pidx = pidx(1);


% deletion of the element pidx
ha.PTSidx(pidx) = []; 
ha.PTSx(pidx)   =[];

ha.PTSdata(pidx, :) = [];
% [ha.PTSx, idx] = sort(ha.PTSx);
% ha.PTSidx      = ha.PTSidx(idx);

guidata(ha.figure1, ha);

PointsList1(ha);
Plott(ha)


% --- Executes on button press in pbLoadBl.
function pbLoadBl_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
handles.BL = hand.src;
if min(hand.src.y) <=0
    handles.BL.y = hand.src.y - min(hand.src.y) + 1e-15; % to prevent "Divide by zero" crap
end
[path,name,ext] = fileparts(safeget(handles.BL, 'filename', 'noname.xxx'));
% n = size(hand.src.y, 2);
% str{1} = 'All';
% for ii=1:n, str{ii+1}=num2str(ii); end
% % set(handles.pmSlice, 'String', str, 'Value', 1);
% while length(handles.set) < n, handles.set(end+1) = handles.set(1); end    
set(h, 'BackgroundColor', [0.64 0.8 0.64]);
set(h, 'String', ['BL: ', name, ext])
guidata(handles.figure1, handles);
% pmSlice_Callback(h, eventdata, handles);
Plott(handles);

% --- Executes on button press in chkUseBl.
function chkUseBl_Callback(hObject, eventdata, handles)
Plott(handles);


function ednPol_Callback(hObject, eventdata, handles)
% hObject    handle to ednPol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ednPol as text
%        str2double(get(hObject,'String')) returns contents of ednPol as a double
Plott(handles);

% --- Executes on button press in chkFilterHF.
function chkFilterHF_Callback(hObject, eventdata, handles)
% hObject    handle to chkFilterHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plott(handles);

% Hint: get(hObject,'Value') returns toggle state of chkFilterHF


function edFilterBeg_Callback(hObject, eventdata, handles)
% hObject    handle to edFilterBeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plott(handles);

% Hints: get(hObject,'String') returns contents of edFilterBeg as text
%        str2double(get(hObject,'String')) returns contents of edFilterBeg as a double

function edFilterEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edFilterEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edFilterEnd as text
%        str2double(get(hObject,'String')) returns contents of edFilterEnd as a double
Plott(handles);


% --- Executes on button press in pbShowFFT.
function pbShowFFT_Callback(hObject, eventdata, ha)
tmpy = ha.Ref.y;
ffty = ifft(tmpy);
    lentty = floor(length(ffty)/2);
    x = [1:lentty]';
    f = x*0;
    x01 = str2num(get(ha.edFilterBeg, 'String'));
    x01 = max(x01, 1);
    f(:, 1) = -1./(1+exp((x-x01)/0.2));
    x02 = str2num(get(ha.edFilterEnd, 'String'));
    x02 = min(lentty, x02);
    f(:, 1) =  1-f(:, 1) - (1./(1+exp((x-x02)/0.2)));
%     
%     yy = spline([x01; x02], real(ffty([x01, x02], 1)), [x01:x02]');
%     yy = yy + i*spline([x01; x02], imag(ffty([x01, x02], 1)), [x01:x02]');
    
% x = [1:(length(ffty)/2)]';
% x0 = str2num(get(ha.edFilterBeg, 'String'));
% x0 = max(x0, 1);
% % f(:, 1) = 1-1./(1+exp((x-x0)/2));
% x0 = str2num(get(ha.edFilterEnd, 'String'));
% x0 = min(length(ffty)/2, x0);
% f(:, 1) = f(:, 1) + (1./(1+exp((x-x0)/2))-1);
% % tmpy = fft(ffty.*f(:, yones));

figure(234345);
plot(x, abs(ffty(1:end/2))/max(abs(ffty)), 'b', x, f, 'r');
title('fft(ref)');
xlabel('points');
legend('Data', 'Filter');

function edBLscale_Callback(hObject, eventdata, handles)
Plott(handles);

function edPtsYaxis_Callback(hObject, eventdata, handles)
% datas = str2num(get(hObject, 'String'));
datas = sscanf(get(hObject, 'String'), '%10f');

if get(handles.cbPointBased, 'Value') % real data
    handles.PTSdata(handles.currdata, 2) = datas;
    
    guidata(handles.figure1, handles);
Plott(handles);
PointsList1(handles);
end


function edPtsXaxis_Callback(hObject, eventdata, handles)
datas = str2num(get(hObject, 'String'));
axx = handles.Ref.ax.x;

if get(handles.cbPointBased, 'Value') % real data
    handles.PTSdata(handles.currdata, 1) = datas;
    
    [mm,cidx] = min(abs(axx - datas));
    handles.PTSidx(handles.currdata) =cidx;
    handles.PTSx(handles.currdata)   =axx(cidx);
else % based on a poin
    cidx = min(max(floor(datas), 1), size(axx, 1));
    handles.PTSidx(handles.currdata) =cidx;
    handles.PTSx(handles.currdata)   =axx(cidx);
    handles.PTSdata(handles.currdata, 1) = axx(cidx);
end

[handles.PTSdata(:, 1), idx] = sort(handles.PTSdata(:, 1));

handles.PTSdata(:, 2) = handles.PTSdata(idx, 2);

[handles.PTSx, idx] = sort(handles.PTSx);
handles.PTSidx      = handles.PTSidx(idx);

ixz = find( handles.PTSidx == cidx);
handles.currdata = ixz(1);
guidata(handles.figure1, handles);
PointsList1(handles);
set(handles.lbPoints, 'Value', ixz(1));
Plott(handles);

function slPtsXaxis_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
shift = get(handles.slPtsXaxis, 'Value');
set(handles.slPtsXaxis, 'Value', 0.5);
num = get(handles.slPtsXaxis, 'UserData');
dim = num;

CurVal = str2num(get(handles.edPtsXaxis, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.edPtsXaxis, 'String', num2str(Val));
edPtsXaxis_Callback(handles.edPtsXaxis, eventdata, handles);


function slPtsYaxis_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
shift = get(handles.slPtsYaxis, 'Value');
set(handles.slPtsYaxis, 'Value', 0.5);
num = get(handles.slPtsYaxis, 'UserData');
dim = num;

CurVal = str2num([get(handles.edPtsYaxis, 'String'), '*1e7']);

Val = CurVal + dim*1e7*(shift - 0.5)*2;
str = sprintf('%10.6f', Val*1e-7);
set(handles.edPtsYaxis, 'String', str);
edPtsYaxis_Callback(handles.edPtsYaxis, eventdata, handles);

function SliderMenu1(hObject, eventdata, handles)
% hObject    handle to em6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hhaa = get(handles.Slider1, 'Children');
set(hhaa, 'Checked', 'off');
num = str2num(get(hObject, 'Label'));
set(handles.slPtsYaxis, 'UserData', num);
set(hObject, 'Checked', 'on');

function SliderMenu2(hObject, eventdata, handles)
hhaa = get(handles.Slider2, 'Children');
set(hhaa, 'Checked', 'off');
num = str2num(get(hObject, 'Label'));
set(handles.slPtsXaxis, 'UserData', num);
set(hObject, 'Checked', 'on');


% --- Executes on button press in chkDiff.
function chkDiff_Callback(hObject, eventdata, handles)
Plott(handles);

function edDiffOrder_Callback(hObject, eventdata, handles)
Plott(handles);

% --- Executes during object creation, after setting all properties.
function edDiffOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edDiffAmp_Callback(hObject, eventdata, handles)
Plott(handles);

% --- Executes during object creation, after setting all properties.
function edDiffAmp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbShowDiff.
function pbShowDiff_Callback(hObject, eventdata, ha)
% Plott(handles);

order = str2num(get(ha.edDiffOrder, 'String'));
amp = str2num(get(ha.edDiffAmp, 'String'));

if get(ha.chkUseBl, 'Value')
    len = size(ha.Ref.y, 2);
    if ~isfield(ha, 'BL'), warning('No baseline is loaded'); 
    elseif isempty(ha.BL), warning('No baseline is loaded');
    else
        tmpy = -1*log10((ha.Ref.y)./(ha.BL.y(:, ones(len, 1))));
    end
else
    tmpy = ha.Ref.y;
end
diff = Local_makediff(tmpy, order, amp);

figure(65332); plot(ha.Ref.ax.x, diff);


function diffs = Local_makediff(y, order, amp)
tty = y;
tty = kv_mvgavg(tty, 2, 'savgol', 1);
for ii = 1:order
    tty = diff(tty, 1);
end
len = size(tty, 2);
diffs = [zeros(ceil(order/2), len); tty; zeros(floor(order/2), len)]*amp;

diffs = kv_mvgavg(diffs, 2, 'savgol', 1);

function amps = Local_autodiff(y, order, amps)
global Data norder

MaxIterations = 500;
nunits = 100;
norder = [];
norder = order;
Data = real(y);
% 
% bounds(:, 1) = -abs(amps);
% bounds(:, 2) = +abs(amps);
% 
%  initPop=initializega(nunits,bounds, 'FTIRBox(''local_lsqlsq1'')');
%         [sol,endPop,bPop,traceInfo]=ga(bounds,'FTIRBox(''local_lsqlsq1'')', [],initPop,[1e-6 1 0],'maxGenTerm',MaxIterations);
% amps = sol(1:end-1) ; 
sol = fminsearch(@local_lsqlsq1, amps)
        
amps = sol;

        
function  val = local_lsqlsq1(sol);
global Data norder

tData = Data;
for ii = 1:length(sol)
    diffs = Local_makediff(tData, norder, sol(ii));
    tData = tData - diffs;
end
bly = kv_mvgavg(Data - tData, 100, 'savgol', 1);
y = tData + bly;
%Res = [h(:, 1)+h(:, 2)+0.1,  h(:, 3), h(:, 4)];
val = 1*sum(abs(Data(10:end) - y(10:end)));


% --------------------------------------------------------------------
function AutoDiffAmps_Callback(hObject, eventdata, ha)

hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 2)']);
if length(x)<2; return; end

hand = guidata(ha.hh);
axx = ha.Ref.ax.x;

[mm, idx] = min(abs(axx(:, ones(length(x), 1)) - x(:, ones(length(axx), 1))'), [], 1);

% for ii=1:length(x)
%     [mm,idx] = min(abs(axx-x(ii)));
%     ha.PTSidx(end+1) =idx;
%     ha.PTSx(end+1)   =axx(idx);
% end

tt = hand.out{end};

idx = sort(idx);
y = tt.y(idx(1):idx(2));

order = str2num(get(ha.edDiffOrder, 'String'));
amp = str2num(get(ha.edDiffAmp, 'String'));

resamps = Local_autodiff(y, order, amp);

set(ha.edDiffAmp, 'String', num2str(resamps));


function edNPFFT_Callback(hObject, eventdata, handles)
if get(handles.chkUseFFT, 'Value')
    Local_makefft(handles);
end

function chkUseFFT_Callback(hObject, eventdata, handles)
 handles = guidata(handles.figure1);
 if get(hObject, 'Value')
     % to backup
     if ~isfield(handles.Ref, 'ty')
         handles.Ref.ty = handles.Ref.y;
         handles.Ref.ax.tx = handles.Ref.ax.x;
     end
     if isfield(handles, 'BL');
         if ~isfield(handles.BL, 'ty')
             handles.BL.ty = handles.BL.y;
             handles.BL.ax.tx = handles.BL.ax.x;
         end
     end
     guidata(handles.figure1, handles);
     Local_makefft(handles)
 else
     if isfield(handles.Ref, 'ty')
         handles.Ref.y = handles.Ref.ty;
         handles.Ref.ax.x = handles.Ref.ax.tx;
     end
     if isfield(handles, 'BL');
         if ~isfield(handles.BL, 'ty')
             handles.BL.y = handles.BL.ty;
             handles.BL.ax.x = handles.BL.ax.tx;
         end
     end
     guidata(handles.figure1, handles);
     Plott(handles);
end

function Local_makefft(handles)
handles = guidata(handles.figure1);

Npoints = str2num(get(handles.edNPFFT, 'String'));
% to backup
% if ~isfield(handles.Ref, 'ty')
%     handles.Ref.ty = handles.Ref.y;
%     handles.Ref.ax.tx = handles.Ref.ax.x;
% end

[dataX, dataY] = Local_fftstuff(handles.Ref.ty, handles.Ref.dsc, Npoints);

handles.Ref.y = dataY;
handles.Ref.ax.x = dataX;

if get(handles.chkUseBl, 'Value')
    if ~isfield(handles, 'BL'), warning('No baseline is loaded'); 
    elseif isempty(handles.BL), warning('No baseline is loaded');
    else
        if ~isfield(handles.BL, 'ty')
            handles.BL.ty = handles.BL.y;
            handles.BL.ax.tx = handles.BL.ax.x;
        end       
        [blX, blY] = Local_fftstuff(handles.BL.ty, handles.BL.dsc, Npoints);
        % to backup
       
        handles.BL.y= blY;
        handles.BL.ax.x =blX;

    end
end

guidata(handles.figure1, handles);

Plott(handles)

function [x, y] = Local_fftstuff(iny, dsc1, Npoints)

iny = iny - mean([iny(1:10); iny((-10:-1:0)+end)]);

BPeakLoc = sscanf(dsc1.BackwardPeakLocation, '%f');
HFlim = sscanf(dsc1.HighFoldingLimit, '%f cm-1');
if Npoints < 10 % it is not "NPoints", but a resolution in [cm-1]
    Npoints = 2*floor(HFlim/Npoints);
end

if BPeakLoc %%% Backward and forward
    
    tmpy1 = iny(1:(end/2), :);
    tmpy2 = -iny((end/2+1):end, :);
    

    [xx, mind2] = max(tmpy2);
    [xx, miind2] = min(tmpy2);
    BPeakLoc = floor(mean([mind2, miind2]));

    % tmpy1 = tmpy1/max(tmpy1);
    % tmpy2 = tmpy2/max(tmpy2);

    [xx, mind1] = max(tmpy1);
    [xx, miind1] = min(tmpy1);
    PeakLoc = floor(mean([mind1, miind1]));
    if Npoints >= length(tmpy1)/2
        Npoints = PeakLoc;
        warning('Npoints > lenght(y)/2....Npoints = length(y)/2-1');
    end

    nzeros = 1;
    NewP = 2^(ceil(log(Npoints)/log(2))+nzeros);
    Zersss = zeros(NewP-Npoints, size(tmpy1, 2));

    % Blackman-Harris 3-term window
    tx = [0:(Npoints*2)].';
    Appod = 0.42323 - 0.49755*cos(2*pi*tx/(Npoints*2)) + 0.07922*cos(4*pi*tx/(Npoints*2));
    data1 = [Zersss; Appod(:, ones(size(tmpy1, 2), 1)).*tmpy1([-Npoints:Npoints]+PeakLoc+1, :); Zersss];
    data2 = [Zersss; Appod(:, ones(size(tmpy1, 2), 1)).*tmpy2([-Npoints:Npoints]+BPeakLoc+1, :); Zersss];

    % data = data + [Zersss; Appod.*tmpy2([-Npoints:Npoints]+BPeakLoc+1); Zersss];

    % figure; plot([data1, data2]);

    y = abs(fft(data1, [], 1))+ abs(fft(data2, [], 1));

else
    tmpy1 = iny;
    [xx, mind1] = max(tmpy1);
    [xx, miind1] = min(tmpy1);
    

    PeakLoc = floor(mean([mind1, miind1]));
    if Npoints >= length(tmpy1)/2
        Npoints = PeakLoc;
        disp('Npoints > lenght(y)/2....Npoints = length(y)/2-1');
    end
    nzeros = 1;
    NewP = 2^(ceil(log(Npoints)/log(2))+nzeros);
    Zersss = zeros(NewP-Npoints, size(tmpy1, 2));

    % Blackman-Harris 3-term window
    tx = [0:(Npoints*2)].';
    Appod = 0.42323 - 0.49755*cos(2*pi*tx/(Npoints*2)) + 0.07922*cos(4*pi*tx/(Npoints*2));
    data1 = [Zersss; Appod(:, ones(size(tmpy1, 2), 1)).*tmpy1([-Npoints:Npoints]+PeakLoc+1, :); Zersss];
    
    y = abs(fft(data1, [], 1));
end

x = linspace(0, HFlim*2, length(y)).';
y = y(1:(end/2), :);
x = x(1:(end/2));

function Local_step(sizey, x01, x02)
    lentty = floor(sizey/2);
    x = [1:lentty]';
    
    x01 = str2num(get(ha.edFilterBeg, 'String'));
    x01 = max(x01, 1);
    f(:, 1) = -1./(1+exp((x-x01)/0.02));
    x02 = str2num(get(ha.edFilterEnd, 'String'));
    x02 = min(lentty, x02);
    f(:, 1) = 1-f(:, 1) - (1./(1+exp((x-x02)/0.02)));
