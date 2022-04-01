function varargout = HFdataRead(varargin)
if nargin == 0  % LAUNCH GUI
    token = 'HFDATAREAD_open';
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
    hahdles.hh = 0;
    handles.mode = 3;
    handles.BLdata = [1, 1, 1, 1];
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
% End initialization code - DO NOT EDIT


% --- Executes just before HFdataRead is made visible.
function HFdataRead_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HFdataRead (see VARARGIN)

% Choose default command line output for HFdataRead
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HFdataRead wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HFdataRead_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function varargout = FileUpdate(varargin)

function eRe1_Callback(hObject, eventdata, handles)

function eRe2_Callback(hObject, eventdata, handles)

function eIm1_Callback(hObject, eventdata, handles)

function eIm2_Callback(hObject, eventdata, handles)

function eFld_Callback(hObject, eventdata, handles)

% phase correction objects
% ------------------------
function ePhReIm_Callback(hObject, eventdata, handles)
Plott(handles);
function ePh_Callback(hObject, eventdata, handles)
Plott(handles);
function ePh1_Callback(hObject, eventdata, handles)
Plott(handles);

function slPhReIm_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.ePhReIm);
set(handles.ePhReIm, 'String', num2str(val));
Plott(handles);

function slPhRe_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.ePh);
set(handles.ePh, 'String', num2str(val));
Plott(handles);

function slPhIm_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.ePh1);
set(handles.ePh1, 'String', num2str(val));
Plott(handles);
% ------------------------

function pbCalculate_Callback(hObject, eventdata, handles)
Plott(handles)

function Plott(ha)
h = guidata(ha.figure1);
hand = guidata(h.hh);
hand.out = {};

dim2 = length(safeget(hand.src.ax, 'z', 1));

dada = 1:dim2; % future possibility to select datasets

str = get(ha.eRe1,'String');
nRe1 = str2num(str);
str = get(ha.eRe2,'String');
nRe2 = str2num(str);
str = get(ha.eIm1,'String');
nIm1 = str2num(str);
str = get(ha.eIm2,'String');
nIm2 = str2num(str);
str = get(ha.eFld,'String');
nFld = str2num(str);

% nRe1d = nRe1*(dim2-1) + dada;
% nRe2d = nRe2*(dim2-1) + dada;
% nIm1d = nIm1*(dim2-1) + dada;
% nIm2d = nIm2*(dim2-1) + dada;
% 
nRe1d = nRe1;
nRe2d = nRe2;
nIm1d = nIm1;
nIm2d = nIm2;

str = get(ha.ePhReIm,'String');
PhReIm = str2num(str) * pi/180;
str = get(ha.ePh,'String');
Ph = str2num(str) * pi/180;
str = get(ha.ePh1,'String');
Ph1 = str2num(str) * pi/180;

Re1 = hand.src.y(:,nRe1d);
Re2 = hand.src.y(:,nRe2d);
Im1 = hand.src.y(:,nIm1d);
Im2 = hand.src.y(:,nIm2d);

if nFld
    nFldd = (nFld-1)*dim2 + dada;
    Fld = hand.src.y(:,nFldd);
else
    Fld(:, dada) = [0:(length(Re1)-1)]';
end

% if get(ha.chkLinFld, 'Value') |get(ha.chkSplFld, 'Value')
%     % sort Field axes points
%     [Fld, ind_F] = sort(Fld);
% else
%     ind_F = 1:length(Fld);
% end

if get(ha.chkSplFld, 'Value') %| dim2==1
    tc = str2num(get(ha.edTimeConst, 'String'));
    pol = str2num(get(ha.edPol, 'String'));    
    Fld = kv_mvgavg(Fld, tc, 'savgol', pol);
end
if dim2>1 % force linearisation of the field
    miF = max(min(Fld));
    maF = min(max(Fld));
    neField = linspace(miF, maF, length(Fld(:, 1))).';
    for ii = 1:dim2
%         [newx, Re1(:, ii)] = kv_spaceman(Fld(:, ii), Re1(:, ii), length(neField));
        [ttx, tx_ind] = sort(Fld(:, ii));
        Re1(:, ii) = interp1q(ttx, Re1(tx_ind, ii), neField);
%         [newx, Re2(:, ii)] = kv_spaceman(Fld(:, ii), Re2(:, ii), length(neField));
        Re2(:, ii) = interp1q(ttx, Re2(tx_ind, ii), neField);
%         [newx, Im1(:, ii)] = kv_spaceman(Fld(:, ii), Im1(:, ii), length(neField));
        Im1(:, ii) = interp1q(ttx, Im1(tx_ind, ii), neField);
%         [newx, Im2(:, ii)] = kv_spaceman(Fld(:, ii), Im2(:, ii), length(neField));
        Im2(:, ii) = interp1q(ttx, Im2(tx_ind, ii), neField);
    end
    Fld = neField;
end

hand.out{1} = [];
hand.out{1}.ax = hand.src.ax;
hand.out{1}.ax.x = Fld;
hand.out{1}.ax.xlabel = 'Magnetic Field';
hand.out{1}.ax.y = dada;
hand.out{1}.ax.ylabel = '';

Re  = exp(-i*Ph)*(Re1 + i*Re2);
Im  = exp(-i*Ph1)*(Im1 + i*Im2);
RRe = real(Re); RIm = real(Im);
ReIm  = exp(-i*PhReIm)*(RRe + i*RIm);



makeBL =get(ha.chkBL, 'Value');
if makeBL & (ha.mode<4) % if show field axes, do nothing
    nn = length(Re1);
    BLindex = [ha.BLdata(1):ha.BLdata(2), (nn - ha.BLdata(3) + 1):(nn - ha.BLdata(4) + 1)];
    polval = get(ha.popBLtype, 'Value')-1;
    [p]=polyfit(Fld(BLindex), Re1(BLindex),polval);
    bRe1 = polyval(p, Fld);
    [p]=polyfit(Fld(BLindex), Re2(BLindex),polval);
    bRe2 = polyval(p, Fld);
    [p]=polyfit(Fld(BLindex), Im1(BLindex),polval);
    bIm1 = polyval(p, Fld);
    [p]=polyfit(Fld(BLindex), Im2(BLindex),polval);
    bIm2 = polyval(p, Fld);

    bRe  = exp(-i*Ph)*(bRe1 + i*bRe2);
    bIm  = exp(-i*Ph1)*(bIm1 + i*bIm2);
    bRRe = real(bRe); bRIm = real(bIm);
    bReIm  = exp(-i*PhReIm)*(bRRe + i*bRIm);
    
    if ~get(ha.chkShowBL, 'Value')
        Re = Re - bRe;
        Im = Im - bIm;
        ReIm = ReIm - bReIm;
    end
end

switch(ha.mode)
    case 1, yy = Re;
    case 2, yy = Im;
    case 3, yy = ReIm;
    case 4, yy = Fld; hand.out{1}.ax.x = [1:length(Fld)]'; hand.our{1}.ax.xlabel ='Points';
end



if get(ha.chkLinFld, 'Value') & (dim2==1)
        [newx, yy] = kv_spaceman(Fld, yy, length(Fld));
        if ha.mode~=4,
            hand.out{1}.ax.x = newx;
        end
end
if dim2>1
    
end
if get(ha.chkInt, 'Value'), % make integration
    hand.out{1}.y = cumsum(yy, 1);
else
    hand.out{1}.y = yy;
end

if makeBL
    if get(ha.chkShowBL, 'Value') & (ha.mode<4)
        hand.out{2} = hand.out{1};
        switch(ha.mode)
            case 1, yy = bRe;
            case 2, yy = bIm;
            case 3, yy = bReIm;
        end
        
        dyy = yy([ha.BLdata(1), ha.BLdata(2), (nn - ha.BLdata(3) + 1), (nn - ha.BLdata(4) + 1)]);
        dxx = Fld([ha.BLdata(1), ha.BLdata(2), (nn - ha.BLdata(3) + 1), (nn - ha.BLdata(4) + 1)]);
        
        if get(ha.chkLinFld, 'Value')
            [newx, yy] = kv_spaceman(Fld, yy, length(Fld));
            if ha.mode~=4,
                hand.out{2}.ax.x = newx;
            end
        end
        
        if get(ha.chkInt, 'Value'), % make integration
            hand.out{2}.y = cumsum(yy);
        else
            hand.out{2}.y = yy;
        end
        
        hand.out{2}.ax.Color = [1, 0, 0];
        
        hand.out{3} = hand.out{2};
        hand.out{3}.y =dyy;
        hand.out{3}.ax.x = dxx;
        hand.out{3}.ax.LineStyle = 'none';
        hand.out{3}.ax.Marker = '*';
        hand.out{3}.ax.Color = [0, 1, 0];
        
    end
end
guidata(h.hh, hand);
[path,name,ext] = fileparts(get(h.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --- Executes on button press in tbRe.
function tbRe_Callback(hObject, eventdata, handles)
group = [handles.tbRe, handles.tbIm, handles.tbFinal, handles.tbField];
set(group, 'Value', 0);
set(hObject, 'Value', 1)

handles.mode = find(group == hObject);
guidata(handles.figure1, handles);
Plott(handles)


function chkSplFld_Callback(hObject, eventdata, handles)
if get(handles.chkSplFld, 'Value')
    str = 'on';
else
    str = 'off';
end
set(handles.edTimeConst, 'Enable', str);
set(handles.edPol, 'Enable', str);
set([handles.text14, handles.text15], 'Enable', str);

Plott(handles);


function chkLinFld_Callback(hObject, eventdata, handles)
Plott(handles);

function chkInt_Callback(hObject, eventdata, handles)
Plott(handles);

function edTimeConst_Callback(hObject, eventdata, handles)
Plott(handles);
function edPol_Callback(hObject, eventdata, handles)
Plott(handles);

function chkBL_Callback(hObject, eventdata, handles)
if get(handles.chkBL, 'Value')
    str = 'on';
else
    str = 'off';
end
set(handles.eBfirst, 'Enable', str);
set(handles.eEfirst, 'Enable', str);
set(handles.eBlast, 'Enable', str);
set(handles.eElast, 'Enable', str);
set(handles.slBfirst, 'Enable', str);
set(handles.slEfirst, 'Enable', str);
set(handles.slBlast, 'Enable', str);
set(handles.slElast, 'Enable', str);

set(handles.chkShowBL, 'Enable', str);
set(handles.popBLtype, 'Enable', str);
 Plott(handles);


function eBfirst_Callback(hObject, eventdata, handles)
try
    handles.BLdata(1) = str2num(get(hObject, 'String'));
    Plott(handles);
catch
    msgbox('Input parameter is not a number');
end


function slBfirst_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.eBfirst);
handles.BLdata(1) = val;
guidata(handles.figure1, handles);
Plott(handles);

function eBlast_Callback(hObject, eventdata, handles)
try
    handles.BLdata(2) = str2num(get(hObject, 'String'));
    Plott(handles);
catch
    msgbox('Input parameter is not a number');
end

function slBlast_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.eBlast);
handles.BLdata(2) = val;
guidata(handles.figure1, handles);
Plott(handles);

function eEfirst_Callback(hObject, eventdata, handles)
try
    handles.BLdata(3) = str2num(get(hObject, 'String'));
    Plott(handles);
catch
    msgbox('Input parameter is not a number');
end

function slEfirst_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.eEfirst);
handles.BLdata(3) = val;
guidata(handles.figure1, handles);
Plott(handles);

function eElast_Callback(hObject, eventdata, handles)
try
    handles.BLdata(4) = str2num(get(hObject, 'String'));
    Plott(handles);
catch
    msgbox('Input parameter is not a number');
end

function slElast_Callback(hObject, eventdata, handles)
val = sliderwork(hObject, handles.eElast);

handles.BLdata(4) = val;
guidata(handles.figure1, handles);
Plott(handles);

function chkShowBL_Callback(hObject, eventdata, handles)
Plott(handles);

function popBLtype_Callback(hObject, eventdata, handles)
Plott(handles);

function val = sliderwork(sliderhand, edithand)
shift =get(sliderhand, 'Value');
set(sliderhand, 'Value', 0.5);
try
    CurVal = str2num(get(edithand, 'String'));
    val = CurVal + floor((shift - 0.5)*2);
    if val <1, val=1; end
    set(edithand, 'String', num2str(val));

catch
    msgbox('wrong data');
end
