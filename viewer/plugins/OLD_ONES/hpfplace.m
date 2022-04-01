function varargout = hpfplace(varargin)
% HPFPLACE Application M-file for hpfplace.fig
%    FIG = HPFPLACE launch hpfplace GUI.
%    HPFPLACE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 23-Feb-2005 12:40:14

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.HPFx = [];
    handles.HPFm = [];
    handles.HPFc = [];
    handles.HPFlw = [];
    [pathstr, a, b] = fileparts(get(fig, 'FileName'));
    outstr = inimanaging([pathstr, '\kazan.ini']);
    try handles.currdir = outstr.Common_For_Plugins.currDir; catch handles.currdir = pwd; end
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
function FileUpdate(handles)



% --------------------------------------------------------------------
function varargout = lbList_Callback(h, eventdata, ha, varargin)

pos = get(ha.lbList, 'Value');
set(ha.ePos, 'String', num2str(ha.HPFx(pos)));
set(ha.eMagn, 'String', num2str(ha.HPFm(pos)));
set(ha.eCoupling, 'String', num2str(ha.HPFc(pos)));
set(ha.eLW, 'String', num2str(ha.HPFlw(pos)));


% --------------------------------------------------------------------
function varargout = pbAdd_Callback(h, eventdata, ha, varargin)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 88888)']);

hand = guidata(ha.hh);
axx = hand.src.ax.x;

for ii=1:length(x)
    ha.HPFx(end+1) = x(ii);
    ha.HPFm(end+1) = 1;
    ha.HPFc(end+1) = 42;
    ha.HPFlw(end+1) = 34;
end

[ha.HPFx, idx] = sort(ha.HPFx);
ha.HPFm        = ha.HPFm(idx);
ha.HPFc        = ha.HPFc(idx);
ha.HPFlw       = ha.HPFlw(idx);
guidata(ha.figure1, ha);

HPFList(ha);
Plot(ha)
lbList_Callback(h, [], ha)


% --------------------------------------------------------------------
function varargout = pbRemove_Callback(h, eventdata, ha, varargin)
pos = get(ha.lbList,'Val');

idx = 1:length(ha.HPFx);
idx = idx~=pos;

ha.HPFx=ha.HPFx(idx);
ha.HPFm=ha.HPFm(idx);
ha.HPFc=ha.HPFc(idx);
ha.HPFlw=ha.HPFlw(idx);
guidata(ha.figure1, ha);

set(ha.lbList,'Val',min(pos,length(ha.HPFx)));

HPFList(ha);
Plot(ha)
lbList_Callback(h, [], ha)
% --------------------------------------------------------------------
function varargout = ePos_Callback(h, eventdata, handles, varargin)
pbChange_Callback(h, [], handles)

% --------------------------------------------------------------------
function varargout = eSpin_Callback(h, eventdata, handles, varargin)
pbChange_Callback(h, [], handles)

% --------------------------------------------------------------------
function varargout = eLW_Callback(h, eventdata, handles, varargin)
pbChange_Callback(h, [], handles)

% --------------------------------------------------------------------
function varargout = eMagn_Callback(h, eventdata, handles, varargin)
pbChange_Callback(h, [], handles)

% --------------------------------------------------------------------
function varargout = eCoupling_Callback(h, eventdata, handles, varargin)
pbChange_Callback(h, [], handles)

% --------------------------------------------------------------------
function varargout = pbChange_Callback(h, eventdata, ha, varargin)
pos = get(ha.lbList, 'Value');
ha.HPFx(pos) = str2num(get(ha.ePos, 'String'));
ha.HPFm(pos) = str2num(get(ha.eMagn, 'String'));
ha.HPFc(pos) = str2num(get(ha.eCoupling, 'String'));
ha.HPFlw(pos) = str2num(get(ha.eLW, 'String'));

guidata(ha.figure1, ha);
HPFList(ha);
Plot(ha)
% --------------------------------------------------------------------

function HPFList(ha)
Str = {};

for ii=1:length(ha.HPFx)
    Str{end+1} = sprintf('x=%5.3f  c=%4.2f m=%4.2f lw=%4.1f',ha.HPFx(ii), ha.HPFc(ii), ha.HPFm(ii), ha.HPFlw(ii));
end
val=get(ha.lbList,'Value')
set(ha.lbList,'String',Str,'Value',max([min([val,length(Str)]),1]));
% --------------------------------------------------------------------

function Plot(ha)
handles = guidata(ha.figure1);

hand = guidata(handles.hh);
hand.out = {};
hand.out{1} = hand.src;
x = hand.src.ax.x;
[sz1,sz2]=size(hand.src.y);

lc = eval(['[',get(handles.eSpin, 'String'),']']);
lp = -(length(lc)-1)/2:(length(lc)-1)/2;

res = zeros(length(x),1);
lw = str2num(get(handles.eLW, 'String'));
for ii=1:length(ha.HPFx)
    ml = zeros(length(x),1);
    for ptn = 1:length(lc)
        lsh = lshape(x,ha.HPFx(ii) + lp(ptn)*ha.HPFc(ii), ha.HPFlw(ii), 1);
        ml  = ml + lsh*lc(ptn);
    end
    ml  = ml./max(ml);
    res = res + ml * ha.HPFm(ii);
end

sh = str2num(get(handles.eXshift, 'String'));

hpf.y         = res(:,ones(1,sz2));
hpf.ax.x      = x;
hpf.ax.dx     = sh;
hpf.ax.xlabel = 'Magnetic Field, ';
hpf.ax.Color     = 'r';
hpf.ax.y      = 0;
hpf.title = 'hpf';
hand.out{2} = hpf;

[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function varargout = pbCopy_Callback(h, eventdata, ha, varargin)
Str = get(ha.lbList,'String');

out = '';
for kk = 1:length(Str)
    out = sprintf( '%s%s\n', out, Str{kk});
end
clipboard('copy',out);
disp(['hpfplace: list copied to the clipboard.']);


% --------------------------------------------------------------------
function varargout = eXshift_Callback(h, eventdata, ha, varargin)
Plot(ha)

% --------------------------------------------------------------------
function mFileLoad_Callback(hObject, eventdata, ha)
[fname,pname] = uigetfile([ha.currdir, '\*.m'],'Open List');
if fname == 0, return; end
ha.currdir = pname;
filename = [pname '\' fname];
fid = fopen(filename, 'r');
if fid == -1, disp('File not found'); return; end

ha.HPFx = [];
ha.HPFm = [];
ha.HPFc = [];
ha.HPFlw = [];

while feof(fid) < 1,
    s = fscanf(fid,'%f,%f,%f,%f\n',4);
    ha.HPFx(end+1)=s(1); ha.HPFm(end+1)=s(2);
    ha.HPFc(end+1)=s(3); ha.HPFlw(end+1)=s(4);
end
fclose(fid);

guidata(ha.figure1, ha);
HPFList(ha);
lbList_Callback(hObject, eventdata, ha)
% --------------------------------------------------------------------
function mFileSave_Callback(hObject, eventdata, ha)
[fname,pname] = uiputfile([ha.currdir, '\*.m'],'Save List');
fid = fopen([pname '\' fname], 'w');
if fid == -1, disp('File not found'); return; end
if fname == 0, return; end
ha.currdir = pname;
guidata(ha.figure1, ha);

for ii=1:length(ha.HPFx)
    fprintf(fid,'%f,%f,%f,%f\n',ha.HPFx(ii),ha.HPFm(ii),ha.HPFc(ii),ha.HPFlw(ii));
end
fclose(fid);

