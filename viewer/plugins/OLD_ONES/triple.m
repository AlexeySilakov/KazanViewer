function varargout = triple(varargin)
% TRIPLE Application M-file for triple.fig
%    FIG = TRIPLE launch triple GUI.
%    TRIPLE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 18-Jun-2005 18:35:54

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
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
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

PlotResult(handles);

% --------------------------------------------------------------------
function PlotResult(handles)
handles = guidata(handles.figure1);

hand = guidata(handles.hh);
ax.x = hand.src.ax.x;
ax.y = hand.src.ax.y;
if ~isfield(hand.src, 'y'), return; end;
y = hand.src.y;

hand.out = {};
hand.out{1}.ax = hand.src.ax;

bl1 = eval(['[',get(handles.eDim1, 'String'),'];']);
bl2 = eval(['[',get(handles.eDim2, 'String'),'];']);

hand.out{1}.y = tripleprocess(hand.src.y, bl1, bl2);

if get(handles.cbInverse,'Value')
  hand.out{1}.y = -hand.out{1}.y;  
%   hand.out{1}.y(hand.out{1}.y<0)=0;
end

guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);


% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)

hand = guidata(handles.hh);
sz = size(hand.src.y);

blfld = str2num(get(handles.eBlfld, 'String'));

set(handles.eDim1, 'String', ['1:', num2str(blfld), ', ', num2str(sz(2)-blfld),':',num2str(sz(2))]);
set(handles.eDim2, 'String', ['1:', num2str(blfld), ', ', num2str(sz(1)-blfld),':',num2str(sz(1))]);


% --------------------------------------------------------------------
function varargout = eBlfld_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = pnReference_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
if get(handles.tbTranspose,'Value')
  handles.ref = hand.src.y(hand.Selection,:).';
else
  handles.ref = hand.src.y(:,hand.Selection);
end
set(handles.eSLice, 'String', num2str(hand.Selection));
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = pbProcess_Callback(h, eventdata, handles, varargin)

handles = guidata(handles.figure1);
hand = guidata(handles.hh);
if ~isfield(hand.src, 'y'), return; end;

hand.out = {};
if get(handles.tbTranspose,'Value') 
    [ax,y] = kv_trans(hand.src.ax, hand.src.y);
else
    y = hand.src.y;
    ax = hand.src.ax;
end
hand.out{1}.ax = ax;
hand.out{1}.y = -y+handles.ref(:,ones(1,size(y,2)));

if get(handles.cbInverse,'Value')
  hand.out{1}.y = -hand.out{1}.y;  
end
guidata(handles.hh, hand);

[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function varargout = cbInverse_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = tbTranspose_Callback(h, eventdata, handles, varargin)

function [rax, y] = kv_trans(ax, y)
y = y.';
rax.x = ax.y;
rax.y = ax.x;
rax.xlabel = safeget(ax,'ylabel', '?,');
rax.ylabel = safeget(ax,'xlabel', '?,');


% --------------------------------------------------------------------
function varargout = eSLice_Callback(h, eventdata, handles, varargin)

