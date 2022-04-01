function varargout = transient(varargin)
% TRANSIENT Application M-file for transient.fig
%    FIG = TRANSIENT launch transient GUI.
%    TRANSIENT('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 22-Jul-2005 16:17:03

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


% --------------------------------------------------------------------
function varargout = Transient_DeleteFcn(h, eventdata, handles, varargin)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.transient, handles.hh)']);

function FileUpdate(handles)


% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

plothis(handles);

% --------------------------------------------------------------------
function plothis(handles)
handles = guidata(handles.Transient);

hand = guidata(handles.hh);
ax.x = hand.src.ax.x;
ax.y = hand.src.ax.y;
if ~isfield(hand.src, 'y'), return; end;
y = hand.src.y;

hand.out = {};
hand.out{1}.ax = hand.src.ax;

if get(handles.pbTrans,'Value')==1, 
    y = y.'; 
    hand.out{1}.ax.x = hand.out{1}.ax.y;
    hand.out{1}.ax.y = [0];
end

hand.out{1}.ax.xlabel = hand.out{1}.ax.ylabel;
hand.out{1}.ax.y = 1;
hand.out{1}.ax.ylabel = '';

BLstart = get(handles.eBLstart, 'String');
BLend = get(handles.eBLend, 'String');
Sstart = get(handles.eSstart, 'String');
Send = get(handles.eSend, 'String');
Lstart = get(handles.eLstart, 'String');
Lend = get(handles.eLend, 'String');

BL = [str2num(BLstart):str2num(BLend)];

LL = [str2num(Lstart):str2num(Lend)];
    
lb = sum(y(LL, :), 1)./size(LL, 2);

% get whole data, alsi 29.05.07 %%%%%%%%%%%%%
if str2num(Sstart)==0 & str2num(Send)==0
%hand.out{1}.y = (sum(y(:, SI), 2)- sum(lb(SI), 2))/size(SI, 2) -...
%    (sum(y(:, BL), 2)-sum(lb(BL), 2))/size(BL, 2);
    ty = (y - lb(ones(size(y, 1), 1), :));
    ds = sum(ty(:, BL), 2)/length(BL);
    hand.out{1}.y = ty - ds(:, ones(size(ty, 2), 1));
    
    if get(handles.pbTrans,'Value')==1, 
        hand.out{1}.ax.x = hand.src.ax.y;
        hand.out{1}.ax.xlabel = hand.src.ax.ylabel;
        hand.out{1}.ax.y = hand.src.ax.x;
        hand.out{1}.ax.ylabel = hand.src.ax.xlabel;
    else
        hand.out{1}.ax.x = hand.src.ax.x;
        hand.out{1}.ax.xlabel = hand.src.ax.xlabel;
        hand.out{1}.ax.y = hand.src.ax.y;
        hand.out{1}.ax.ylabel = hand.src.ax.ylabel;
    end

else
   SI = [str2num(Sstart):str2num(Send)];
    hand.out{1}.y = (sum(y(:, SI), 2)- sum(lb(SI), 2))/size(SI, 2) -...
    (sum(y(:, BL), 2)-sum(lb(BL), 2))/size(BL, 2);
end

if get(handles.cbInvert,'Value')==1
    hand.out{1}.y = -hand.out{1}.y;
end

guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);



% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = pbTrans_Callback(h, eventdata, handles, varargin)


% --- Executes on button press in cbInvert.
function cbInvert_Callback(hObject, eventdata, handles)
% hObject    handle to cbInvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbInvert


