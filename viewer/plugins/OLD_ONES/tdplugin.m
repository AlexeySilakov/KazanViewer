function varargout = tdplugin(varargin)
% TDPLUGIN Application M-file for tdplugin.fig
%    FIG = TDPLUGIN launch tdplugin GUI.
%    TDPLUGIN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 06-Dec-2004 10:49:19

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
function varargout = slider2_Callback(h, eventdata, handles, varargin)
slice = str2num(get(handles.eSlice,'String'));
step = str2num(get(handles.eStep,'String'));

shift =get(handles.slChange, 'Value');
set(handles.slChange, 'Value', 0.5);
Val = step*(shift - 0.5)*2;

handles = guidata(handles.figure1);
hand = guidata(handles.hh);

hand.src.y(slice,:) = hand.src.y(slice,:)*(1+Val);
guidata(handles.hh, hand);

PlotThis(handles);

% --------------------------------------------------------------------
function varargout = slSlice_Callback(h, eventdata, handles, varargin)
slice = str2num(get(handles.eSlice,'String'));
shift =get(handles.slSlice, 'Value');
set(handles.slSlice, 'Value', 0.5);
Val = (shift - 0.5)*2;
set(handles.eSlice,'String', num2str(Val+slice));

% --------------------------------------------------------------------
function varargout = togglebutton1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = eS1_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = eS2_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = eS3_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = eS4_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = eS5_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = eJ_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function varargout = pbRefresh_Callback(h, eventdata, handles, varargin)
PlotThis(handles);

% --------------------------------------------------------------------
function PlotThis(handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
hand.out = {};
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));

J = str2num(get(handles.eJ,'String'));
x(1) = str2num(get(handles.eS1,'String'));
x(2) = str2num(get(handles.eS2,'String'));
x(3) = str2num(get(handles.eS3,'String'));
x(4) = str2num(get(handles.eS4,'String'));
x(5) = str2num(get(handles.eS5,'String'));

tr   = JFITfunc(x,hand.src.ax.x, J);
hand.out{end+1} = hand.src;
dat = eval([name,'(''GetSourceSlice'', hand)']);
norm2= max(sum(tr,2));
[norm,idx] = max(abs(dat.y));
norm = dat.y(idx);
hand.out{end}.y = hand.out{end}.y;
hand.out{end}.ax.Marker = '.';
hand.out{end}.title     = ['Data'];

nn  = norm / norm2;
hand.out{end+1}.y = sum(tr,2)*nn; 
hand.out{end}.ax.x = hand.src.ax.x;
hand.out{end}.ax.y = [0];
hand.out{end}.ax.Color = 'm';
hand.out{end}.ax.xlabel = 'Temperature, K';
hand.out{end}.ax.ylabel = 'Spin State, []';
hand.out{end}.title     = ['Sum'];
hand.out{end}.ax.LineWidth  = 2;

for kk=1:5
   hand.out{end+1}.y = tr(:,kk)*nn; 
   hand.out{end}.ax.x = hand.src.ax.x;
   hand.out{end}.ax.y = [kk];
   hand.out{end}.ax.xlabel = 'Temperature, K';
   hand.out{end}.ax.ylabel = 'Spin State, []';
   hand.out{end}.title     = ['State ',num2str(kk)];
end


guidata(handles.hh, hand);
eval([name '(''SetDataSource'', 2, hand)']);
% --------------------------------------------------------------------

function res = JFITfunc(x,T,J)
S = [1:5]; 
Si = [0:5];

sz  = length(T);
idx = ones(1,length(S));
SS(1:sz,1)=1;SS(1:sz,2)=2;SS(1:sz,3)=3;SS(1:sz,4)=4;SS(1:sz,5)=5;
for kk=1:6, SSi(1:sz,kk)=Si(kk); end


amp = x;

nn = sum(((2*SSi+1).*exp(-SSi.*(SSi+1)*J./T(:, [1,1,1,1,1,1]))),2);
res = amp(ones(sz,1),:).*exp(-SS.*(SS+1).*J./T(:, idx))./nn(:,idx)./T(:, idx);


% --------------------------------------------------------------------
function varargout = eSlice_Callback(h, eventdata, handles, varargin)





