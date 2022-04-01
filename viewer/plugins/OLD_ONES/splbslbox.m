function varargout = splbslbox(varargin)
% SPLBSLBOX Application M-file for splbslbox.fig
%    FIG = SPLBSLBOX launch splbslbox GUI.
%    SPLBSLBOX('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 13-Jan-2005 13:10:59

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    handles.PTSidx = [];
    handles.PTSx = [];
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
function varargout = FileUpdate(varargin)


% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, ha, varargin)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 88888)']);

hand = guidata(ha.hh);
axx = hand.src.ax.x;

for ii=1:length(x)
    [mm,idx] = min(abs(axx-x(ii)));
    ha.PTSidx(end+1) =idx;
    ha.PTSx(end+1)   =axx(idx);
end

[ha.PTSx, idx] = sort(ha.PTSx);
ha.PTSidx      = ha.PTSidx(idx);

guidata(ha.figure1, ha);

PointsList(ha);
Plot(ha)

% --------------------------------------------------------------------
function varargout = Plot(handles)
handles = guidata(handles.figure1);

hand = guidata(handles.hh);
hand.out = {};
tt = hand.src;

if get(handles.tbAvSubtracrt, 'Value')
    tt.ax.filt     = 'savgol';
    tt.ax.filtpar1 = '1D';
    tt.ax.tc       = str2num(get(handles.eAverage, 'String'));
    tt.ax.pol      = 2;
    tt.ax.Color     = 'b';
    [tt1.y, tt.ax] = processing(tt.y, tt.ax);
    tt.y = tt.y - tt1.y;
    tt.title = 'AvDiff';
end

hand.out{1} = tt;

if length(handles.PTSidx) > 0
    spl.y    = tt.y(handles.PTSidx,:);
    spl.ax.x = tt.ax.x(handles.PTSidx);
    spl.ax.Marker    = '*';
    spl.ax.Color     = 'r';
    spl.ax.LineStyle = 'none';
    spl.title = 'Marker';
    hand.out{end+1} = spl;
    
    avery = (tt.y(handles.PTSidx,:) + tt.y(handles.PTSidx-1,:) + tt.y(handles.PTSidx+1,:))/3;
    
    for ii=1:size(tt.y,2)
        switch get(handles.pmMethod, 'Value')
        case 1
            [p]=polyfit(tt.ax.x(handles.PTSidx), avery(:,ii),5);
            spl1.y(:,ii) = polyval(p, tt.ax.x);
        case 2
            
        case 3
            spl1.y(:,ii) = spline(tt.ax.x(handles.PTSidx), ...
                avery(:,ii), tt.ax.x);
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
    spldif.ax.xlabel =     tt.ax.xlabel;
    spldif.ax.ylabel =     tt.ax.ylabel; 
    spldif.y = tt.y - spl1.y;
    hand.out{end+1} = spldif;
end

[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);



% --------------------------------------------------------------------
function varargout = pbAvSubtracrt_Callback(h, eventdata, handles, varargin)
Plot(handles)

% --------------------------------------------------------------------
function varargout = eAverage_Callback(h, eventdata, handles, varargin)
Plot(handles)

% --------------------------------------------------------------------
function varargout = lbPoints_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = pbClear_Callback(h, eventdata, handles, varargin)
handles.PTSidx = [];
handles.PTSx = [];
guidata(handles.figure1, handles);
PointsList(handles)
Plot(handles)

function PointsList(handles)
Str = {};

for ii=1:length(handles.PTSidx)
    Str{end+1} = sprintf('[%5d]   %5.3f',handles.PTSidx(ii),handles.PTSx(ii));
end


set(handles.lbPoints,'String',Str);

% --------------------------------------------------------------------
function varargout = pbChange_Callback(h, eventdata, ha, varargin)
hand = guidata(ha.hh);
[path,name,ext] = fileparts(get(ha.hh, 'FileName'));

[x,y]=eval([name '(''GetPoint'', hand, 2)']);

if length(x)~=2, return; end;

hand = guidata(ha.hh);
axx = hand.src.ax.x;

[mm,pidx] = min(abs(ha.PTSx-x(1)));

[mm,idx] = min(abs(axx-x(2)));
ha.PTSidx(pidx) =idx;
ha.PTSx(pidx)   =axx(idx);

[ha.PTSx, idx] = sort(ha.PTSx);
ha.PTSidx      = ha.PTSidx(idx);

guidata(ha.figure1, ha);

PointsList(ha);
Plot(ha)

% --------------------------------------------------------------------
function varargout = cbPointBased_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pmMethod_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pbAutoPeaks_Callback(h, eventdata, ha, varargin)

hand = guidata(ha.hh);
axx = hand.src.ax.x;

tt = hand.src;

if get(ha.tbAvSubtracrt, 'Value')
    tt.ax.filt     = 'savgol';
    tt.ax.filtpar1 = '1D';
    tt.ax.pol      = 2;
    tt.ax.tc       = str2num(get(ha.eAverage, 'String'));
    [tt1.y, tt.ax] = processing(tt.y, tt.ax);
    tt.y = tt.y - tt1.y;
end

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

PointsList(ha);
Plot(ha)



% --------------------------------------------------------------------
function varargout = eTolerance_Callback(h, eventdata, handles, varargin)
Plot(handles)

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
