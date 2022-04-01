function varargout = description(varargin)
% DESCRIPTION Application M-file for description.fig
%    FIG = DESCRIPTION launch description GUI.
%    DESCRIPTION('callback_name', ...) invoke the named callback.

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.0 16-Dec-2003 16:46:47

if nargin == 0  % LAUNCH GUI
    token = 'DESCRIPTION_open';
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
%| 'Description_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.Description, handles.slider2. This
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
function Description_DeleteFcn(hObject, eventdata, handles)
[path,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''DeletePlugins'', handles.Description, handles.hh)']);

% --------------------------------------------------------------------
function varargout = lbDescription_Callback(h, eventdata, handles, varargin)

function FileUpdate(handles)
hand = guidata(handles.hh);
if isfield(hand.src, 'dsc')
    names = fieldnames(hand.src.dsc);
    str = {};
    for k = 1:size(names)
        str{end+1}=sprintf('%-16s',[names{k} ':'],getfield(hand.src.dsc, names{k}));
    end
    set(handles.lbDescription, 'String', str, 'Value', 1);
end 

% --------------------------------------------------------------------
function varargout = pbSeq_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
if isfield(hand.src, 'dsc')
    if isfield(hand.src.dsc, 'Psd1')
        figure
        pos = 2;
        slength = 35;
        shift = 68;
        hold on;
        idx = findstr(getfield(hand.src.dsc, 'Psd5'), '}'); % for new Bruker software
        if isempty(idx)
            idx = 1;
        end
        idx = idx+1;
        str =  getfield(hand.src.dsc, 'Psd5');
        val = str2num(str(idx:end));
        maxval = max(val(1:slength))+max(val(slength:slength*2))+50;
        
        for i=1:50 
            %%%%%%%% number of channels is not defined. But it is less than 50 for sure

            if ~isfield(hand.src.dsc, ['Psd' num2str(i)])
                break;
            end
           start = 0;
           str =  getfield(hand.src.dsc, ['Psd' num2str(i)]);
           val1 = str2num(str(idx:end));
        
%             val = getfield(hand.src.dsc, ['Psd' num2str(i)]);
%             val1 = str2num(val);
            for j=1:10
                if val1(slength+j) > 0
                    plot([start val1(pos+j)], i+[0 0]);
                    plot([val1(pos+j) val1(pos+j)], i+[0 .5]);
                    start = val1(pos+j)+val1(slength+j);
                    start1=val1(pos+j);
                    plot([start1 start], i+[0.5 0.5]);
                    text(start1,i+0.7,[num2str(start-start1),'ns'])
                    plot([start start], i+[0.5 0]);
                end
            end
            if start > 0
                plot([start maxval], i+[0 0]);
                %         text(20, -2*i+.3, num2str(i))  
            end
        end
        
        
        xlabel('Sequence, ns');
        ylabel('Channels');
        axis tight
    elseif isfield(hand.src.dsc, 'XPD1')
        figure; hold on
        for i=1:16
            if isfield(hand.src.dsc, ['XPD',num2str(i)])               
                start = 0;
                val = getfield(hand.src.dsc, ['XPD' num2str(i)]);
                val1 = str2num(val);
                if length(val1)>4
                    ppulse=floor((length(val1)+3)/4)-1;
                    val1(end+1:end+4) = 0;
                    for j=1:ppulse
                        start1 = val1(1+j*4)*8;
                        plot([start,start1], i+[0 0]);
                        text(start+10,i+0.2,[num2str(start1-start),'ns'])
                        start = start1;
                        plot([start,start], i+[0 .5]);
                        start1 = 8*(val1(2+j*4)-1)+start;
                        plot([start,start1], i+[0.5 0.5]);
                        text(start,i+0.7,[num2str(start1-start),'ns'])
                        plot([start1 start1], i+[0.5 0]);
                        start=start1;
                    end
                end
            end
        end
        xlabel('Sequence, ns');
        ylabel('Channels');
        title(hand.src.filename);
        axis tight        
    end
end
