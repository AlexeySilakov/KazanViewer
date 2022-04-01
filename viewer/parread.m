% PARREAD Read data from ascii files
%
%   data = parread(filename, ...)
%   [ax,data] = parread(filename, ...);
%   [ax,data,dsc] = parread(filename, ...);
%
%   Reads data and axis information from dat
%   files. Aditional arguments are always 
%   coming in pairs field, value. 
%   ax contains fields x,y,z depends on data 
%   dimension and may contain other fields 
%   as well. dsc is a cell array of description 
%   strings. 

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% Silakov Alexey, 09-dec-09, MPI  


function varargout = parread(filename, segments, representation);
fid = fopen(filename, 'r');
if fid == -1, error('sorry, file not found'); return; end

k = 1; a='';
strings = {};
actionin = 0;
actioncount = 0;
segcount = 0;
totpoints = 0;
expcount = 0;
expin = 0;
datanames = [];


while ~feof(fid)
    a = fgets(fid);
    if isempty(a), continue; end
    %%% start a new description
    if sum(strfind(a, '<Action')), actioncount = actioncount+1; actionin = actioncount;
        points(actioncount) = 0;
        while actionin
            a = fgets(fid);
            if sum(strfind(a, 'Total Points')), points(actioncount) = sscanf(a, 'Total Points=%d');  end
            if sum(strfind(a, '</Action')), actionin = 0; break; end
            if sum(strfind(a, 'Name')), action_names{actioncount} = a(6:(end-2));  end
        end
    end
    %%%% Experiment section to get time and date
    % TimeAcquired=4:47:30 PM
    % DateAcquired=Wednesday, August 21, 2019

    if sum(strfind(a, '<Experiment')), expcount = expcount+1; expin = expcount;
        Time{expcount} = {};
        Date{expcount} = {};
        while expcount
            a = fgets(fid);
            if sum(strfind(a, 'TimeAcquired')), Time{expcount} = strtrim(a(numel('TimeAcquired='):end)); end
            if sum(strfind(a, 'DateAcquired')), Date{expcount} = strtrim(a(numel('DateAcquired='):end)); end
            if sum(strfind(a, '</Experiment')), expin = 0; break; end
%             if sum(strfind(a, 'Name')), action_names{expcount} = a(6:(end-2));  end
        end
    end    
    %%%% Data block
    if sum(strfind(a, '<Segment'))
        segcount = segcount+1;
        % Type=2
        % Version=1
        % Definition=
        a = fgets(fid);
        type(segcount) = sscanf(a, 'Type=%d');
        a = fgets(fid);
        Version(segcount) = sscanf(a, 'Version=%d');
        %%%%%%%%%%% get data column names
        a = fgets(fid);
        poscom = strfind(a(1:end), ',');
        poscom = [9, poscom-1, length(a)];
        ndatas = (length(poscom)-1);
        for ii = 1:(length(poscom)-1)
            namev = a((poscom(ii)+3):poscom(ii+1));
            nnn = strfind(namev, '(');
            namev(nnn) = '_';
            nnn = strfind(namev, ')');
            namev(nnn) = '_';
            nnn = strfind(namev, ' ');
            namev(nnn) = '_';
            nnn = strfind(namev, '#');
            namev(nnn) = 'n';
            if strcmp(namev(1), '0'), ndatas = ndatas-1; continue; end
            nnn = strfind(namev, '0');
            namev(nnn) = 'u';
            datanames = setfield(datanames, namev, ii);
        end
        insegment = 1;
%         F = textscan(fid, '%f', ndatas*sum(points) , 'delimiter', ',');
        F = textscan(fid, '%f', 'delimiter', ',');
        In = reshape(F{1}, ndatas, size(F{1}, 1)/ndatas).';
    end
end
fclose(fid);

% dlmread('filename', delimiter, R, C)

switch representation
    case 'I(E)'
        datanum = datanames.I_A_;
        axnum = datanames.E_V_;
        xlabel = 'E, V';
        ylabel = 'I, A';
    case 'E(T)'
        datanum = datanames.E_V_;
        axnum = datanames.Elapsed_Time_s_;
        xlabel = 'Time, s';
        ylabel = 'E, V';
    case 'I(T)'
        datanum = datanames.I_A_;
        axnum = datanames.Elapsed_Time_s_;
        xlabel = 'Time, s';
        ylabel = 'I, A';
end

if segments>=0
    if min(segments)>max(In(:, datanames.Segment_n)), segments = max(In(:, datanames.Segment_n)); end
    datapoints_choose = [];
    for ii = 1:length(segments)
        fff = find(In(:, datanames.Segment_n)==segments(ii) );
        datapoints_choose = [datapoints_choose; fff];
    end
    ax.x = In(datapoints_choose, axnum);
    Y = In(datapoints_choose, datanum);

else
    ax.x = In(:, axnum);
    Y = In(:, datanum);
    if isfield(datanames, 'Segment_n')
        ax.Segment = In(:, datanames.Segment_n);
    else
        ax.Segment = [];
    end
end
    
ax.type = 'data';
ax.xlabel = xlabel;
ax.title = '';
dsc.ASCII='PAR loader for KAZAN viewer';
% if fstcol,dsc.FirstColumn='X';else dsc.FirstColumn='Y';end
dsc.Representation = ['Treated as ', representation];
dsc.size_x = num2str(size(Y,1));
dsc.size_y = num2str(size(Y,2));
dsc.NActions = num2str(length(points));
for ii = 1:length(Time)
    dsc = setfield(dsc, ['Time', num2str(ii-1)], Time{ii});
    dsc = setfield(dsc, ['Date', num2str(ii-1)], Date{ii});
end

for ii = 1:length(points)
    dsc = setfield(dsc, ['DataSize_A', num2str(ii-1)], num2str(points(ii)));
    dsc = setfield(dsc, ['DataName_A', num2str(ii-1)], action_names{ii});
end
dsc.ActionSize = num2str(sum(points));
% if ~isempty(strings)
%     for n = 1:size(strings, 1)
%       dsc = setfield(dsc, ['c', num2str(n)], strings{n});
%     end
% end

% assign output depending on number of output arguments
% requested
switch nargout
  case 1,
    varargout = {Y};
  case 2,
    varargout = {ax, Y};
  case 3,
    varargout = {ax, Y, dsc};
end