% D01READ Read data from ESP380 transient hack
%
%   data = ESPtransread(filename,...)
%   [ax,data] = ESPtransread(filename,...);
%   [ax,data,dsc] = ESPtransread(filename,...);
%
%   Reads data and axis information from Bruker
%   ESP380 transient hack *.dat.   Aditional 
%   arguments  are  always  coming  in pairs 
%   field, value. Fields are:
%     'Dataset', [real dataset, imag dataset] 
%   ax contains fields x,y,z  depends on data 
%   dimension and may contain other fields as 
%   well. dsc is a cell array  of description 
%   strings.

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de

% boep,   20-jan-04, MPI  

function varargout=ESPtransread(filename, varargin)

% [h,msg] = fopen(filename, 'r', 'ieee-le');
% fread(h, 4,'char')
% p = fread(h, 100,'float');
% p(abs(p)<1e5)
% fclose(h);

[h,msg] = fopen(filename, 'r', 'ieee-be');
error(msg);

uint8 Dims;

% read Field
CenterField = fread(h, 1, 'double');
FieldSpan   = fread(h, 1, 'double');
fldresol    = fread(h, 1, 'int32');

wavedesc      = char(fread(h, 16, 'char')).';
template_name = char(fread(h, 16, 'char')).';

enum_byte_word = fread(h, 1, 'int16');

enum_hi_lo_first = fread(h, 1, 'int16');

wavedecriptor = fread(h, 1, 'long');
usertext = fread(h, 1, 'long');
res_desc1 = fread(h, 1, 'long');

trigtime_array = fread(h, 1, 'long');
ristime_array = fread(h, 1, 'long');
res_array1  = fread(h, 1, 'long');
wave_array1 = fread(h, 1, 'long');
wave_array2 = fread(h, 1, 'long');
res_array2 = fread(h, 1, 'long');
res_array3 = fread(h, 1, 'long');

instrument_name = char(fread(h, 16, 'char')).';
instrument_number = fread(h, 1, 'long');

trace_label = char(fread(h, 16, 'char')).';
reserved = fread(h, 2, 'int16');

wave_array_count = fread(h, 1, 'long');
pts_per_screen   = fread(h, 1, 'long');
first_valid_point = fread(h, 1, 'long');
last_valid_point = fread(h, 1, 'long');
first_point = fread(h, 1, 'long');
sparcing_factor = fread(h, 1, 'long');
segment_index = fread(h, 1, 'long');
sub_array_count = fread(h, 1, 'long');
sweeps_per_acq = fread(h, 1, 'long');
obs = fread(h, 1, 'long');
vertical_gain = fread(h, 1, 'float');
vertical_offset = fread(h, 1, 'float');
max_value = fread(h, 1, 'float');
min_value = fread(h, 1, 'float');
nominal_bits = fread(h, 1, 'int16');
nom_subarray_count = fread(h, 1, 'int16');
horiz_interval = fread(h, 1, 'float');
horiz_offset = fread(h, 1, 'double');
pixel_offset = fread(h, 1, 'double');
vertunit = char(fread(h, 48, 'char')).'; 
horunit = char(fread(h, 48, 'char')).';
res34 = fread(h, 2, 'int16');
trigger_time = fread(h, 16, 'char');
acq_duration = fread(h, 1, 'float');
rec_type = fread(h, 1, 'int16');
processing_done = fread(h, 1, 'int16');
res5 = fread(h, 1, 'int16');
RIS_sweep = fread(h, 1, 'int16');
timebase  = fread(h, 1, 'int16');
vert_coupling = fread(h, 1, 'int16');
probe_att = fread(h, 1, 'float');
fixed_vert_gain = fread(h, 1, 'int16');
bandwidth_limit = fread(h, 1, 'int16');
vertical_vernie = fread(h, 1, 'float');
acq_vert_offset = fread(h, 1, 'float');
wave_source     = fread(h, 1, 'int16');

if usertext > 0
    usertextstr = char(fread(h, usertext, 'char')).';
else
    usertextstr = '';
end

fseek(h,wavedecriptor+20, -1);

ax.x = linspace(CenterField-FieldSpan/2, CenterField+FieldSpan/2, fldresol).';
ax.y = ([1:wave_array_count]-1).' * horiz_interval + horiz_offset;

% the data start 0x16e (366)
Data = fread(h, fldresol*wave_array_count,'int16')*vertical_gain+vertical_offset;

spec = reshape(Data, wave_array_count, fldresol).';

ax.xlabel = 'Field, G';
ax.ylabel = 'Time, s';
ax.type = 'data';

% close data file
fclose(h);

dsc.c0 = ['Field pts: ', num2str(fldresol)];
dsc.c1 = ['Transient pts: ', num2str(wave_array_count)];

% assign output depending on number of output arguments
% requested
switch nargout
 case 1,
   varargout = {spec};
 case 2,
   varargout = {ax, spec};
 case 3,
   varargout = {ax, spec, dsc};
end