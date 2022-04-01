% TECMAGREAD Read TecMag data files
%
%   data = TecMagread(filename,...)
%   [ax,data] = TecMagread(filename,...);
%   [ax,data,dsc] = TecMagread(filename,...);
%
%   Reads data and axis information from Bruker
%   files  DSC/DTA  and  PAR/SPC.   Aditional 
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

% bme,   21-dec-03, MPI  

function varargout=TecMagread(filename, varargin)

fid=fopen(char(filename),'r', 'ieee-be');

dsc = [];
nrecords = fread(fid, 1, 'int');
dsc.nrecords = num2str(nrecords);
if double(fread(fid, 1, 'char'))
  dsc.dtype = 'Frequency, Hz';
else
  dsc.dtype = 'Time, us';
end
fseek(fid,50,'bof');
npoints1D = fread(fid, 1, 'int');
dsc.npoints1D = num2str(npoints1D);
fread(fid, 1, 'int');
dsc.sample = char(fread(fid, 32, 'char'))';
dsc.sequence = char(fread(fid, 32, 'char'))';
fread(fid, 32, 'char');
dsc.mfield = [num2str(fread(fid, 1, 'double')), ' T'];
ax.freq1 = fread(fid, 1, 'double');
dsc.frequency1 = [num2str(ax.freq1), ' MHz'];
swidth = fread(fid, 1, 'double');
dsc.spwidth = [num2str(swidth), ' Hz'];
dw = fread(fid, 1, 'double')*1E6;
dsc.dw = [num2str(dw), ' us'];
dsc.acqtime = [num2str(fread(fid, 1, 'double')), ' s'];
fread(fid, 1, 'double');
ax.freq2 = fread(fid, 1, 'double');
dsc.frequency2 = [num2str(ax.freq2), ' MHz'];

fseek(fid,1024,'bof');
comment=fread(fid, 1024,'char');
comment = char(comment(comment~=0))';
data = fread(fid, 2*nrecords*npoints1D,'float');
fclose(fid);

data1 = reshape(data, 2, npoints1D).';

spec = data1(:,2)+j*data1(:,1);

ax.x = dw*[1:size(spec,1)].';
ax.xlabel = dsc.dtype;
k = 1;
while ~isempty(comment)
  [line, comment] = strtok(comment, char([10,13]));
  dsc = setfield(dsc, ['c', num2str(k)], trim(line));
  k = k + 1;
  comment = trim(comment(2:end));
end
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




% # MacNMR binary scripts
% #
% # Tecmag MacNMR 5.0 binary files have no naming convention and parameters 
% # are contained in the header and tailer. The 1kB Notes string (contents of 
% # the MacNMR Notes window) in the header is also examined for user-defined 
% # parameters. Sometimes I add parameters to the Notes window called 'lc6' or
% # 'ept' for the number of points per echo and 'tau' for the time between 
% # pulses (or echoes).
% 
% [Caption]
% /* Use this section to supply a file dialog caption */
% // fdlog.dlgName$="Your caption";
% return;
% 
% [Main]
% file.nBytes=qdnBytes;      /* bytes per value from file type settings */
% file.ByteOrder=qdbOrder;  /* byte order */
% file.signed=qdsigned;     /* signed or unsigned */
% file.open(%B%P);
% qdal2=file.read();      /* The number of records is given at the beginning of the header. */
% %S="Count 2D"; %A=$(qdal2); logPar f;
% qdft1=file.read(1);     /* domain: 0=time, >0=freq */
% file.setpos(50);
% qdal1=file.read();        /* points per record */
% %S="Points 1D"; %A=$(qdal1); logPar f;
% qdnpt=qdal1*qdal2;
% file.read();
% file.cread(M,32);
% file.cread(M,32);
% %S=Sequence; %A=%M; logPar f; /* pulse sequence name */
% file.cread(M,32);
% %S=Field; %A=$(file.dread()); logPar f; /* magnetic field value */
% qdsf=file.dread();
% %S="Obs Freq"; %A=$(qdsf); logPar f; /* spectrometer frequency */
% %S="SW +/-"; %A=$(file.dread()); logPar f; /* spectral width (+/-) */
% %S="Dwell 1D"; qdx=file.dread(); qdw1=qdx; logParu f; /* D1 dwell time */
% %S="Acq. Time"; qdx=file.dread(); logParu f; /* D1 acquisition time */
% 
% qdlb1=file.dread();
% %S="LB 1D"; %A=$(qdlb1); logPar f; /* D1 line broadening */
% %S="F2 freq"; %A=$(file.dread()); logPar f; /* Decoupling frequency */
% %S="Shim Units"; %A=$(file.dread()); logPar f; /* Shim units */
% %S="Last Delay"; qdx=file.dread(); logParu f;/* sequence delay */
% qdnn=2*qdnBytes*qdnpt;
% qdnn+=2474; file.setpos(qdnn);
% qdac=file.read();
% %S="Scan Count"; %A=$(qdac); logPar f; /* scan count */
% qdnn+=94; file.setpos(qdnn);
% %S="Acq. Points"; %A=$(file.read()); logPar f;
% qdnn=2*qdnBytes*qdnpt;
% qdnn+=2610;
% file.setpos(qdnn);
% qdw2=file.dread();
% %S="Dwell 2D"; qdx=qdw2; logParu f;
% file.setpos(1024);
% file.cread(Z,1024);             /* Load notes into %Z */
% file.close();
% /* Sort through notes for parameter declarations in <name>=<value> format */
% if(%Z!="") {
%   loop (i,1,%[%Z]-1) {
%     %M=%[%Z,#i]; %A="";
%     loop (j,1,%[%M]) {
%       if((%[%M,j:j]=="=") && (j<%[%M])) {%A=%[%M,#2,=]; break};
%     };
%     if(%A!="") {
%       %S=%[%M,#1,=]; logPar f;
%     };
%   };
% };
% if(%Z!="") run.section(,WriteCom);
% 
% [WriteCom]
% /* Put all comments into a label that will appear on the data worksheet */
% nn=0; qd0=0;
% loop (i,1,%[%Z]) {
%   if(asc(%[%Z,i:i])==13) {
%     nn++; qd$(nn)=i;
%     label -s -n com$(nn) %[%Z,qd$(nn-1)+1:qd$(nn)-1];
%   };
% };
% if(nn<%[%Z]) {
%   nn++; label -s -n com$(nn) %[%Z,qd$(nn-1)+1:%[%Z]];
% };
% loop (i,1,nn) {
%   %R=com$(i).text$;
%   if(i==1) %Z=%R;
%   else {
%     %Z=%Z
% %R;
%   };
%   label -r com$(i);
% };
% run.section(nmrtools\Data,WriteCom);
% 
% [Custom]
% /* This section runs just before the data is prepared for plotting.
%    Use it to for custom parameter replacements or data manipulation. */
% // return;
% %B=lc6; getValue f;
% if(%B!="lc6") {
%   qdal1=%B; qdal2=int(qdnpt/qdal1);
%   %B=tau; getValue f;
%   if(%B!="tau") {%A=%B; Parse; qdw2=qdx};
% };

