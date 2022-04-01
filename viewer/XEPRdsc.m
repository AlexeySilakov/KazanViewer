function [par] = XEPRdsc(filename)
h = fopen(filename);
if h < 0
    error(['File not found: ',filename])
end
olda = '';
par = [];
forbidden = '~!@#$%^&*()./\';
while feof(h)<1
    s = fgetl(h);
    [a,s]=strtok(s);

    if ~isempty(a) & sum(forbidden == a(1))==0
        try 
            if sum('0123456789' == a(1)) == 0
                olda = a;
                s = trim(s(2:end));
                if(s(1)=='''' & s(end)=='''')
                    s = s(2:end-1);
                end
                par = setfield(par, a, s);
            else
                fld = getfield(par, olda);
                fld = trim(strcat(fld, s));
                par = setfield(par, olda, fld);
            end
        end
    end
end
fclose(h);