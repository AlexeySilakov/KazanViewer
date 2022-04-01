function dy = svd_tmpfunc(t, y) 
global k0 k12 k23; 
dy = zeros(3,1); 
dy(1)=-k12*y(1) ; 
dy(2)=+k12*y(1)-k23*y(2) ; 
dy(3)=k23*y(2) ; 
