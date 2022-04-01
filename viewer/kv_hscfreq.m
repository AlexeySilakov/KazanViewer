function Freq = hscfreq(Sys, Exp, phi, theta)
% Freq = hscfreq(Sys, Exp, phi, theta)
%
% Calculates HYSCORE crosspeacks in FD for system: S=1/2 I=1/2, 1
% Freq = [w1, w2]
% 
% Sys: gn, g, A[MHz], Apa.
% Exp: Field[mT], tau
switch Sys.I
   case 0.5
      Bo = Exp.Field*1e-3; % T
      nu_i = nmagn/planck*Bo*Sys.gn(1, 1); % Hz
      gx = Sys.g(1);
      gy = Sys.g(2);
      gz = Sys.g(3);
      
      gn = Sys.gn(1, 1);
      if isfield(Sys, 'TA')
          A = Sys.TA;
      else
          R = erot(Sys.Apa(1, :));
          A = R*diag(Sys.A(1, :)*1e6)*R.'; % Hz
      end
      %tau = Exp.tau;
      k = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];         
      
      ha = nu_i*k + k*A/2;
      hb = nu_i*k - k*A/2;
      
      nu_a = norm(ha)*1e-6;
      nu_b = norm(hb)*1e-6; 
      % 
      % nu_p = nu_a + nu_b;
      % nu_m = nu_b - nu_a;
      % 
      % c = sqrt(abs(nu_i^2 - 1/4*(nu_a - nu_b)^2)/(nu_a*nu_b));
      % s = sqrt(abs(nu_i^2 - 1/4*(nu_a + nu_b)^2)/(nu_a*nu_b));
      % k = 4*c^2*s^2;
      % 
      % if tau==0
      %     amp1 = 1;
      % else
      %     amp1 =k*sin(pi*nu_a*tau)*sin(pi*nu_b*tau);
      % end
      
      w1 = nu_a; w2 = nu_b;% amp1 = c^2*amp1; amp2 = s^2*amp1;
      Freq = [w1 w2];
   case 1
      E = levchick(Sys,phi,theta, Exp.Field);    
      nu_a(1, 1) = abs(E(1) - E(2));
      nu_a(2, 1) = abs(E(2) - E(3));
      nu_a(3, 1) = abs(E(1) - E(3));
      
      nu_b(1, 1) = abs(E(4) - E(5));
      nu_b(2, 1) = abs(E(5) - E(6));
      nu_b(3, 1) = abs(E(4) - E(6));  
      Freq = [nu_a, nu_b];
   otherwise
      Freq = local_exact(Sys,phi,theta, Exp.Field); 
      
      return;
end
% function E = directdiag(Sy, phi, theta, Field)


function E=levchick(Sy, phi, theta, Field)
planck = 6.6261e-034;
nmagn = 5.0508e-027;
pi = 3.141593;
if isfield(Sy, 'TA')
  A = Sy.TA;
else

  R = erot(safeget(Sy, 'Apa', [0,0,0]));
  A = R * diag(Sy.A)*R.';
end


ctheta = cos(theta);
stheta = sin(theta);

cphi = cos(phi);
sphi = sin(phi);

angle(1,1) = stheta*cphi;
angle(2,1) = stheta*sphi;
angle(3,1) = ctheta;
sangle = stheta;

omega_l = Field/planck /1E9 * Sy.gn * nmagn;

B = omega_l * angle;

Gx=angle(1,1)*Sy.g(1);
Gy=angle(2,1)*Sy.g(2);
Gz=angle(3,1)*Sy.g(3);

geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

Axyz(1,:)=(Gx*A(1,1) + Gy*A(1,2) + Gz*A(1,3))./geff;
Axyz(2,:)=(Gx*A(2,1) + Gy*A(2,2) + Gz*A(2,3))./geff;
Axyz(3,:)=(Gx*A(3,1) + Gy*A(3,2) + Gz*A(3,3))./geff;

sq2 = sqrt(2);

Hzeeman(1,:) = -B(3,:);                        % (0,0)
Hzeeman(2,:) = 0;                              % (1,1)
Hzeeman(3,:) = B(3,:);                         % (2,2)
Hzeeman(4,:) = (-B(1,:) - j*B(2,:))/2*sq2;     % (1,0)
Hzeeman(5,:) = 0;                              % (2,0)
Hzeeman(6,:) = (-B(1,:) - j*B(2,:))/2*sq2;     % (2,1)

Hhpf(1,:) = Axyz(3,:);
Hhpf(2,:) = 0;
Hhpf(3,:) = -Axyz(3,:);
Hhpf(4,:) = Axyz(1,:)/2*sq2 + j*Axyz(2,:)/2*sq2;
Hhpf(5,:) = 0;
Hhpf(6,:) = Axyz(1,:)/2*sq2 + j*Axyz(2,:)/2*sq2;

Hhpf=Hhpf.*0.5;  %% e-spin

% not tested yet 
sz= size(angle,2);
if isfield(Sy,'TQ')
    Q = Sy.TQ;
    HQuad(1,1:sz) = Q(3,3)/2;
    HQuad(2,1:sz) = -Q(3,3);
    HQuad(3,1:sz) = Q(3,3)/2;
    HQuad(4,1:sz) = (Q(1,3)+j*Q(2,3))*sq2/2;
    HQuad(5,1:sz) = (Q(1,1)-Q(2,2))/2 + j*Q(1,2);
    HQuad(6,1:sz) = -(Q(1,3)+j*Q(2,3))*sq2/2;
    H =  Hzeeman + Hhpf + HQuad;
    H1 = Hzeeman - Hhpf + HQuad;
elseif isfield(Sy,'Q')
    
    RQ = erot(safeget(Sy, 'Qpa', [0,0,0]));
    Q = RQ * diag(Sy.Q)*RQ.';
    HQuad(1,1:sz) = Q(3,3)/2;
    HQuad(2,1:sz) = -Q(3,3);
    HQuad(3,1:sz) = Q(3,3)/2;
    HQuad(4,1:sz) = (Q(1,3)+j*Q(2,3))*sq2/2;
    HQuad(5,1:sz) = (Q(1,1)-Q(2,2))/2 + j*Q(1,2);
    HQuad(6,1:sz) = -(Q(1,3)+j*Q(2,3))*sq2/2;
    H =  Hzeeman + Hhpf + HQuad;
    H1 = Hzeeman - Hhpf + HQuad;

else
    H =  Hzeeman + Hhpf;
    H1 = Hzeeman - Hhpf;
end

e1 = p3eq(H);
e2 = p3eq(H1);
% wa(1,:) = abs(e1(1,:)-e1(2,:));
% wa(2,:) = abs(e1(2,:)-e1(3,:));
% wa(3,:) = abs(e1(1,:)-e1(3,:));
% 
% wb(1,:) = abs(e2(1,:)-e2(2,:));
% wb(2,:) = abs(e2(2,:)-e2(3,:));
% wb(3,:) = abs(e2(1,:)-e2(3,:));
E = [e1; e2];
return

% input is (0,0)  (1,1)  (2,2)  (1,0)  (2,0)  (2,1)
%            1      2      3      4      5      6
% see Muha article for details
function w = p3eq(H)
sh = (H(1,:)+H(2,:)+H(3,:))/3;
H(1,:) = H(1,:) - sh;
H(2,:) = H(2,:) - sh;
H(3,:) = H(3,:) - sh;
W1=H(4,:).*conj(H(4,:));
W2=H(5,:).*conj(H(5,:));
W3=H(6,:).*conj(H(6,:));
Q=real(H(1,:)).*real(H(2,:)).*real(H(3,:))-...
    real(H(1,:)).*W3-real(H(2,:)).*W2-real(H(3,:)).*W1+...
    2*(real(H(4,:)).*(real(H(5,:)).*real(H(6,:))+imag(H(5,:).*imag(H(6,:)))))+...
    2*(imag(H(4,:)).*(imag(H(5,:)).*real(H(6,:))-real(H(5,:).*imag(H(6,:)))));
P=real(H(1,:)).*real(H(2,:))+real(H(1,:)).*real(H(3,:))+...
    real(H(2,:)).*real(H(3,:)) -W1-W2-W3;
ABS_P=abs(P);
C=(3./ABS_P).^1.5.*Q/2;
acosC=acos(C); 
P1=sqrt(4*ABS_P/3);
w=[P1.*cos(acosC/3)+sh;P1.*cos((acosC+4*pi)/3)+sh;P1.*cos((acosC+2*pi)/3)+sh];

function Freq=local_exact(Sys, phi, theta, Field)
Bo = Field*1e-3; % T
nu_i = nmagn/planck*Bo*Sys.gn(1, 1); % Hz

R = erot(safeget(Sys, 'Apa', [0,0,0]));
A = R * diag(Sys.A*1e6)*R.';

Ix = sop(Sys.I, 'x');
Iy = sop(Sys.I, 'y');
Iz = sop(Sys.I, 'z');

lx = sin(theta)*cos(phi);
ly = sin(theta)*sin(phi);
lz = cos(theta);

Gx=lx*Sys.g(1);
Gy=ly*Sys.g(2);
Gz=lz*Sys.g(3);

Zeeman_A = -nu_i*(lx*Ix + ly*Iy + lz*Iz);


geff=sqrt(Gx.*Gx+Gy.*Gy+Gz.*Gz);

HFI_A = (Ix*A(1, 1) + Iy*A(2, 1) + Iz*A(3, 1))*Gx./geff*1/2+...
        (Ix*A(1, 2) + Iy*A(2, 2) + Iz*A(3, 2))*Gy./geff*1/2++...
        (Ix*A(1, 3) + Iy*A(2, 3) + Iz*A(3, 3))*Gz./geff*1/2;
%       Quadrupole coupling 
Quad = 0;
if isfield(Sys, 'Q');
    if sum(abs(Sys.Q)) ~= 0

        Q = diag(Sys.Q)*1e6; %[Hz]
        if isfield(Sys, 'Qpa'),
            RMQ = erot(Sys.Qpa);
        else
            RMQ = diag(ones(3, 1));
        end
        Q = RMQ*Q*RMQ.';
        Quad = Q(1,1)*Ix*Ix + Q(1, 2)*Ix*Iy + Q(1, 3)*Ix*Iz +...
               Q(2,1)*Iy*Ix + Q(2, 2)*Iy*Iy + Q(2, 3)*Iy*Iz +...
               Q(3,1)*Iz*Ix + Q(3, 2)*Iz*Iy + Q(3, 3)*Iz*Iz;
    end
end
    
Ham_A = Zeeman_A - HFI_A + Quad;
Ham_B = Zeeman_A + HFI_A + Quad;

[V, D] = eig(Ham_A);
dd = diag(D);
oness = dd*0+1;
Ea = dd(:, oness)  - dd(:, oness).';

[V, D] = eig(Ham_B);
dd = diag(D);
Eb = dd(:, oness)  - dd(:, oness).';
Freq = [];
nI = (2*Sys.I+1);
for ii = 1:(nI-1)
    for jj = (ii+1):nI
        Freq(end+1, 1) = Ea(ii, jj);
        Freq(end, 2) = Eb(ii, jj);
    end
end
