% make discrete time PID controller from matrices
% supplied by Jo Lister.
% A S Sharma 2001

load dina_13333_500
load ps

dt = 1e-4;
method = 'tustin';

% generate the P,I,D and give them the correct dimensions

% integrator
[a,b,c,d]=tf2ss([0 1],[1 0]);
sysdi = c2d(ss(a,b,c,d),dt,method);
sysdii = sysdi;
for i=2:ncp sysdii = append(sysdii,sysdi); end

% proportional
[a,b,c,d]=tf2ss([0 1],[0 1]);
sysdp = c2d(ss(a,b,c,d),dt,method);
sysdpp = sysdp;
for i=2:ncp sysdpp = append(sysdpp,sysdp); end

% derivative
tau=1e-5;
[a,b,c,d]=tf2ss([1 0],[tau 1]);
sysdd = c2d(ss(a,b,c,d),dt,method);
sysddd = sysdd;
for i=2:ncp sysddd = append(sysddd,sysdd); end

mmat = [M_mat(1:18,1:18);zeros(1,18)];
if(shot==19288);
 mmat = M_mat(1:19,1:22);     % for 19288
end
% modification so that 13333 still runs / JBL 20/3/3
g_i2 = g_i(1:size(mmat,2),:);
g_p2 = g_p(1:size(mmat,2),:);
g_d2 = g_d(1:size(mmat,2),:);

% Modification to lower zIp proportional gain / AS 24/3/03
%g_p2(:,4) = g_p2(:,4)/1.2;
%g_d2(:,4) = g_d2(:,4)*1.5;

sysd = parallel(mmat*g_i2*sysdii,mmat*g_p2*sysdpp);
sysd = parallel(sysd,mmat*g_d2*sysddd);

Ak = sysd.a;
Bk = sysd.b;
Ck = sysd.c;
Dk = sysd.d;

save PID_discrete *k