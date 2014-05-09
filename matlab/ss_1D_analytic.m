function ss_1D_analytic(u, uw, s_in)

% Assume swit = 0 and swor = 1

% Define rock boundaries
N = 2;   % Nr of domains (rocks)
n = 100; % Nr of discrete points in each domain

% Parameters
mu = [1e-3 0.1]; % [water, oil]
mD = 9.869e-16;
k  = 1*mD;
pc_entry = [1e2, 1e3];

% Fluid functions
krw    = @(s)s;
kro    = @(s)(1-s);
pc     = @(s,r)pcfun(pc_entry, r, s);
pcinv  = @(p,r)pcinvfun(pc, r, 1000, p);
dpc    = @(s,r)dpcfun(pc_entry, r, s);

% Derived functions in Dale et. al. (eq. 5 and 6)
Phi  = @(s)Phifun(krw, kro, mu, s);
Psi  = @(s,r)Psifun(krw, kro, dpc, mu, k, u, s, r);

Delta_12     = @(s)(pcinv(pc(s, 1), 2));
integrand{1} = @(s)integrandfun(Psi, Phi, uw/u, 1, s);
integrand{2} = @(s)integrandfun(Psi, Phi, uw/u, 2, s);

sat_vec = 0:0.001:1;
plot(sat_vec, integrand{1}(sat_vec), '-r');

%pc_vec  = -1e8:1e6:1e8;

%plot(sat_vec, Phi(sat_vec), '-r');
%plot(sat_vec, Psi(sat_vec,1), '-b');
%plot(sat_vec, Psi(sat_vec,2), '-g');
%plot(sat_vec, pc(sat_vec,1));
%hold on
%plot(sat_vec, pc(sat_vec,2));
%axis([0,1,-1e5,1e5])
%plot(sat_vec, dpc(sat_vec,1));
%hold off

% Main loop
f  = uw/u;
Sa = calcSa(Phi, f);
Sl = s_in; % Saturation at left boundary (xl)
S  = zeros(N,n); % Solution
for i = 1:N
    % Homogeneous flow
    if Sl == Sa
        S(i) = ones(1,n)*Sa;
        continue
    elseif (Sl == 0) && (Psi(0,i) == -Inf)
        S(i) = zeros(1,n);
        continue    
    elseif (Sl == 1) && (Psi(0,i) == -Inf)
        S(i) = ones(1,n);
        continue
    end
    x = zeros(1,n);
    s = 0:(1/(n-1)):1;
    for j = 1:n
        x(j) = integral(integrand{i}, Sl, s(j));
    end
    xmax = max(x);
    S(i,:) = interp1(x, s, (i-1):(1/(n-1)):i);
    Sl   = Delta_12(S(i,n));
end

plot([S(1,:), S(2,:)])

end


function pc = pcfun(pc_entry, r, s)
assert(r==1 || r==2);
assert(numel(pc_entry) == 2);
pc = 1./s.^2 - 1./(1-s).^2;
pc = pc*pc_entry(r);
%pc(s==0) =  0.1/eps;
%pc(s==1) = -0.1/eps;
end

function dpc = dpcfun(pc_entry, r, s)
assert(r==1 || r==2);
assert(numel(pc_entry) == 2);
dpc = 1./s.^3 + 1./(1-s).^3;
dpc = -2*pc_entry(r)*dpc;
%dpc(s==0) = -0.1/eps;
%dpc(s==1) = -0.1/eps;
end

function sat = pcinvfun(pc_fun, r, n, pc)
% Invert function handle pc_fun (by linear interpolation with n points)
assert(r==1 || r==2);
s_vec  = 0:(1/n):1;
pc_vec = pc_fun(s_vec, r);
pc_vec(pc_vec==Inf)  =  0.1/eps;
pc_vec(pc_vec==-Inf) = -0.1/eps;
sat = interp1(pc_vec, s_vec, pc);
end

function Phi = Phifun(krw, kro, mu, s)
Phi = mu(2)*krw(s) ./ (mu(2)*krw(s) + mu(1)*kro(s));
end

function Psi = Psifun(krw, kro, dpc, mu, k, u, s, r)
assert(r==1 || r==2);
Psi = k*krw(s).*kro(s).*dpc(s, r) ./ (u * (mu(2)*krw(s) + mu(1)*kro(s)) );
Psi(isnan(Psi)) = -inf;
end

function val = integrandfun(Psi, Phi, f, r, s)
val = Psi(s,r) ./ (f - Phi(s));
val(val==Inf)  =  1e5;
val(val==-Inf) = -1e5;
end

function Sa = calcSa(Phi, f)
max_it = 1000;
tol = 1e-5;
Sa = 0.5;
i = 0;
I = [0, 1];
while ( abs(Phi(Sa) - f) > tol ) && ( i < max_it )
    if Phi(Sa) > f
        I  = [I(1), Sa];
    else
        I  = [Sa, I(2)];
    end
    Sa = sum(I)/2;
end
if ( abs(Phi(Sa) - f) > tol )
    error('Failed to calculate Sa');
end

end

