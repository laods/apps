function ss_1D_analytic(u, uw, s_in)

% Assume swit = 0 and swor = 1

% Define rock boundaries
x = 0:2;

% Parameters
mu = [1e-3 0.1]; % [water, oil]
mD = 9.869e-16;
k  = 1*mD;
pc_entry = [1e2, 1e4];

% Fluid functions
krw    = @(s)s;
kro    = @(s)(1-s);
pc     = @(s,r)pcfun(pc_entry, r, s);
pcinv  = @(pc,r)pcinvfun(pc, r, 100, pc);
dpc    = @(s,r)dpcfun(pc_entry, r, s);

% Derived functions in Dale et. al. (eq. 5 and 6)
Phi  = @(s)Phifun(krw, kro, mu, s);
Psi  = @(s,r)Psifun(krw, kro, dpc, mu, k, u, s, r);

sat_vec = 0:0.01:1;
pc_vec  = -1e8:1e6:1e8;


plot(sat_vec, Phi(sat_vec), '-r');
hold on
plot(sat_vec, Psi(sat_vec,1), '-b');
plot(sat_vec, Psi(sat_vec,2), '-g');
%plot(sat_vec, pc(sat_vec,1));
%plot(sat_vec, dpc(sat_vec,1));
hold off

% All values up till here are verified!

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