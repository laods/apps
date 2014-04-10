function generateRockFile(points, output, swir, swor, krw0)
% Generate rock files with
%   k_rw = S^2
%   k_ro = (1-S)^2
%   J    = 1/S^2 - 1/(1-S)^2
% (possibly scaled if swir and swor given

if nargin < 3
  swir = 0;
  swor = 0;
  krw0 = 1;
end
if nargin < 2
  output = 0;
end

S   = swir:(1-swor-swir)/(points-1):1-swor;
Swn = (S-swir) ./ (1-swir-swor);
krw = krw0*Swn.^2;
kro = (1-Swn).^2;
J   = ( 1 ./ (Swn.^2) - 1 ./ ((1-Swn).^2) );
J(1) = 10*J(2);
J(end) = 10*J(end-1);

%[S' Swn' krw' kro' J']

if output ~= 0
  fid = fopen(output, 'w');
  if fid == -1
    error(sprintf('Could not open %s', output));
  end
  fprintf('Writing to %s...\n', output);
  fprintf(fid, '-- Rock file generated by %s\n', mfilename('fullpath'));
  fprintf(fid, '-- Sw\tkrw\t\tkro\t\tJ\n');
  for i=1:points
    fprintf(fid, '%1.4f\t%1.3e\t%1.3e\t%1.3e\n', S(i), krw(i), kro(i), J(i));
  end
  fclose(fid);
end

close all
plot(S, J, 'g');
axis([0 1 -100 100]);
figure
plot(S, krw, 'b');
hold on
plot(S, kro, 'r');
hold off

end








