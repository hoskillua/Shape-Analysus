function runMVU(k, n)

if nargin < 2
    n = 1000;
end

[x, theta] = swissroll(n);
y = mvu(x, k);
figure; scatter3(y(1, :), y(2, :), y(3, :), [], theta, '.'); view(3); axis image vis3d off;

function [x, theta] = swissroll(n)
    theta = 1.5*pi + 3 * pi * rand(1, n);
    z = 10 * rand(1, n);
    x = [theta .* [cos(theta); sin(theta)]; z];
end

end