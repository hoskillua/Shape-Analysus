function elasticRods(bendModulus, twistModulus, totalTwist)

arguments
    bendModulus (1, 1) double = 1
    twistModulus (1, 1) double = 1
    totalTwist (1, 1) double = pi
end

nSamples = 100;
dt = 0.001;
curveFunction = @(t) [cos(2 * pi * t); sin(2 * pi * t); 0.3 * sin(4 * pi * t)];
curveData.verts = curveFunction(linspace(0, 1, nSamples + 1));
curveData.verts = curveData.verts(:, 1:nSamples);
curveData.totalTwist = totalTwist;
curveData.velocities = zeros(3, nSamples);
curveData = prepare(curveData);

% Set up Bishop frame
u0 = cross(curveData.tangentsR(:, 1), [0; 0; 1]);
u0 = u0 ./ sqrt(sum(u0.^2));
v0 = cross(curveData.tangentsR(:, 1), u0);
v0 = v0 ./ sqrt(sum(v0.^2));
curveData.bishopFrame = propagateBishopFrame(curveData, [u0 v0]);
curveData = updateMaterialFrame(curveData);

links = [(1:nSamples).' circshift((1:nSamples).', -1)];
DCii = repmat((1:nSamples).', 1, 6);
DCjj = mod(3 * ((1:nSamples).' - 1) + (0:5), 3 * nSamples) + 1;

figure;
patch('Faces', links, 'Vertices', curveData.verts.', 'LineWidth', 3, 'EdgeColor', 'black'); hold on;
quiver3(curveData.midpointsR(1, :), curveData.midpointsR(2, :), curveData.midpointsR(3, :),...
        curveData.u(1, 1:end-1), curveData.u(2, 1:end-1), curveData.u(3, 1:end-1),...
        'color','red','linewidth',1);
view(3); axis image vis3d manual off;
ax = gca;
ax.Clipping = 'off';

for i = 1:10000
    bendForce = computeBendForce(curveData);
    twistForce = computeTwistForce(curveData);
    totalForce = bendModulus * bendForce + twistModulus * twistForce;
    totalAcceleration = totalForce ./ curveData.dualLengths;

    [verts, velocities] = symplecticEuler(curveData.verts, curveData.velocities, totalAcceleration, dt);
    [verts, velocities] = fastProjection(curveData.verts, verts, curveData.dualLengths, curveData.edgeLengthsR, dt);
    newCurveData.verts = verts;
    newCurveData.velocities = velocities;
    newCurveData = prepare(newCurveData);
    newCurveData = updateTwist(newCurveData, curveData);
    newCurveData = updateMaterialFrame(newCurveData);
    curveData = newCurveData;
    
    if mod(i - 1, 40) == 0
        pause(0.001);
        cla;
        patch('Faces', links, 'Vertices', curveData.verts.', 'LineWidth', 3, 'EdgeColor', 'black'); hold on;
        quiver3(curveData.midpointsR(1, :), curveData.midpointsR(2, :), curveData.midpointsR(3, :),...
                curveData.u(1, 1:end-1), curveData.u(2, 1:end-1), curveData.u(3, 1:end-1),...
                'color','red','linewidth',1);
    end
end

function curveData = prepare(curveData)
    % Basic closed curve
    curveData.edgesR = circshift(curveData.verts, -1, 2) - curveData.verts;
    curveData.edgesL = circshift(curveData.edgesR, 1, 2);
    curveData.midpointsR = 0.5 * (circshift(curveData.verts, -1, 2) + curveData.verts);
    curveData.edgeLengthsR = sqrt(sum(curveData.edgesR.^2, 1));
    curveData.edgeLengthsL = circshift(curveData.edgeLengthsR, 1, 2);
    curveData.dualLengths = 0.5 * (curveData.edgeLengthsR + curveData.edgeLengthsL);
    curveData.totalLength = sum(curveData.dualLengths);
    
    curveData.tangentsR = curveData.edgesR ./ curveData.edgeLengthsR;
    curveData.tangentsL = circshift(curveData.tangentsR, 1, 2);
    curveData.curvatureBinormalsDenom = curveData.edgeLengthsL .* curveData.edgeLengthsR ...
                                      + dot(curveData.edgesL, curveData.edgesR, 1);
    curveData.curvatureBinormals = 2 * cross(curveData.edgesL, curveData.edgesR, 1) ...
                                 ./ curveData.curvatureBinormalsDenom;
	curveData.curvatureSquared = sum(curveData.curvatureBinormals.^2, 1);
end

function [totalEnergy, bendingEnergy, twistEnergy, kineticEnergy] = computeEnergy(curveData)
    bendingEnergy = bendModulus * sum(sum(curveData.curvatureBinormals.^2, 1) ./ (2 * curveData.dualLengths));
    twistEnergy = twistModulus * curveData.totalTwist.^2 / (2 * curveData.totalLength);
    kineticEnergy = 0.5 * sum(curveData.dualLengths .* sum(curveData.velocities.^2, 1));
    totalEnergy = bendingEnergy + twistEnergy + kineticEnergy;
end

function bishopFrame = propagateBishopFrame(curveData, startFrame)
    %%% PROBLEM 3(a) - YOUR CODE HERE TO PARALLEL TRANSPORT AROUND THE CURVE
    parallelTransport = zeros(3, 3, nSamples);
    %%% END HOMEWORK PROBLEM
    
    bishopFrame = zeros(3, 2, nSamples + 1);
    bishopFrame(:, :, 1) = startFrame;
    for j=2:(nSamples + 1)
        ptIdx = mod(j - 1, nSamples) + 1;
        bishopFrame(:, :, j) = parallelTransport(:, :, ptIdx) * bishopFrame(:, :, j - 1);
    end
end

function newCurveData = updateTwist(newCurveData, oldCurveData)
    timePtAxis = cross(oldCurveData.tangentsR(:, 1), newCurveData.tangentsR(:, 1));
    timePtAngle = atan2(sqrt(sum(timePtAxis.^2, 1)), dot(oldCurveData.tangentsR(:, 1), newCurveData.tangentsR(:, 1)));
    timePt = axang2rotm([timePtAxis ./ sqrt(sum(timePtAxis.^2, 1)); timePtAngle].');
    
    newCurveData.bishopFrame = propagateBishopFrame(newCurveData, timePt * oldCurveData.bishopFrame(:, :, 1));
    
    holonomy = (timePt * oldCurveData.bishopFrame(:, :, nSamples + 1))' * newCurveData.bishopFrame(:, :, nSamples + 1);
    holonomyAngle = atan2(holonomy(2, 1), holonomy(1, 1));
    newCurveData.totalTwist = oldCurveData.totalTwist - holonomyAngle;
end

function curveData = updateMaterialFrame(curveData)
    cumTwist = curveData.totalTwist * cumsum([0 circshift(curveData.dualLengths, -1)]) ./ curveData.totalLength;
    curveData.u = cos(cumTwist) .* squeeze(curveData.bishopFrame(:, 1, :)) + ...
                  sin(cumTwist) .* squeeze(curveData.bishopFrame(:, 2, :));
    curveData.v = cos(cumTwist) .* squeeze(curveData.bishopFrame(:, 2, :)) - ...
                  sin(cumTwist) .* squeeze(curveData.bishopFrame(:, 1, :));
end


%%% PROBLEM 3(c) Part I - YOUR CODE HERE
function bendForce = computeBendForce(curveData)
    bendForce = zeros(size(curveData.verts));
end
%%% END HOMEWORK PROBLEM

%%% PROBLEM 3(c) Part II - YOUR CODE HERE
function twistForce = computeTwistForce(curveData)
    twistForce = zeros(size(curveData.verts));
end
%%% END HOMEWORK PROBLEM

function [verts, velocities] = symplecticEuler(verts0, velocities0, totalAcceleration, dt)
    % Update velocities and positions via Symplectic Euler Integration
    velocities = velocities0 + dt * totalAcceleration;
    verts = verts0 + dt * velocities;
end

function [verts, velocities] = fastProjection(verts0, verts, mass, lengthsR, dt)
    % Project positions to satisfy length constraints
    M = spdiags(repelem(reshape(mass, nSamples, 1), 3, 1), 0, 3 * nSamples, 3 * nSamples);
    edgesR = circshift(verts, -1, 2) - verts;
    constraint = (sum(edgesR.^2, 1) - lengthsR.^2).';
    while max(abs(constraint)) > 1e-10
        constraintGrad = 2 * sparse(DCii, DCjj, [-edgesR.' edgesR.'], nSamples, 3 * nSamples);
        MinvDC = M \ constraintGrad';
        DCMinvDC = constraintGrad * MinvDC;
        dLambda = DCMinvDC \ constraint;
        dx = -reshape(MinvDC * dLambda, 3, nSamples);
        verts = verts + dx;
        edgesR = circshift(verts, -1, 2) - verts;
        constraint = (sum(edgesR.^2, 1) - lengthsR.^2).';
    end
    velocities = (verts - verts0) / dt;
end

end