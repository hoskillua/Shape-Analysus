module ElasticRods

export elasticRods

using LinearAlgebra
using SparseArrays
using Rotations
using WGLMakie

function multicross(x, y)
    reduce(hcat, cross.(eachcol(x), eachcol(y)))
end

mutable struct CurveData
    verts :: Array{Float64}
    velocities :: Array{Float64}
    nv :: Int
    edgesR :: Array{Float64}
    edgesL :: Array{Float64}
    midpointsR :: Array{Float64}
    edgeLengthsR :: Array{Float64}
    edgeLengthsL :: Array{Float64}
    dualLengths :: Array{Float64}
    massMatrix :: Array{Float64}
    totalLength :: Float64
    tangentsR :: Array{Float64}
    tangentsL :: Array{Float64}
    curvatureBinormalsDenom :: Array{Float64}
    curvatureBinormals :: Array{Float64}
    curvatureSquared :: Array{Float64}
    bishopFrame :: Array{Float64}
    totalTwist :: Float64
    u :: Array{Float64}

    function CurveData(verts, velocities)
        curveData = new(verts, velocities)

        curveData.nv = size(verts, 2)
        curveData.edgesR = circshift(curveData.verts, (0, -1)) .- curveData.verts
        curveData.edgesL = circshift(curveData.edgesR, (0, 1))
        curveData.midpointsR = 0.5 .* (circshift(curveData.verts, (0, -1)) .+ curveData.verts)
        curveData.edgeLengthsR = sqrt.(sum(curveData.edgesR.^2, dims=1))
        curveData.edgeLengthsL = circshift(curveData.edgeLengthsR, (0, 1))
        curveData.dualLengths = 0.5 .* (curveData.edgeLengthsR .+ curveData.edgeLengthsL)
        curveData.totalLength = sum(curveData.dualLengths)
        curveData.massMatrix = kron(Diagonal(curveData.dualLengths[:]), Diagonal(ones(3)))
        
        curveData.tangentsR = curveData.edgesR ./ curveData.edgeLengthsR
        curveData.tangentsL = circshift(curveData.tangentsR, (0, 1))
        curveData.curvatureBinormalsDenom = curveData.edgeLengthsL .* curveData.edgeLengthsR .+
                                            sum(curveData.edgesL .* curveData.edgesR, dims=1)
        curveData.curvatureBinormals = 2 * multicross(curveData.edgesL, curveData.edgesR) ./
                                       curveData.curvatureBinormalsDenom
        curveData.curvatureSquared = sum(curveData.curvatureBinormals.^2, dims=1)
        return curveData
    end
end

function elasticRods(bendModulus = 1, twistModulus = 1, totalTwist = pi)
    nSamples = 100
    dt = 0.001
    curveFunction(t) = [cos.(2pi .* t); sin.(2pi .* t); 0.3 .* sin.(4pi .* t)]
    verts = curveFunction(range(0, stop=1, length=nSamples + 1)')[:, 1:nSamples]
    curveData = CurveData(verts, zeros(3, nSamples))

    # Set up Bishop frame
    u0 = cross(curveData.tangentsR[:, 1], [0;0;1])
    u0 = u0 ./ sqrt(sum(u0.^2))
    v0 = cross(curveData.tangentsR[:, 1], u0)
    v0 = v0 ./ sqrt(sum(v0.^2))
    curveData.bishopFrame = propagateBishopFrame(curveData, [u0 v0])
    curveData.totalTwist = totalTwist
    updateMaterialFrame!(curveData)

    links = [(1:nSamples) circshift((1:nSamples), -1)]
    DCii = repeat((1:nSamples), 1, 6)
    DCjj = mod.(3 .* ((1:nSamples) .- 1) .+ (0:5)', 3 .* nSamples) .+ 1

    # Observables to be plotted
    loop = Node(Point3.(eachcol([verts verts[:, 1]])))
    bases = Node(Point3.(eachcol(curveData.midpointsR)))
    vecs = Node(Point3.(eachcol(curveData.u[:, 1:end-1])))

    fig = Figure(resolution=(1000, 1000))
    ax = LScene(fig)
    lines!(loop, linewidth=30)
    arrows!(bases, vecs, linewidth=10, linecolor=:red, arrowcolor=:red, arrowsize=0.01, lengthscale=0.1)
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)

    for i = 1:10000
        bendForce = computeBendForce(curveData)
        twistForce = computeTwistForce(curveData)
        totalForce = bendModulus .* bendForce .+ twistModulus .* twistForce
        totalAcceleration = totalForce ./ curveData.dualLengths

        verts, velocities = symplecticEuler(curveData.verts, curveData.velocities, totalAcceleration, dt)
        verts, velocities = fastProjection(curveData.verts, verts, curveData.massMatrix, curveData.edgeLengthsR, dt, DCii, DCjj)
        newCurveData = CurveData(verts, velocities)
        updateTwist!(newCurveData, curveData)
        updateMaterialFrame!(newCurveData)
        curveData = newCurveData
        
        if mod(i - 1, 40) == 0
            sleep(0.001)
            loop[] = Point3.(eachcol([verts verts[:, 1]]))
            bases[] = Point3.(eachcol(curveData.midpointsR))
            vecs[] = Point3.(eachcol(curveData.u[:, 1:end-1]))
        end
    end
end

function computeEnergy(curveData, bendModulus, twistModulus)
    bendingEnergy = bendModulus .* sum(sum(curveData.curvatureBinormals.^2, dims=1) ./ (2 .* curveData.dualLengths))
    twistEnergy = twistModulus .* curveData.totalTwist.^2 ./ (2 * curveData.totalLength)
    kineticEnergy = 0.5 .* sum(curveData.dualLengths .* sum(curveData.velocities.^2, dims=1))
    totalEnergy = bendingEnergy + twistEnergy + kineticEnergy
    return totalEnergy, bendingEnergy, twistEnergy, kineticEnergy
end

function propagateBishopFrame(curveData, startFrame)
    ### PROBLEM 3(a) - YOUR CODE HERE
    parallelTransport = fill(one(RotMatrix{3, Float64}), curveData.nv)
    ### END HOMEWORK PROBLEM
    
    bishopFrame = zeros(3, 2, curveData.nv + 1)
    bishopFrame[:, :, 1] = startFrame
    for j = 2:(curveData.nv + 1)
        ptIdx = mod(j - 1, curveData.nv) + 1
        bishopFrame[:, :, j] = parallelTransport[ptIdx] * bishopFrame[:, :, j - 1]
    end
    return bishopFrame
end

function updateTwist!(newCurveData, oldCurveData)
    timePtAxis = cross(oldCurveData.tangentsR[:, 1], newCurveData.tangentsR[:, 1])
    timePtAngle = atan(sqrt(sum(timePtAxis.^2)),
                       sum(oldCurveData.tangentsR[:, 1] .* newCurveData.tangentsR[:, 1]))
    timePt = RotMatrix(AngleAxis(timePtAngle, timePtAxis[1], timePtAxis[2], timePtAxis[3]))

    newCurveData.bishopFrame = propagateBishopFrame(newCurveData, timePt * oldCurveData.bishopFrame[:, :, 1])

    holonomy = (timePt * oldCurveData.bishopFrame[:, :, end])' * newCurveData.bishopFrame[:, :, end]
    holonomyAngle = atan(holonomy[2, 1], holonomy[1, 1])
    newCurveData.totalTwist = oldCurveData.totalTwist - holonomyAngle
end

function updateMaterialFrame!(curveData)
    cumTwist = curveData.totalTwist * cumsum([0 circshift(curveData.dualLengths, -1)], dims=2) ./ curveData.totalLength
    curveData.u = cos.(cumTwist) .* curveData.bishopFrame[:, 1, :] .+ sin.(cumTwist) .* curveData.bishopFrame[:, 2, :]
end

### PROBLEM 3(c) Part I - YOUR CODE HERE
function computeBendForce(curveData)
    bendForce = zeros(size(curveData.verts));
    return bendForce
end
### END HOMEWORK PROBLEM

### PROBLEM 3(c) Part II - YOUR CODE HERE
function computeTwistForce(curveData)
    twistForce = zeros(size(curveData.verts));
    return twistForce
end
### END HOMEWORK PROBLEM

function symplecticEuler(verts0, velocities0, totalAcceleration, dt)
    # Update velocities and positions via Symplectic Euler Integration
    velocities = velocities0 + dt * totalAcceleration;
    verts = verts0 + dt * velocities;
    return verts, velocities
end

function fastProjection(verts0, verts, massMatrix, lengthsR, dt, DCii, DCjj)
    # Project positions to satisfy length constraints
    edgesR = circshift(verts, (0, -1)) - verts;
    constraint = (sum(edgesR.^2, dims=1) .- lengthsR.^2)';
    while maximum(abs.(constraint)) > 1e-10
        constraintGrad = 2 * sparse(DCii[:], DCjj[:], vec([-edgesR' edgesR']), length(constraint), length(verts));
        MinvDC = massMatrix \ constraintGrad';
        DCMinvDC = constraintGrad * MinvDC;
        dLambda = DCMinvDC \ constraint;
        dx = -reshape(MinvDC * dLambda, size(verts));
        verts = verts + dx;
        edgesR = circshift(verts, (0, -1)) .- verts;
        constraint = (sum(edgesR.^2, dims=1) .- lengthsR.^2)';
    end
    velocities = (verts .- verts0) ./ dt;
    return verts, velocities
end

end