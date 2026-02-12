# WENO5 FUNCTION

function polys(fmm,fm,fc,fp,fpp)
    p1 = (2*fmm .- 7*fm .+ 11*fc )/6 # polynomials
    p2 = (-fm .+ 5*fc .+ 2*fp )/6
    p3 = (2*fc .+ 5*fp .- fpp)/6

    b1 = 13/12*(fmm .- 2*fm .+ fc ).^2 .+ 1/4*(fmm .- 4*fm .+ 3*fc ).^2 # smoothenss
    b2 = 13/12*(fm  .- 2*fc .+ fp ).^2 .+ 1/4*(fm .- fp ).^2
    b3 = 13/12*(fc  .- 2*fp .+ fpp).^2 .+ 1/4*(3*fc  .- 4*fp .+   fpp).^2

    g = (0.1, 0.6, 0.3)
    eps = 1e-6
    wp1 = g[1] ./ (b1 .+ eps)^2
    wp2 = g[2] ./ (b2 .+ eps)^2
    wp3 = g[3] ./ (b3 .+ eps)^2

    flux_face = (wp1 .* p1 .+ wp2 .* p2 .+ wp3 .* p3) ./ (wp1 .+ wp2 .+ wp3)
    
    return flux_face
end

function stencil(f, dim, I)

    i, j = I

    if dim == "x"
        fmmm = f[i-3, j]
        fmm  = f[i-2, j]
        fm   = f[i-1, j]
        fc   = f[i,   j]
        fp   = f[i+1, j]
        fpp  = f[i+2, j]
        fppp = f[i+3, j]

        fppos = polys(fmmm, fmm, fm, fc, fp) # need to define these stencils for the weno5 scheme - will need to be adapted to work with the grid structure and kernel launches
        fpneg = polys(fppp, fpp, fp, fc, fm)
        fmpos = polys(fmmm, fmm, fm, fc, fp)
        fmneg = polys(fpp, fp, fc, fm, fmm)

    elseif dim == "y"
        fmmm = f[i, j-3]
        fmm  = f[i, j-2]
        fm   = f[i, j-1]
        fc   = f[i, j]
        fp   = f[i, j+1]
        fpp  = f[i, j+2]
        fppp = f[i, j+3]

        fppos = polys(fmmm, fmm, fm, fc, fp) # need to define these stencils for the weno5 scheme - will need to be adapted to work with the grid structure and kernel launches
        fpneg = polys(fppp, fpp, fp, fc, fm)
        fmpos = polys(fmmm, fmm, fm, fc, fp)
        fmneg = polys(fpp, fp, fc, fm, fmm)
    end

    return fppos, fpneg, fmpos, fmneg

end

function face_velocity(V, dim)

    if dim == "x"
        vx = V.x[i,j]

        v_pos = max(vx, 0f0)
        v_neg = min(vx, 0f0)

    elseif dim == "y"
        vy = V.y[i,j]

        v_pos = max(vy, 0f0) 
        v_neg = min(vy, 0f0)
    end

    return v_pos, v_neg

end
