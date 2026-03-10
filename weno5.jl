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

function stencil(f, dim, I, nx, ny)

    i, j = Tuple(I)
        
    boundary = (1,1,0,0)

    left = Val(boundary[1])
    right = Val(boundary[2])
    bottom = Val(boundary[3])
    top = Val(boundary[4])

    illl = minus_index(i,3,nx, left)
    ill  = minus_index(i,2,nx, left)
    il   = minus_index(i,1,nx, left)
    ir   = plus_index(i,1,nx, right)
    irr  = plus_index(i,2,nx, right)
    irrr = plus_index(i,3,nx, right)

    jbbb = minus_index(j,3,ny, bottom)
    jbb  = minus_index(j,2,ny, bottom)
    jb   = minus_index(j,1,ny, bottom)
    jt   = plus_index(j,1,ny, top)
    jtt  = plus_index(j,2,ny, top)
    jttt = plus_index(j,3,ny, top)


    if dim == "x"

        fmmm = f[illl, j]
        fmm  = f[ill, j]
        fm   = f[il, j]
        fc   = f[i,   j]
        fp   = f[ir, j]
        fpp  = f[irr, j]
        fppp = f[irrr, j]

        fppos = polys(fmmm, fmm, fm, fc, fp) # need to define these stencils for the weno5 scheme - will need to be adapted to work with the grid structure and kernel launches
        fpneg = polys(fppp, fpp, fp, fc, fm)
        fmpos = polys(fmmm, fmm, fm, fc, fp)
        fmneg = polys(fpp, fp, fc, fm, fmm)

    elseif dim == "y"
        fmmm = f[i, jbbb]
        fmm  = f[i, jbb]
        fm   = f[i, jb]
        fc   = f[i, j]
        fp   = f[i, jt]
        fpp  = f[i, jtt]
        fppp = f[i, jttt]

        fppos = polys(fmmm, fmm, fm, fc, fp) # need to define these stencils for the weno5 scheme - will need to be adapted to work with the grid structure and kernel launches
        fpneg = polys(fppp, fpp, fp, fc, fm)
        fmpos = polys(fmmm, fmm, fm, fc, fp)
        fmneg = polys(fpp, fp, fc, fm, fmm)
    end

    return fppos, fpneg, fmpos, fmneg

end

function face_velocity(V, dim, i, j)

    if dim == "x"

        vx = V.x[i,j]
        vxm = V.x[i-1,j]
        vxp = V.x[i+1,j]

        vmpos = 0.5 * (vxm + abs(vxm))
        vmneg = 0.5 * (vxm - abs(vxm))

        vppos = 0.5 * (vxp + abs(vxp))
        vpneg = 0.5 * (vxp - abs(vxp))

    elseif dim == "y"

        vy = V.y[i,j]
        vym = V.y[i,j-1]
        vyp = V.y[i,j+1]

        vmpos = 0.5 * (vym + abs(vym))
        vmneg = 0.5 * (vym - abs(vym))

        vppos = 0.5 * (vyp + abs(vyp)) 
        vpneg = 0.5 * (vyp - abs(vyp))
    end

    return vmpos, vmneg, vppos, vpneg

end




# neumann (mirrored boundary value)
function minus_index(i, d, nx, ::Val{0})
    return max(i-d,1)
end

# periodic (wrap around)
function minus_index(i, d, nx, ::Val{1})
    return mod1(i-d, nx)
end

function plus_index(i, d, nx, ::Val{0})
    return min(i+d, nx)
end

function plus_index(i, d, nx, ::Val{1})
    return mod1(i+d, nx)
end


