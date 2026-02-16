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

function stencil(f, dim, i, j, I)

    if dim == "x"
        println("f: ", size(f))
        fmmm_a = leftx(f, I...)
        println("fmmm_a: ", size(fmmm_a))
        fmmm_b = leftx(fmmm_a, I...)
        
        "almost there with this ^
        
        instead i need to define bLx = Val(boundary[1]) and set the i,j = Tuple(I) 
        then do the illl, ill, il, ir, irr, irrr with illl = left_index(i,3,nx,bLx)
        u1 = u[illl,j] where u1 is the fmmm since u is the field f"
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
