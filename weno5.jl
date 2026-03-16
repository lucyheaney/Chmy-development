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

function advection_fct!(q_adv, f, V, dx, nx, ny, g::StructuredGrid, I, O)
    
    i, j = Tuple(I)

    vxmpos, vxmneg, vxppos, vxpneg = face_velocity(V, "x", i, j) # get the face velocities for the x direction - will need to be adapted to work with the grid structure and kernel launches vypmpos, vypmneg, vyppos, vypneg = face_velocity(V, "y") fxppos, fxpneg, fxmpos, fxmneg = stencil(f, "x", I) # get the stencils for the x direction - will need to be adapted to work with the grid structure and kernel launches fyppos, fypneg, fympos, fymneg = stencil(f, "y", I) # compute the WENO5 fluxes for x and y directions using the stencils and face velocities wp1 = g[1] ./ (b1 .+ eps)^2 wp2 = g[2] ./ (b2 .+ eps)^2
    vympos, vymneg, vyppos, vypneg = face_velocity(V, "y", i, j) # get the face velocities for the y direction - will need to be adapted to work with the grid structure and kernel launches # compute the WENO5 fluxes for x and y directions using the stencils and face velocities F_xphalf = vxmpos * fxppos + vxmneg * fxpneg F_xmhalf = vxmpos * fxmpos + vxmneg * fxmneg q_adv.x[I...] = (F_xphalf - F_xmhalf) / dx F_yphalf = vympos * fyppos + vymneg * fypneg F_ymhalf = vympos * fympos + vymneg * fymneg q_adv.y[I...] = (F_yphalf - F_ymhalf) / dx end

    fxppos, fxpneg, fxmpos, fxmneg = stencil(f, "x", I, nx, ny) # get the stencils for the x direction - will need to be adapted to work with the grid structure and kernel launches
    fyppos, fypneg, fympos, fymneg = stencil(f, "y", I, nx, ny)
    

    qxp = vxppos .* fxppos + vxpneg .* fxpneg
    qxm = vxmpos .* fxmpos + vxmneg .* fxmneg
    q_adv.x[I...] = (qxp - qxm) / dx #???? need to revisit

    qyp = vyppos .* fyppos + vypneg .* fypneg
    qym = vympos .* fympos + vymneg .* fymneg
    q_adv.y[I...] = (qyp - qym) / dx

    #adv[I...] = (qxp - qxm) / dx + (qyp - qym) / dx # change this to correct q_advs
    return q_adv.x[I...] + q_adv.y[I...] #rate of change of scalar field due to advection
    #println(size(adv), ": adv")

end

function diffusion_fct!(q_diff, f, κ, g::StructuredGrid, I, O)
    
    i, j = Tuple(I)

    q_diff.x[I...] = -κ * ∂x(f, g, I...)
    q_diff.y[I...] = -κ * ∂y(f, g, I...)

    return -(q_diff.x[I...] + q_diff.y[I...]) # this is the divergence of the diffusive flux - need to check the signs and make sure its correct for the update of f in the update_thermal kernel

    #println(size(diff), ": diff")
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


@kernel inbounds = true function update_old_stokes!(τ,τ_old, O) # 
    I = @index(Global, NTuple) # this sets up the indexes depending on 2D 3D etc
    I = I + O # O is an offset tuple to shift the index
    τ_old.xx[I...] = τ.xx[I...] # storing those values
    τ_old.yy[I...] = τ.yy[I...]
    τ_old.xy[I...] = τ.xy[I...]
end

@kernel inbounds = true function update_old_fields!(f, f_old, O) # 
    I = @index(Global, NTuple) # this sets up the indexes depending on 2D 3D etc
    I = I + O # O is an offset tuple to shift the index
    f_old[I...] = f[I...] # I... means split the tuple into individual indexes 
end