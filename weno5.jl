# WENO5 FUNCTION

struct dim{D} end

function polys(fmm,fm,fc,fp,fpp)

    # Prefactors
    a1 = 2f0
    a2 = 5f0
    a3 = 13f0/12f0
    a4 = 1f0/4f0

    # p1 = (2*fmm .- 7*fm .+ 11*fc )/6 # polynomials
    p1 = (a1*fmm .- 7f0*fm .+ 11f0*fc )*inv(6f0) # polynomials
    p2 = (-fm .+ a2*fc .+ a1*fp )*inv(6f0)
    p3 = (a1*fc .+ a2*fp .- fpp)*inv(6f0)

    b1 = a3*(fmm .- a1*fm .+ fc ).^2f0 .+ a4*(fmm .- 4f0*fm .+ 3f0*fc ).^2f0 # smoothenss
    b2 = a3*(fm  .- a1*fc .+ fp ).^2f0 .+ a4*(fm .- fp ).^2f0
    b3 = a3*(fc  .- a1*fp .+ fpp).^2f0 .+ a4*(3f0*fc .- 4f0*fp .+ fpp).^2f0

    g = (0.1f0, 0.6f0, 0.3f0)
    eps = 1f-6
    wp1 = g[1] ./ (b1 .+ eps)^2f0
    wp2 = g[2] ./ (b2 .+ eps)^2f0
    wp3 = g[3] ./ (b3 .+ eps)^2f0

    flux_face = (wp1 .* p1 .+ wp2 .* p2 .+ wp3 .* p3) ./ (wp1 .+ wp2 .+ wp3)
    
    return flux_face
end


function stencil(f, ::dim{0}, I, nx, ny) # x direction

    i, j = Tuple(I)
        
    boundary_x = (1,1) # left, right

    left = Val(boundary_x[1])
    right = Val(boundary_x[2])
    

    illl = minus_index(i,3,nx, left)
    ill  = minus_index(i,2,nx, left)
    il   = minus_index(i,1,nx, left)
    ir   = plus_index(i,1,nx, right)
    irr  = plus_index(i,2,nx, right)
    irrr = plus_index(i,3,nx, right)

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

    return fppos, fpneg, fmpos, fmneg

end

function stencil(f, ::dim{1}, I, nx, ny) # y direction

    i, j = Tuple(I)
        
    boundary_y = (0,0) # bottom, top

    bottom = Val(boundary_y[1])
    top = Val(boundary_y[2])

    jbbb = minus_index(j,3,ny, bottom)
    jbb  = minus_index(j,2,ny, bottom)
    jb   = minus_index(j,1,ny, bottom)
    jt   = plus_index(j,1,ny, top)
    jtt  = plus_index(j,2,ny, top)
    jttt = plus_index(j,3,ny, top)

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

    return fppos, fpneg, fmpos, fmneg

end

function face_velocity(V, ::dim{0}, i, j)

    vxm = V.x[i-1,j]
    vxp = V.x[i+1,j]

    vmpos = 0.5f0 * (vxm + abs(vxm))
    vmneg = 0.5f0 * (vxm - abs(vxm))
    vppos = 0.5f0 * (vxp + abs(vxp))
    vpneg = 0.5f0 * (vxp - abs(vxp))

    return vmpos, vmneg, vppos, vpneg

end

function face_velocity(V, ::dim{1}, i, j)

    vym = V.y[i,j-1]
    vyp = V.y[i,j+1]

    vmpos = 0.5f0 * (vym + abs(vym))
    vmneg = 0.5f0 * (vym - abs(vym))
    vppos = 0.5f0 * (vyp + abs(vyp)) 
    vpneg = 0.5f0 * (vyp - abs(vyp))

    return vmpos, vmneg, vppos, vpneg

end

function advection_fct!(q_adv, f, V, dx, nx, ny, g::StructuredGrid, I, O)
    
    i, j = Tuple(I)

    vxmpos, vxmneg, vxppos, vxpneg = face_velocity(V, dim{0}(), i, j) # get the face velocities for the x direction - will need to be adapted to work with the grid structure and kernel launches
    vympos, vymneg, vyppos, vypneg = face_velocity(V, dim{1}(), i, j) # get the face velocities for the y direction - will need to be adapted to work with the grid structure and kernel launches

    fxppos, fxpneg, fxmpos, fxmneg = stencil(f, dim{0}(), I, nx, ny) # get the stencils for the x direction - will need to be adapted to work with the grid structure and kernel launches
    fyppos, fypneg, fympos, fymneg = stencil(f, dim{1}(), I, nx, ny)
    

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