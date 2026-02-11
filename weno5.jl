# WENO5 FUNCTION

function polys(fmm,fm,fc,fp,fpp)
    p1 = (2*fmm - 7*fm + 11*fc )/6
    p2 = ( -fm  + 5*fc +  2*fp )/6
    p3 = (2*fc  + 5*fp -    fpp)/6

    b1 = 13/12*(fmm - 2*fm + fc ).^2 + 1/4*(  fmm - 4*fm + 3*fc ).^2
    b2 = 13/12*(fm  - 2*fc + fp ).^2 + 1/4*(  fm  -          fp ).^2
    b3 = 13/12*(fc  - 2*fp + fpp).^2 + 1/4*(3*fc  - 4*fp +   fpp).^2

    g = (0.1, 0.6, 0.3)
    eps = 1e-6
    wp = g[1] ./ (b1 + eps)^2
    wp = g[2] ./ (b2 + eps)^2
    wp = g[3] ./ (b3 + eps)^2

    flux_face = (wp1 .* p1 + wp2 .* p2 + wp3 .* p3) ./ (wp1 + wp2 + wp3)
end

# set the stencil here where we need to define fmm , fmm, fm, f, fp, fpp, fppp

fppos = polys(fmm, fm, f, fp, fpp)
fpneg = polys(fppp, fpp, fp, f, fm)

fmpos = polys(fmmm, fmm, fm, f, fp)
fmneg = polys(fpp, fp, f, fm, fmm)

# remember can do .@
fmm = 2
fm = 2
fc = 2
p1 = (2*fmm - 7*fm + 11*fc )/6