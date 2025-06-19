module MapBField2Surface
using OffsetArrays

# XXX refactor: see this as a ODE, where the r-steps are done by the B-field,
# and add callbacks to not go into the sun?
# for field-line tracing, see https://book.magneticearth.org/geomag-obs-models/03a_magnetic-field-line-tracing

function mapb2s(B, ∇B, nr, nθ, nφ)

    # magnetic field at 1 R_sun
    B_s = zeros(0:nr, 0:nθ, 0:nφ, 2)

    r_min = zeros(3)
    r0 = zeros(3)
    r = zeros(3)
    b = zeros(3)

    t = 0
    dr = 0.01
    dt = 0.005

    radii = zeros(0:nr)
    radii .= 1.00001 + axes(radii, 1) * 0.01
    θs = zeros(0:nθ)
    θs .= axes(θs, 1) * DEG_TO_RAD
    φs = zeros(0:nφ)
    φs .= axes(φs, 1) * DEG_TO_RAD

    for (k, r_iter) in enumerate(radii),
        (i, θ_iter) in enumerate(θs), (j, φ_iter) in enumerate(φs)

        r0 .= [r_iter, θ_iter, φ_iter]
        r = r0
        b = magfield(r, nr, nθ, nφ)
        if b[1] > 0
            pol = 1.0
        else
            pol = -1.0
        end
        n = 0

        bmag0 = norm(b)
        b1rs[k,i,j,2] = bmag0
        bmag = bmag0
        if k == 1 && i == 1 && j == 1
            @debug b
        end
        r_min[1] = 100000
        while r[1] > 1.000009
            r[2] = acos(cos(r[2]))
            r[3] = atan2(sin(r[3]), cos(r[3]))
            if r[3] < 0
                r[3] += 2π
            end
            r1 = r
            b = magfield(r, nr, nθ, nφ)
            bmag = norm2(b)
            #v = -pol*b/bmag
            #v[2] = v[2]/r1[1]
            #sinθ = sin(r1[2])
            #if sinθ == 0
            #   sinθ = 1e-6
            #end
            #v[3] = v[3] / (r1[1]*sinθ)
            dt = dr * 0.1
            #r = rk4(r1, v, 3, t, dt, vfunc, pol)
            r = rk4(vfunc, t, r1, dt, pol)
            n = n + 1

            if r[2] < 0
                r[2] = -r[2]
                r[3] = r[3] + π
                if r[3] > 2π
                    r[3] -= 2π
                end
            end

            if r[1] < r_min[1]
                r_min .= r
            end
            if (r[1] - r_min[1]) > 2.5 || n > 10000 # wrong direction?
                pol = -pol #try back with opposite polarity
                n1 = 0
                r = r0
                r_min1[1] = 100000.0
                while r[1] > 1.000009
                    r[2] = acos(cos(r[2]))
                    r[3] = atan2(sin(r[3]), cos(r[3]))
                    if r[3] < 0
                        r[3] += 2π
                    end
                    r1 = r
                    b = magfield(r, nr, nθ, nφ)
                    bmag = norm2(b)
                    # v = -pol * b / bmag
                    # v[2] = v[2] / r1[1]
                    # sinθ = sin(r1[2])
                    # if sinθ == 0
                    #     sinθ = 1e-6
                    # end
                    # v[3] /= r1[1] * sinθ
                    dt = dr * 0.1
                    #r = rk4(r1, v, 3, t, dt, vfunc, pol)
                    r = rk4(vfunc, t, r1, dt, pol)
                    n1 = n1 + 1
                    if r[2] < 0
                        r[2] = -r[2]
                        r[3] = r[3]+π
                        if r[3] > 2π
                            r[3] = r[3] - 2π
                        end
                    end
                    if r[1] < r_min1[1]
                        r_min1 .= r
                    end

                    # double open field or stuck
                    if (r[1] - r_min1[1]) > 2.5 || n1 > 10000
                        if r_min[1] < r_min1[1]
                            pol = -pol
                            r = r_min # use r_min to remap
                        else
                            pol = pol
                            r = r_min1
                        end
                        lr = floor((r[1] - 1) / 0.01)
                        lθ = floor(r[2] / π* 180)
                        lφ = floor(r[3] / π* 180)
                        if ((lr*(nθ+1)+lθ) * (nφ+1) + lφ) < ((k*(nθ+1) + i) * (nφ+1) + j)
                            map[k,i,j] = (lr*(nθ+1)+lθ) * (nφ+1) + lφ
                        else
                            @info("New x point at", lr, lθ, lφ)
                            map[k,i,j] = ((lr-1)*(nθ+1)+lθ) * (nφ+1) + lφ
                        end
                        bmag = 0.0
                        break
                    end
                end
                break
            end
        end

        #nan1 = -3.0
        #nan = sqrt(nan1)
        b1rs[k,i,j,1] = pol * bmag
        map[k,i,j] = pol * map(k,i,j)
    end

end

@doc raw"""
    vfunc(t, r, pol)

Velocity function. Find dr/dt given r and t (t not really needed here).

```math
v = -\frac{\sgn(B_r)}{|B|} \left(B_r \hat{r} + \frac{B_θ}{r} \hat{θ} + \frac{B_φ}{r \sin θ} \hat{φ}\right)
```
"""
function vfunc(t, r, pol)
    b = magfield(r, nr, nθ, nφ)
    bmag = norm2(b)
    if bmag == 0
        return zeros(3)
    end
    v = -pol * b / bmag
    v[2] /= r[1]
    sinθ = sin(r[2])
    if sinθ == 0
        sinθ = 1.0e-6
    end
    v[3] /= r[1] * sinθ
    return v
end

"""
    magfield(rvec, nr, nθ, nφ)

Get ``\\mathbf{B} = \\mathbf{B}(\\mathbf{r})``
"""
function magfield(rvec, nr, nθ, nφ)
    r = rvec[1]
    θ = rem(rvec[2], π, RoundToZero)
    φ = mod2pi(rvec[3]) # make sure it's between 0 and 2π

    #  find the grid cell
    ir = floor((r-1) / (RSS-1) * n_r)
    if ir >= nr
        ir = nr - 1
    end

    iθ = floor(θ / π * nθ)
    if iθ >= nθ
        iθ = nθ - 1
    end

    iφ = floor(φ / 2π * nφ)
    if iφ >= nφ
        iφ = nφ - 1
    end

    if ir < 0
        # going into the sun, stop at the surface
        return magfieldgrid[0, iθ, iφ, :]
    end

    # relative displacement from lower grids
    px = [(r-1.0) / (RSS-1.0) * nr - ir,
          θ / pi * N_θ - iθ,
          φ / twopi * N_φ - iφ]

    fc = zeros(2,2,2)
    b = zeros(3)
    for m in 1:3
        fc[1,1,1] = magfieldgrid[ir,   iθ,   iφ,   m]
        fc[2,1,1] = magfieldgrid[ir+1, iθ,   iφ,   m]
        fc[1,2,1] = magfieldgrid[ir,   iθ+1, iφ,   m]
        fc[2,2,1] = magfieldgrid[ir+1, iθ+1, iφ,   m]
        fc[1,1,2] = magfieldgrid[ir,   iθ,   iφ+1, m]
        fc[2,1,2] = magfieldgrid[ir+1, iθ,   iφ+1, m]
        fc[1,2,2] = magfieldgrid[ir,   iθ+1, iφ+1, m]
        fc[2,2,2] = magfieldgrid[ir+1, iθ+1, iφ+1, m]

        # interpolate the magnetic field
        b[m] = trilinear(fc, px)
    end

    return b
end

"""
    trilinear(fc, x)

Trilinear interpolation

### Arguments
- `fc`: the value of f at the corner of cubic box of side 1 (3-d array)
- `r`:  location inside the cube (0≤x≤1) or outside x<0 x>1
"""
function trilinear(fc, r)
    x, y, z = r

    f = (fc[1,1,1] * (1-x) * (1-y) * (1-z)
       + fc[2,1,1] *    x  * (1-y) * (1-z)
       + fc[1,2,1] * (1-x) *    y  * (1-z)
       + fc[1,1,2] * (1-x) * (1-y) *    z
       + fc[2,1,2] *    x  * (1-y) *    z
       + fc[1,2,2] * (1-x) *    y  *    z
       + fc[2,2,1] *    x  *    y  * (1-z)
       + fc[2,2,2] *    x  *    y  *    z)

    return f
end

end
