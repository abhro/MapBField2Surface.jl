module MapBField2Surface

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

    for k in 0:nr
        r0[1] = 1.00001 + k*0.01
        for i in 0:nθ
            r0[2] = i * DEG_TO_RAD
            for j = 0:nφ
                r0[3] = j * DEG_TO_RAD
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
                rmin[1] = 100000
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
                    #   sinθ = 1.0d-6
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

                    if r[1] < rmin[1]
                        rmin = r
                    end
                    if (r[1] - rmin[1]) > 2.5 || n > 10000 # wrong direction?
                      pol = -pol #try back with opposite polarity
                      n1 = 0
                      r = r0
                      rmin1[1] = 100000.0
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
                        # sinθ = sin(r1(2))
                        # if sinθ == 0
                        #     sinθ = 1.0d-6
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
                        if r[1] < rmin1[1]
                            rmin1 = r
                        end

                        # double open field or stuck
                        if (r[1] - rmin1[1]) > 2.5 || n1 > 10000
                          if rmin[1] < rmin1[1]
                            pol = -pol
                            r = rmin # use rmin to remap
                          else
                            pol = pol
                            r = rmin1
                          end
                          lr = floor((r[1] - 1.0d0) / 0.01)
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
    end

end

# velocity function. find dr/dt given r and t (t not really needed)
# v = -sgn(B_r)/|B| (B_r r_hat + B_θ / r θ_hat + B_φ / (r sin θ) φ_hat)
function vfunc(t, r, pol)
    b = magfield(r, nr, nθ, nφ)
    bmag = norm2(b)
    if bmag == 0
        return zeros(3)
    end
    v = -pol * b / bmag
    v[2] /= r(1)
    sinθ = sin(r[2])
    if sinθ == 0
        sinθ = 1.0e-6
    end
    v[3] /= r[1] * sinθ
    return v
end

# get \vec{B} = \vec{B}(\vec{r})
function magfield(rvec, nr, nθ, nφ)
    r = rvec[1]
    θ = rvec[2]
    φ = rvec[3]
    θ = acos(cos(θ)) # XXX replace with rem2pi or mod2pi
    φ = atan2(sin(φ), cos(φ)) # XXX replace with rem2pi or mod2pi
    if φ < 0.0
        φ += 2π
    end

    #  find the grid cell
    irr = floor((r-1) / (RSS-1) * n_r)
    if irr >= nr
        irr = nr - 1
    end

    iθ = floor(θ / π * nθ)
    if iθ >= nθ
        iθ = nθ - 1
    end

    iφ = floor(φ / 2π * nφ)
    if iφ >= nφ
        iφ = nφ - 1
    end

    if irr < 0
      # going into the sun, stop at the surface
      return magfieldgrid[0, iθ, iφ, :]
    end

    #  relative displacement from lower grids
    px = [(r-1.0) / (RSS-1.0) * nr - irr,
          θ / pi * N_θ - iθ,
          φ / twopi * N_φ - iφ]

    fc = zeros(2,2,2)
    b = zeros(3)
    for m in 1:3
      fc[1,1,1] = magfieldgrid[irr,   iθ,   iφ,   m]
      fc[2,1,1] = magfieldgrid[irr+1, iθ,   iφ,   m]
      fc[1,2,1] = magfieldgrid[irr,   iθ+1, iφ,   m]
      fc[2,2,1] = magfieldgrid[irr+1, iθ+1, iφ,   m]
      fc[1,1,2] = magfieldgrid[irr,   iθ,   iφ+1, m]
      fc[2,1,2] = magfieldgrid[irr+1, iθ,   iφ+1, m]
      fc[1,2,2] = magfieldgrid[irr,   iθ+1, iφ+1, m]
      fc[2,2,2] = magfieldgrid[irr+1, iθ+1, iφ+1, m]

      # interpolate the magnetic field
      b[m] = trilinear(fc, px)
    end

    return b
end

# trilinear interpolation
function trilinear(
        fc,     # fc the value of f at the corner of cubic box of side 1
        x)      # location inside the cube (0<=x<=1) or outside x<0 x>1

    f = (fc[1,1,1] * (1-x[1]) * (1-x[2]) * (1-x[3])
       + fc[2,1,1] *    x[1]  * (1-x[2]) * (1-x[3])
       + fc[1,2,1] * (1-x[1]) *    x[2]  * (1-x[3])
       + fc[1,1,2] * (1-x[1]) * (1-x[2]) *    x[3]
       + fc[2,1,2] *    x[1]  * (1-x[2]) *    x[3]
       + fc[1,2,2] * (1-x[1]) *    x[2]  *    x[3]
       + fc[2,2,1] *    x[1]  *    x[2]  * (1-x[3])
       + fc[2,2,2] *    x[1]  *    x[2]  *    x[3])

    return f
end


end
