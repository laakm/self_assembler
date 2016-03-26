using Winston

xgrid = linspace(-3.0, 3.0, 500)
ygrid = xgrid
xy = zeros(length(xgrid), length(ygrid))

function circle(x0, y0, R)
    xyt = zeros(length(xgrid), length(ygrid))
    for j = 1:length(ygrid)
        y = ygrid[j]
        if (y-y0)^2 < R^2
            for i = 1:length(xgrid)
                x = xgrid[i]
                if (x-x0)^2 + (y-y0)^2 < R^2
                    xyt[j,i] += 1
                end
            end
        end
    end

    return xyt
end

function annulus(x0, y0, R1, R2)
    xyt = zeros(length(xgrid), length(ygrid))
    
    xyt += circle(x0, y0, R1)
    xyt -= circle(x0, y0, R2)

    return xyt
end

function asym_circle(x0, y0, ang, R1, R2, a)
    xyt = zeros(length(xgrid), length(ygrid))

    xyt += circle(x0, y0, R1)
    xyt -= circle(x0+a*sin(ang), y0+a*cos(ang), R2)

    return xyt
end


dxdy = (xgrid[2]-xgrid[1])*(ygrid[2]-ygrid[1])
function compute_area(xy) 
    if maximum(xy) == 1.0
        return 0.0
    end

    return countnz(clamp(xy .- maximum(xy)*0.999, 0.0, 1.0))*dxdy
end
   


#initial background shape
#annulus case
#xyt = annulus(0.0, 0.0, 1.0, 0.4) #initial non-moving annulus
#mshape(x,y,ang) = annulus(x, y, 1.0, 0.4) #second moving annulus

#asymmetric circles
xyt = asym_circle(0.0, 0.0, 0.0, 1.0, 0.3, 0.3)
mshape(x,y,ang) = asym_circle(x, y, ang, 1.0, 0.3, 0.3)

#initial location for second obj
xx = 1.99
yy = 0.0
chi = 0.0

drs = 0.1
dangs = 0.1

Nsteps = 20
areas = zeros(Nsteps)
for k = 1:Nsteps
     
     #search maximum area increment direction
     mdarea = 0.0
     bang = 0.0
     for ang = linspace(0.0, 2pi, 30)
        #println("ang: $ang")
        dxs = drs*sin(ang)
        dys = drs*cos(ang)

        xy = xyt
        xy += mshape(xx+dxs, yy+dys, chi)
        darea = compute_area(xy)
        if darea > mdarea
            #println("  found better dir at: $ang")
            bang = ang
            mdarea = darea
        end
     end 

    #move into best dir
    xx += drs*sin(bang)
    yy += drs*cos(bang)

    #search best angle
    dchi = 0.0
    dchi_lim = 0.1
    mdarea = 0.0
    #for ang = linspace(chi-dchi, chi+dchi, 10)
    while abs(dchi) < dchi_lim
        xy = xyt
        xy += mshape(xx, yy, ang)

        darea = compute_area(xy)
        if darea > mdarea
            #println("  found better dir at: $ang")
            bang = ang
            mdarea = darea
        end
    end
  
    #optimum angle (instant turn)
    chi = bang


    xy = xyt
    xy += mshape(xx, yy, chi)

    areas[k] = compute_area(xy)
    p2 = imagesc(xy, (minimum(xy), maximum(xy)))
    display(p2)
    sleep(0.01)

end
