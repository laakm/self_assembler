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
    xyt -= circle(x0+a*cos(ang), y0+a*sin(ang), R2)

    return xyt
end


function spiral(x0, y0, ang, R, b)
    xyt = zeros(length(xgrid), length(ygrid))

    for j = 1:length(ygrid)
        y = ygrid[j]
        if (y-y0)^2 < (R + (2pi*b))^2
            for i = 1:length(xgrid)
                x = xgrid[i]
                xyang = atan2(y-y0,x-x0)
                if (x-x0)^2 + (y-y0)^2 < (R + b*mod2pi(ang+xyang))^2
                    xyt[j,i] += 1
                end
            end
        end
    end

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
#xyt = asym_circle(0.0, 0.0, 0.0, 1.0, 0.4, 0.24)
#mshape(x,y,ang) = asym_circle(x, y, ang, 1.0, 0.4, 0.24)

xyt = spiral(0.0, 0.0, 0.0, 1.0, 0.1)
mshape(x,y,ang) = spiral(x,y,ang, 1.0, 0.1)

#initial location for second obj
xx = 1.9
yy = 0.0
chi = pi/2

xx = 0.0
yy = 0.0
chi = pi/2

drs = 0.1

mdarea = 0.0 #current biggest area
Nsteps = 50
areas = zeros(Nsteps)
darea_move =Float64[]
darea_rot = Float64[]
for k = 1:Nsteps
     
     #search maximum area increment direction
     make_step = false
     bang = 0.0
     mdarea_step = 0.0
     for ang = linspace(0.0, 2pi, 30)
        #println("ang: $ang")
        dxs = drs*cos(ang)
        dys = drs*sin(ang)

        xy = xyt
        xy += mshape(xx+dxs, yy+dys, chi)
        darea = compute_area(xy)
        if darea > mdarea
            #println("  found better dir at: $ang")
            bang = ang
            mdarea_step = darea - mdarea
            mdarea = darea
            make_step = true
        end
     end 

    #move into best dir
    if make_step    
#    if false
        println("moving to: ",round(drs*cos(bang),3)," ", round(drs*sin(bang),3))
        push!(darea_move, mdarea_step)
        xx += drs*cos(bang)
        yy += drs*sin(bang)
    else
        println("no move")
    end

    #search best angle
    make_rot = false
    dchi = 0.0
    dchi_lim = 0.1 #max rotation that the system makes
    dchis = dchi_lim/10 #max rotation step that the system can make
   
    #rotate to left 
    xy = xyt
    xy += mshape(xx, yy, chi-dchis)
    darea_left = compute_area(xy)
    println("rot left: ",darea_left, " ",mdarea," diff:",mdarea-darea_left)

    #rotate to right
    xy = xyt
    xy += mshape(xx, yy, chi+dchis)
    darea_right = compute_area(xy)
    println("rot right: ",darea_right, " ",mdarea," diff:",mdarea-darea_right)

    if mdarea >= max(darea_right, darea_left)
        println("no rot")
    else 
        make_rot = true

        rot_sign = 0.0
        if darea_right > darea_left
            rot_sign = 1.0
            push!(darea_rot, darea_right - mdarea)
            mdarea = darea_right
        else
            rot_sign = -1.0
            push!(darea_rot, darea_left - mdarea)
            mdarea = darea_left
        end
        chi = chi + rot_sign*dchis

        #rotate to correct dir until step size reached
        dchi = 0.0
        while abs(dchi) < dchi_lim
            dchi = dchi + rot_sign*dchis

            xy = xyt
            xy += mshape(xx, yy, chi+dchi)
            darea = compute_area(xy)
            if darea < mdarea
                println("stopping rotation")
                break
            else
                println("rotating $(round(dchi,3))")
                push!(darea_rot, darea - mdarea)
                mdarea = darea
            end
        end
        chi = chi + rot_sign*dchi
    end

    xy = xyt
    obj_area = countnz(xy)*dxdy
    xy += mshape(xx, yy, chi)

    areas[k] = compute_area(xy)
    p2 = imagesc(xy, (minimum(xy), maximum(xy)))

    #plot derivatives
    p3 = plot(darea_move, "b-")

    p4 = plot(darea_rot, "r-")

    t = Table(3,1)
    t[1,1] = p2
    t[2,1] = p3
    t[3,1] = p4
    display(t)

    println("area frac: ",round(areas[k]/obj_area,3))
    if !(make_step || make_rot)
        println("equilibrium found")
        break
    end
end
