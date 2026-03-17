module Leappart

using LinearAlgebra

    #packages to use 

    

    #------------------------------------------------------------------------------------------------

    # compute the K partition coefficients

    function K(var, cal)

        TmH = cal.T0 .- cal.dTH2O .* (var.H2Om .^ cal.pH2O)
        cal.Tm  .= TmH  .* (1 .+ var.P ./ cal.A) .^ (1 ./ cal.B)

        L = (var.T .+ 273.15) .* cal.Dsx

        cal.Kx[:, 1:end-1] .= exp.( (L ./ cal.r) .* (1 ./ (var.T .+ 273.15) .- 1 ./ (cal.Tm .+ 273.15)) )'

        cal.Kf[:, cal.ncmp] .= 1 / (cal.H2Osat+1e-16)
        return cal
    end


    #------------------------------------------------------------------------------------------------
    # compute Tsolidus

    function Tsolidus(var,cal)

        var.H2Om = cal.H2Osat * (var.H2O>0)
        cal = K(var,cal)
        i = 1:(cal.ncmp - 1)

        if isnan(cal.Tsol) # just to check our guess is there and if its not there then use one inbetween the min and max
            weighted = sum(var.c[i] .* cal.Tm)
            #Tsol = clamp!(weighted, minimum(cal.Tm),maximum(cal.Tm))
            Tsol = max(minimum(cal.Tm),min(maximum(cal.Tm),sum(var.c[i] .* cal.Tm)))
            #Tsol = clamp(weighted, Tmin = minimum(cal.Tm), Tmax = maximum(cal.Tm))
        else
            Tsol = cal.Tsol
        end
        

        var.T = Tsol
        cal = K(var,cal)

        f = var.H2O

        rnorm = 1
        n = 0
        its_tol = 200

        while rnorm > cal.tol

            var.T = Tsol
            cal = K(var,cal)
            TK = Tsol + 273.15
            TmK = cal.Tm .+ 273.15

            E = -(var.c[i] .* exp.((cal.Dsx .* (TK .- TmK)) ./ (TmK .* cal.r))) ./ (f .- 1)
            r = sum(E) .+ var.H2Om .- 1

            drdT = sum( -(cal.Dsx .* var.c[i] .* exp.((cal.Dsx .* (TK .- TmK)) ./ (TmK .* cal.r))) ./ (TmK .* cal.r .* (f.-1)) )

            Tsol = Tsol - r/drdT/1.1

            rnorm = abs(r)/Tsol

            n += 1
            if n == its_tol
                error("!!! Newton solver for solidus T not converged after $its_tol iterations; rnorm = $rnorm !!!")
            end
        end

        return Tsol
    end

    #------------------------------------------------------------------------------------------------

    # use same for liquidus
    function Tliquidus(var,cal)

        var.H2Om = min(var.H2O,cal.H2Osat)
        cal = K(var,cal)
        i = 1:(cal.ncmp - 1)

        if isnan(cal.Tliq) # just to check our guess is there and if its not there then use one inbetween the min and max
            weighted = sum(var.c[i] .* cal.Tm)
            Tliq = max(minimum(cal.Tm),min(maximum(cal.Tm),sum(var.c[i] .* cal.Tm)))
        else
            Tliq = cal.Tliq
        end

        var.T = Tliq
        cal = K(var,cal)

        f = (var.H2O - var.H2Om) ./ (1 - var.H2Om)
    
        rnorm = 1
        n = 0
        its_tol = 100

            while rnorm > cal.tol

                var.T = Tliq
                cal = K(var,cal)
                TK = Tliq + 273.15
                TmK = cal.Tm .+ 273.15

                A = -(var.c[i] .* exp.(-(cal.Dsx .* (TK .- TmK)) ./ (TmK .* cal.r))) ./ (f-1)
                r = sum(A) .- 1

                drdT = sum( (cal.Dsx .* var.c[i] .* exp.(-(cal.Dsx .* (TK .- TmK)) ./ (TmK .* cal.r))) ./ (TmK .* cal.r .* (f.-1)) )

                Tliq = Tliq .- r/drdT/1.1

                rnorm = abs(r) / Tliq

                n += 1
                if n == its_tol
                error("!!! Newton solver for liquidus T not converged after $its_tol iterations; rnorm = $rnorm !!!")
            end

        end

        return Tliq

    end

    #------------------------------------------------------------------------------------------------

    # compute overlall equilibrium(var,cal) 

    function equilibrium(var,cal,Tsol,Tliq)

        var.H2Om = max(0, min(cal.H2Osat, var.H2O ./ (var.m + 1e-16)))
        i = 1:(cal.ncmp - 1)
        if !isnan(Tsol)
            var.T = max(Tsol,min(Tliq, var.T))
        end

        var.H2Om = max(0, min(cal.H2Osat, var.H2O ./ (var.m + 1e-16))) 
        var.cf = [0, 0, 1]
        var.f = max(0,min(var.H2O, var.H2O - cal.H2Osat))
        var.m = max(0,min(1-var.f,var.m))
        var.x = max(0,min(1,1-var.m-var.f))

        rnorm = 1
        n = 0
        its_tol = 100
        eps = 1e-9

        while rnorm > cal.tol
        
            cal = K(var,cal)
            res = sum((var.c .* (1 .- cal.Kx)) ./ (cal.Kx .+ (1 .- cal.Kx) .* var.m .+ (cal.Kf .- cal.Kx) .* var.f .+ 1e-32))

            varp = deepcopy(var) #melt fraction is nudged slightly (plus)
            varp.m = var.m + eps
            varp.H2Om = max(0,min(cal.H2Osat, varp.H2O ./ (varp.m .+ 1e-32)))
            varp.f = max(0,min(varp.H2O, varp.H2O - varp.m .* cal.H2Osat))
            calp = K(varp,deepcopy(cal))
            resp = sum((varp.c .* (1 .- calp.Kx)) ./ (calp.Kx .+ (1 .- calp.Kx) .* varp.m .+ (calp.Kf .- calp.Kx) .* varp.f .+ 1e-32))

            varm = deepcopy(var) #melt fraction is nudged slightly (minus)
            varm.m = var.m - eps
            varm.H2Om = max(0,min(cal.H2Osat, varm.H2O ./ (varm.m .+ 1e-32)))
            varm.f = max(0,min(varm.H2O, varm.H2O .- varm.m .* cal.H2Osat))
            calm = K(varm,deepcopy(cal))
            resm = sum((varm.c .* (1 .- calm.Kx)) ./ (calm.Kx .+ (1 .- calm.Kx) .* varm.m .+ (calm.Kf .- calm.Kx) .* varm.f .+ 1e-32))

            drdm_n = (resp .- resm) ./ (2 * eps)

            upd_m = -cal.alpha .* res ./ drdm_n
            upd_m = max(-var.m/2, min((1 .- var.m)/2, upd_m))
            var.m = max(0, min(1, var.m .+ upd_m))

            var.H2Om = max(0, min(cal.H2Osat, var.H2O ./ (var.m + 1e-32)))
            var.f = max(0, min(var.H2O, var.H2O - var.m * cal.H2Osat))
            var.x = max(0, min(1, 1 - var.m - var.f))

            rnorm = norm(upd_m) / norm(ones(size(upd_m)))

            n += 1
            if n == its_tol
                error("!!! Newton solver for liquidus T not converged after $its_tol iterations; rnorm = $rnorm !!!")
            end
        end

        comp_denom = vec(var.m .+ var.x .* cal.Kx .+ var.f .* cal.Kf .+ 1e-16)
        var.cm = max.(0, min.(1, var.c ./ comp_denom))
        var.cx = max.(0, min.(1, (var.c .* vec(cal.Kx)) ./ comp_denom))

        return var,cal
    end

    #------------------------------------------------------------------------------------------------
    # leppart function

    function leappart(var,cal)

        cal = K(var,cal)

        Tsol = Tsolidus(var,cal)
        Tliq = Tliquidus(var,cal)

        var,cal = equilibrium(var,cal,Tsol,Tliq)

        return var,cal,Tsol,Tliq
    end


    #------------------------------------------------------------------------------------------------
    #ATTEMPTING TO CALL FUNCTIONS -------------------------------------------------------------------

    #create the var and cal structure (ignoring type for now - will just do all in one)

    mutable struct Var 
        H2O::Float64         # water content
        c::Vector{Float64}   # composition
        T::Float64           # temperature
        P::Float64           # pressure
        H2Om::Float64        # water dissolved in melt

        x::Float64           # crystalline fraction
        m::Float64           # melt fraction
        f::Float64           # fluid fraction

        cx::Vector{Float64}  # melt composition
        cm::Vector{Float64}  # crystal composition
        cf::Vector{Float64}  # fluid composition
    end

    mutable struct Cal 
        ncmp::Int # number of components 

        T0::Vector{Float64} # zero-water and -pressure of melting temperatures 
        A::Vector{Float64} # pressure param controlling how strongly pressure increases Tm
        B::Vector{Float64} # pressure param controls curvature of pressure effect 
        r::Vector{Float64} # a scaling constant related to the gas constant 
        Dsx::Float64 # entrpoy for each component 
        pH2O::Float64 # exponent for nonlinear dependence of dissolved H2O effecting melting point
        dTH2O::Vector{Float64} # controls strength of the dissolved H2O on temperature of solidus
        H2Osat::Float64 # how much melt the water can hold (max dissolved water) 

        tol::Float64 # tolerance for loop 
        alpha::Float64 # damping factor for Newton solver 

        Kx::Matrix{Float64} # partition coefficient of the composition phase 
        Kf::Matrix{Float64} # partition coefficient of the fluid phase 
        Tm::Vector{Float64} # melting temperatures of components 
        Tsol::Float64 # temperature of solidus 
        Tliq::Float64 # temperature of liquidus 
    end

    #------------------------------------------------------------------------------------------------
    #setting variables in var

    function set_var(Var)

        H2O = 0.05                # water content
        c = [0.3, 0.7, H2O]       # composition of magma components in fraction form with adjustable water content
        c[1:end-1] .= c[1:end-1]./sum(c[1:end-1]).*(1-H2O) # normalise the composition to account for water
        T = 1200.0                # temperature of system
        P = 1.0                   # pressure of system
        H2Om = 0.00               # initial dissolved water

        x = 0.1                   # initial crystal fraction
        m = 0.6                   # initial melt fraction
        f = 0.3                   # initial fluid fraction

        cx = copy(c)              # melt composition (initial guess)
        cm = copy(c)              # crystal composition (initial guess)
        cf = copy(c)              # fluid composition (initial guess)

        var = Var(H2O, c, T, P, H2Om, x, m, f, cx, cm, cf)
        return var 
    end

    #setting variables in cal

    function set_cal(Cal)
        
        ncmp = 3                  # set three components when including water

        T0   = [1400, 1000]   # dry zero-P melting temperatures of components (°C)
        A    = [5, 3]  
        B    = [5, 3]   
        r    = [10, 6]   
        Dsx = 300
        pH2O = 0.75
        dTH2O = 1400 .* 1200 ./ T0
        H2Osat = 0.05              

        tol = 1e-6                # Newton tolerance
        alpha = 0.5               # damping factor for Newton solver

        Kx = zeros(1,ncmp)   # initialise Kx matrix
        Kf = zeros(1,ncmp)  # initialise Kf matrix
        Tm = [NaN,NaN]  # intialise melting temperature
        Tsol = NaN                # initialise solidus
        Tliq = NaN                # initialise liquidus

        cal = Cal(ncmp, T0, A, B, r, Dsx, pH2O, dTH2O, H2Osat, tol, alpha, Kx, Kf, Tm, Tsol, Tliq)
        return cal 
    end

    #-------------------------------------------------------------------------------------------------- 
    # run functions and print 

    #var = set_var(Var); 
    #cal = set_cal(Cal); 
    #cal = K(var,cal); 
    #Tsol = Tsolidus(var,cal); 
    #Tliq = Tliquidus(var,cal); 
    #var,cal,Tsol,Tliq = leappart(var,cal)


    #println("""
    #        K partition coefficients: Kx = $(cal.Kx), Kf = $(cal.Kf[3])
    #        Tsolidus = $Tsol
    #        Tliquidus = $Tliq
    #        Tm = $(cal.Tm)
    #        Melt fraction m = $(var.m)
    #        Crystal fraction x = $(var.x)
    #        Fluid fraction f = $(var.f)
    #        """)


end # module Leappart
