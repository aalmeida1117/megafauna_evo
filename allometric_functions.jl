module AllometricFunctions

    export find_velocity, bite_size_allo, number_of_chews, chew_rate_allo, chew_allo,
        gut_volume_g, foragingtime, bc, energy_cost, mean_energy_closed, mean_energy_open,
        reactionwidth, reactionheight, indperarea, expectedlifetime, predation_mortality_closed, predation_mortality_open,
        multiplier_closed, prob_pred_closed, prob_pred_open
 

    using Base.Threads
    using UnicodePlots
    using Distributions

    function find_velocity(mass)
        # mass in kg
        # from kramer 2010 "allometric scaling of resource acquisition"
        #Consumer Velocity (meters/second)
        velocity = (0.5 * mass^0.13);
        return velocity
    end

    function bite_size_allo(mass)
        # mass in kg
        # bite size in g
        bite_size = 0.002 * ((mass)^0.969); # g
        return bite_size

    end

    function number_of_chews(mass)
        # shipley 94
    #chew/g (processed to mean particle size (allo) and swallowed)
    # mass in kg
        chews_per_gram = 343.71 * (mass^(-0.83)); 
        return chews_per_gram
    end

    function chew_rate_allo(mass)
        # from "dental functional morphology predicts scaling"
        # mass in kg
        # duration in ms
        chewing_cycle_duration = (228.0* (mass)^0.246) / 1000;
        return 1 / (chewing_cycle_duration )  #[s/chew -> chews/s]
    
    end

    function chew_allo(mass)
        # allometric function for mouth->gut rate
    chew_rate = chew_rate_allo(mass); # chews/s
    chews_per_gram = number_of_chews(mass)   # chew/kg
    chew = chew_rate / chews_per_gram        # chew/s / chew/g -> g/s
    return chew
    end

    function gut_volume_g(mass)
        capacity = (0.030 * (mass)^0.881);
        return capacity * 1000 #[kg -> g]
        #return capacity
    end

    # function indperarea(mass)
    #     #Enter mass in kg
    #     #Convert to grams
    #     massg= mass*1000;
    #     popdensity= (0.0116)*massg^-0.776; #inds/area (from Damuth)
    #     return popdensity
    # end

    function foragingtime(mass)
        #Owen Smith 1988 book (data grabbed using PlotDigitizer in /data/)
        #mass in kg
        #foraging time in % of 24 hours

        forageperc = 21.0885*mass^0.124
        forageperc_upper95 = 27*mass^0.17

        #translate to hours
        foragehrs = (forageperc/100)*24
        foragehrs_upper95 = (forageperc_upper95/100)*24

        return foragehrs, foragehrs_upper95 # hours
    end

    function bc(yvalue, yc, ymax)
        return clamp(yvalue, yc, ymax)
    end

    function energy_cost(mass)
        ## energy in kcal
        energy = 0.15*(mass^0.68)
        return energy
    end

    function mean_energy_closed(mass)
        #enter mass in kg
    #massg=mass*1000 #convert to grams
        mean_edensity_closed=30*(mass^-0.24)#kcal
        return mean_edensity_closed
    end

    function mean_energy_open(mass)
        #enter mass in kg
    #massg=mass*1000 #convert to grams
        mean_edensity_open=20*(mass^-0.24)#kcal
        return mean_edensity_open
    end

    #reaction distance
    function reactionwidth(mass)
        #mass in kg
        #distance in meters
        #see 'fit' with Pawar data in analysis
        #Note it is not a fit, but anchoring the intercept 

        width = 2 * 5 * (mass^(1/3)) #meters

        return width

    end
    # reaction height
    function reactionheight(mass)
        #mass in KG
        #distance in meters ~ shoulder height
        #From Larramendi 2016
        ra = 0.1501*(mass^0.357)
        return ra
    end

    function indperarea(mass)
        # mass in kg
        # Pop density: Damuth 1981
        popdensity = (5.45e-5) * mass^-0.776  # inds/area

        # Foraging bout time (s)
        foragebout_hrs, _ = foragingtime(mass)
        foragebout_s = foragebout_hrs * 3600

        # Foraging velocity (m/s)
        foragevelocity = find_velocity(mass)

        # Area foraged as corridor (m²)
        corridorwidth = 2  # meters
        corridorareaforaged = foragevelocity * foragebout_s * corridorwidth

        # OR: use Owen-Smith homerange model (commented out)
        # HR = 13500 * mass^1.25

        # Total individuals in foraged area
        popinarea = 1.0 + popdensity * corridorareaforaged

        return popdensity, popinarea
    end
    function expectedlifetime(mass)
        #Calder 1984;
        #mass in kg
        #expected lifetime in years
        explifetime = 5.08*mass^0.35;
        explifetime_days = explifetime *365
        return explifetime_days
    end
    # --- OPEN habitat (ajuste baseado nos dados do Serengeti) ---
    function predation_mortality_open(mass)
        a, b, ymin, ymax = 14.9738, -6.02, 0.05, 0.95
        #a, b, ymin, ymax = 20.1, -8.0, 0.05, 0.9
        p = 1 / (1 + exp(-(a + b * log10(mass))))
        return ymin + (ymax - ymin) * p   # fração anual (0–1)
    end

    function prob_pred_open(mass)
        annual = predation_mortality_open(mass)
        base = max(1 - annual, 0.0)
        daily_survival = base^(1/365)
        return 1 - daily_survival
    end
    # --- multiplicador para closed habitat ---
    function multiplier_closed(mass; a=5.0, b=-1.3, low=0.1, high= 0.9)
        sig = 1 / (1 + exp(-(a + b * log10(mass))))
        return low + (high - low) * sig
    end

    # --- CLOSED habitat ---
    function predation_mortality_closed(mass; a=5.0, b=-1.3, low=0.1, high=0.9)
        base = predation_mortality_open(mass)
        mult = multiplier_closed(mass; a=a, b=b, low=low, high=high)
        return clamp(base * mult, 0.0, 0.99)
    end

    function prob_pred_closed(mass; a=5.0, b=-1.3, low=0.1, high=0.9)
        annual = predation_mortality_closed(mass; a=a, b=b, low=low, high=high)
        base = max(1 - annual, 0.0)
        daily_survival = base^(1/365)
        return 1 - daily_survival
    end

end 