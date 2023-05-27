using JuMP, Ipopt
#using NLopt
using DataFrames
#m = Model(NLopt.Optimizer)
#set_optimizer_attribute(m, "algorithm", :LD_MMA)
const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter"   => 9999)






function solveRSMaddCentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, ζ >= 0)

    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pc, (μ + a + b*c + b*α*(hm-hc))/(2b));
    @NLexpression(m, goal, (pc - v - α*hm)*F⁻ - (c-v) );

    @NLconstraint(m, goal <= ζ);
    @NLconstraint(m, goal >= -ζ);

    @objective(m, Min, ζ)

    JuMP.optimize!(m)
    @NLexpression(m, obj, (pc-v-α*hm)*(μ+a-b*(pc+α*hc))-(c-v)*(z+a-b*(pc+α*hc)))

    (;z=value(z), pc=value(pc), ζ = value(ζ),obj=value(obj), status=termination_status(m))
end
#solveRSMaddCentral()


# %%
function solveRSMaddCentral2(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, p >= 0)
    @variable(m, ζ >= 0)
    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));

    @NLobjective(m, Max, (p-v-α*hm)*(μ+a-b*(p+α*hc))-(c-v)*(z+a-b*(p+α*hc)))

    JuMP.optimize!(m)


    (;z=value(z), p=value(p), obj = objective_value(m), status=termination_status(m))
end
#solveRSMaddCentral2()


function solveRSMaddDecentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)

    @NLexpression(m, F, (z-A)/(B-A));
    @NLexpression(m, F⁻, 1-F);
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pd, (μ + a - b*α*hc + b*v*F)/(b*(1 + F)));
    @NLexpression(m, w, r/(b*(1 + F)) * ((μ + a - b*α*hc)*F⁻ + 2b*v*F)  ) ;


    @NLobjective(m, Max, 
        ((1 - r)*(pd - v) - α*hm)*(μ + a - b*(pd + α*hc )) 
          + (w - c + (1 - r)*v)*(z + a - b*(pd + α*hc))
    )

    JuMP.optimize!(m)

    (;z=value(z), pd=value(pd),w=value(w), obj = objective_value(m), status=termination_status(m))
end

function rsMaddDecentral(z;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    #@assert 2A + a - b*c - b*α*hc > 0
    #@assert A + a - b*c - b*α*(hc + 2hm) > 0
    
    F(z) = (z-A)/(B-A)
    F⁻(z) = 1-F(z);
    μ(z) = (B-z)*(B-z)/(2*(A-B))
    pd(z) = (μ(z) + a - b*α*hc + b*v*F(z))/(b*(1 + F(z)))
    w(z) = r/(b*(1 + F(z))) * ((μ(z) + a - b*α*hc)*F⁻(z) + 2b*v*F(z))  

    ((1 - r)*(pd(z) - v) - α*hm)*(μ(z) + a - b*(pd(z) + α*hc )) + (w(z) - c + (1 - r)*v)*(z + a - b*(pd(z) + α*hc))

end

"""
using Alpine, Gurobi

const gurobi = optimizer_with_attributes(Gurobi.Optimizer, 
                                         MOI.Silent() => true,
                                         "Presolve"   => 1) 



const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                         "nlp_solver" => ipopt,
                                         "mip_solver" => gurobi)
solveRSMaddCentral2(alpine)                                         
"""

function solveCRSMaddCentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, ζ >= 0)

    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pc, (μ + a + b*c + b*α*(hm-hc))/(2b));
    @NLexpression(m, goal, (pc - v - α*hm)*F⁻ - (c-v) );

    @NLconstraint(m, goal <= ζ);
    @NLconstraint(m, goal >= -ζ);

    @objective(m, Min, ζ)

    JuMP.optimize!(m)
    @NLexpression(m, obj, (pc-v-α*hm)*(μ+a-b*(pc+α*hc))-(c-v)*(z+a-b*(pc+α*hc)))

    (;z=value(z), pc=value(pc), ζ = value(ζ),obj=value(obj), status=termination_status(m))
end

function solveCRSMaddCentral2(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, p >= 0)
    @variable(m, ζ >= 0)
    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));

    @NLobjective(m, Max, (p-v-α*hm)*(μ+a-b*(p+α*hc))-(c-v)*(z+a-b*(p+α*hc)))
                         
    JuMP.optimize!(m)

    (;z=value(z), p=value(p), obj = objective_value(m), status=termination_status(m))
end

function crsMaddDecentral(z;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    F(z) = (z-A)/(B-A)
    F⁻(z) = 1-F(z);
    μ(z) = (B-z)*(B-z)/(2*(A-B))
    pd(z) = (μ(z) + a - b*α*hc + b*(α * hm + v)*F(z))/(b*(1 + F(z)))
    w(z) = r/(b*(1 + F(z))) * ((μ(z) + a - b*α*hc)*F⁻(z) - b*α*hm*F⁻(z) + 2b*v*F(z))  


    (1 - r)*(pd(z) - v  - α*hm)*(μ(z) + a - b*(pd(z) + α*hc )) + (w(z) - c + (1 - r)*v)*(z + a - b*(pd(z) + α*hc))

end

# %%
function solveCRSMaddDecentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hm=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)

    @NLexpression(m, F, (z-A)/(B-A));
    @NLexpression(m, F⁻, 1-F);
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pd, (μ + a - b*α*hc + b*(α*hm + v)*F)/(b*(1 + F)));
    @NLexpression(m, w, r/(b*(1 + F)) * ((μ + a - b*α*hc)*F⁻  - b*α*hm*F⁻ + 2b*v*F)  ) ;


    @NLobjective(m, Max, 
        (1 - r)*(pd - v - α*hm)*(μ + a - b*(pd + α*hc )) 
          + (w - c + (1 - r)*v)*(z + a - b*(pd + α*hc))
    )

    JuMP.optimize!(m)

    (;z=value(z), pd=value(pd),w=value(w), obj = objective_value(m), status=termination_status(m))
end

#################3
#
#  TEN JEST DULPKATAME solveRSMaddCentral (tylko inaczej sie nazywa parametr h)
#
#

function solveRSRaddCentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hr=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    #@assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, ζ >= 0)

    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pc, (μ + a + b*c + b*α*(hr-hc))/(2b));
    @NLexpression(m, goal, (pc - v - α*hr)*F⁻ - (c-v) );

    @NLconstraint(m, goal <= ζ);
    @NLconstraint(m, goal >= -ζ);

    @objective(m, Min, ζ)

    JuMP.optimize!(m)
    @NLexpression(m, obj, (pc-v-α*hr)*(μ+a-b*(pc+α*hc))-(c-v)*(z+a-b*(pc+α*hc)))

    (;z=value(z), pc=value(pc), ζ = value(ζ),obj=value(obj), status=termination_status(m))
end
#solveRSRaddCentral()


# %%
function solveRSRaddCentral2(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hr=0.1, 
    hc=0.1 ,
    v= c - 0.1)

    @assert 2A + a - b*c - b*α*hc > 0
    #@assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)
    @variable(m, p >= 0)
    @variable(m, ζ >= 0)
    @NLexpression(m, F⁻, 1-(z-A)/(B-A));
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));

    @NLobjective(m, Max, (p-v-α*hr)*(μ+a-b*(p+α*hc))-(c-v)*(z+a-b*(p+α*hc)))

    JuMP.optimize!(m)

    (;z=value(z), p=value(p), obj = objective_value(m), status=termination_status(m))
end
#solveRSRaddCentral2()



# %%
function solveRSRaddDecentral(opt=Ipopt.Optimizer;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hr=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    @assert 2A + a - b*c - b*α*hc > 0
    #@assert A + a - b*c - b*α*(hc + 2hm) > 0

    m = Model(opt)
    @variable(m, A <= z <= B)

    @NLexpression(m, F, (z-A)/(B-A));
    @NLexpression(m, F⁻, 1-F);
    @NLexpression(m, μ, (B-z)*(B-z)/(2*(A-B)));
    @NLexpression(m, pd, (r*(μ + a - b*α*hc) + b*(α*hr+r*v)*F)/(b*r*(1 + F)));
    @NLexpression(m, w, 1/(b*(1 + F)) * (r*(μ + a - b*α*hc)*F⁻ - b*α*hr*F⁻ + 2b*r*v*F)  ) ;


    @NLobjective(m, Max, 
        (1 - r)*(pd - v)*(μ + a - b*(pd + α*hc )) 
          + (w - c + (1 - r)*v)*(z + a - b*(pd + α*hc))
    )

    JuMP.optimize!(m)

    (;z=value(z), pd=value(pd),w=value(w), obj = objective_value(m), status=termination_status(m))
end
#solveRSRaddDecentral()

# %%

function rsRaddRSDecentral(z;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hr=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    @assert 2A + a - b*c - b*α*hc > 0
    @assert A + a - b*c - b*α*(hc + 2hr) > 0

    F(z) = (z-A)/(B-A)
    F⁻(z) = 1-F(z);
    μ(z) = (B-z)*(B-z)/(2*(A-B))
    pd(z) = (r*(μ(z) + a - b*α*hc) + b*(α*hr+r*v)F(z))/(b*r*(1 + F(z)))
    w(z) =  1/(b*(1 + F(z))) * (r*(μ(z) + a - b*α*hc)*F⁻(z) - b*α*hr*F⁻(z) + 2b*r*v*F(z))  

    (1 - r)*(pd(z) - v)*(μ(z) + a - b*(pd(z) + α*hc )) + (w(z) - c + (1 - r)*v)*(z + a - b*(pd(z) + α*hc))

end

#x = -4:0.05:4
#y = rsRaddRSDecentral.(x, A=-4, B=4)
#using Plots
#plot(x, y)

# %%

#
# RADD CRS
#

# %%

function crsRaddDecentral(z;
    A=-1,   # A < 0
    B=-A,
    a=10.0,
    b=0.6,
    c=0.5,
    α=0.05,
    hr=0.1, 
    hc=0.1 ,
    v= c - 0.1,
    r=0.4)

    #@assert 2A + a - b*c - b*α*hc > 0
    #@assert A + a - b*c - b*α*(hc + 2hm) > 0

    F(z) = (z-A)/(B-A)
    F⁻(z) = 1-F(z);
    μ(z) = (B-z)*(B-z)/(2*(A-B))
    pd(z) = (r*(μ(z) + a - b*α*hc) + b*(α*hr+r*v)F(z))/(b*r*(1 + F(z)))
    w(z) =  1/(b*(1 + F(z))) * (r*(μ(z) + a - b*α*hc)*F⁻(z) - b*α*hr*F⁻(z) + 2b*r*v*F(z))  

    (1 - r)*(pd(z) - v - α*hr)*(μ(z) + a - b*(pd(z) + α*hc )) + (w(z) - c + (1 - r)*v)*(z + a - b*(pd(z) + α*hc))

end



"""
wholesalePriceContractM(opt=Ipopt.Optimizer;
    A=-1, B=-A, a=10.0, b=0.6, c=0.5, α=0.05, hm=0.1, hc=0.1 , v= c - 0.1)  = 
    solveRSMaddDecentral(opt;A,B,a,b,c, α, hm, hc, v, r=1)

wholesalePriceContractR(opt=Ipopt.Optimizer;
    A=-1, B=-A, a=10.0, b=0.6, c=0.5, α=0.05, hr=0.1, hc=0.1 , v= c - 0.1)  = 
    solveRSRaddDecentral(opt;A,B,a,b,c, α, hr, hc, v, r=1)
"""
# %%

# pARAMETER SWEEP: 2X PO H PLUS ALPHA


# mmm => solveRSMaddCentral
# wwwM => wholesalePriceContractM
# wwwR => wholesalePriceContractR
# rrrM => solveRSMaddDecentral
# cccM => solveCRSMaddDecentral
# rrrR => solveRSRaddDecentral
# cccR=cccM cost-revenue sharing R= cost-revenue sharing M
# Madd central = Radd central



# model central manufacturer hand
wpcm = []
for hm in 1:0.5:3
    for hc in 2:0.5:4
        for α in 0.05:0.05:0.15
            for c in 0.5:0.5:2
                for v in 0.2:0.1:0.4
                    mmm = solveRSMaddCentral(;hm, hc, α, c, v)
                    res = (;hm, hc, α, c, v, mmm...)
                    push!(wpcm, res)
              end
           end
        end
    end
end
dfwpcm = DataFrame(wpcm)

using CSV
CSV.write("dfwpcm.csv", dfwpcm)

# model central retailer hand

wpcr = []
for hm in 3:0.5:5
    for hc in 1:0.5:2
        for α in 0.05:0.05:0.15
            for c in 0.5:0.5:2
                for v in 0.2:0.1:0.4
                    mmm = solveRSMaddCentral(;hm, hc, α, c, v)
                    res = (;hm, hc, α, c, v, mmm...)
                    push!(wpcr, res)
              end
           end
        end
    end
end
dfwpcr = DataFrame(wpcr)

using CSV
CSV.write("dfwpcr.csv", dfwpcr)

# wholesaleprice manufacturer hand

cm = []
for hm in 1:0.5:3
    for hc in 2:0.5:4
        for α in 0.05:0.05:0.15
            for c in 0.5:0.5:2
                for v in 0.2:0.1:0.4
                    wwwM = wholesalePriceContractM(;hm, hc, α, c, v)
                    res = (;hm, hc, α, c, v, wwwM...)
                    push!(cm, res)
               end
           end
        end
    end
end
dfcm = DataFrame(cm)
using CSV
CSV.write("dfcm.csv", dfcm)

using Plots
d = df[ (df.hc .== 2.0) .&& (df.α .== 0.15), : ]

# wholesaleprice retailer hand

cr = []
for hr in 3:0.5:5
    for hc in 1:0.5:2
        for α in 0.05:0.05:0.15
            for c in 0.5:0.5:2
                for v in 0.2:0.1:0.4
                    wwwR = wholesalePriceContractR(;hr, hc, α, c, v)
                    res = (;hr, hc, α, c, v, wwwR...)
                    push!(cr, res)
               end
           end 
        end
    end
end
dfcr = DataFrame(cr)

using CSV
CSV.write("dfcr.csv", dfcr)

# RS reveue-sharing manufacturer hand
rsm = []
for hm in 1:0.5:2
    for hc in 3:0.5:4
        for α in 0.1:0.05:0.15
            for c in 1:0.5:2
                for v in 0.2:0.1:0.4
                    for r in 0.5:0.5:1
                        rrrM = solveRSMaddDecentral(;hm, hc, α, c, v, r)
                        res = (;hm, hc, α, c, v, r, rrrM...)
                        push!(rsm, res)
                     end 
                end
            end
        end
    end
end
dfrsm = DataFrame(rsm)

using CSV
CSV.write("dfrsm.csv", dfrsm)

# CRS cost-reveue-sharing manufacturer hand

crsm = []
for hm in 1:0.5:2
    for hc in 3:0.5:4
        for α in 0.1:0.05:0.15
            for c in 1:0.5:2
                for v in 0.2:0.1:0.4
                    for r in 0.5:0.5:1
                    cccM = solveCRSMaddDecentral(;hm, hc, α, c, v, r)
                    res = (;hm, hc, α, c, v, r, cccM...)
                    push!(crsm, res)
                  end 
              end
           end
        end
    end
end
dfcrsm = DataFrame(crsm)

using CSV
CSV.write("dfcrsm.csv", dfcrsm)


using Plots
d = df[ (df.hc .== 2.0) .&& (df.α .== 0.15), : ]

#RS retailer hand

rsr = []
for hr in 3:0.5:4
    for hc in 0.5:0.5:1.5
        for α in 0.1:0.05:0.15
            for c in 1:0.5:2
                for v in 0.2:0.1:0.4
                    for r in 0.5:0.5:1
                        rrrR = solveRSRaddDecentral(;hr, hc, α, c, v, r)
                        res = (;hr, hc, α, c, v, r, rrrR...)
                        push!(rsr, res)
                   end 
               end
           end 
        end
    end
end
dfrsr = DataFrame(rsr)

using CSV
CSV.write("dfrsr.csv", dfrsr)

# CRS cost-reveue-sharing retailer hand

crsr = []
for hm in 3:0.5:4
    for hc in 0.5:0.5:1.5
        for α in 0.1:0.05:0.15
            for c in 1:0.5:2
                for v in 0.2:0.1:0.4
                    for r in 0.5:0.5:1
                        cccM = solveCRSMaddDecentral(;hm, hc, α, c, v, r)
                        res = (;hm, hc, α, c, v, r, cccM...)
                        push!(crsr, res)
                    end 
                end
             end
        end
    end
end
dfcrsr = DataFrame(crsr)

using CSV
CSV.write("dfcrsr.csv", dfcrsm)


using Plots
d = df[ (df.hc .== 2.0) .&& (df.α .== 0.15), : ]


