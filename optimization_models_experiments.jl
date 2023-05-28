# load the file with all the implementations
include("optimization_models.jl")


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
go(ParameterSweep(
    "wpcm",
    solveRSMaddCentral,
    [   :hm => 1:0.5:3,
        :hc => 2:0.5:4,
        :α => 0.05:0.05:0.15,
        :c => 0.5:0.5:2,
        :v => 0.2:0.1:0.4]
))

# model central retailer hand
go(ParameterSweep(
    "wpcr",
    solveRSMaddCentral,
    [   :hm => 3:0.5:5,
        :hc => 1:0.5:2,
        :α => 0.05:0.05:0.15,
        :c => 0.5:0.5:2,
        :v => 0.2:0.1:0.4]
))

# wholesaleprice manufacturer hand
# NIE MA TEGO
"""
go(ParameterSweep(
    "cm",
    wholesalePriceContractM,  # BRAK TAKIEJ FUNKCJI
    [   :hm => 1:0.5:3,
        :hc => 2:0.5:4,
        :α => 0.05:0.05:0.15,
        :c => 0.5:0.5:2,
        :v => 0.2:0.1:0.4]
))
"""


# wholesaleprice retailer hand

# NIE MA TEGO
"""
go(ParameterSweep(
    "cr",
    wholesalePriceContractR,  # BRAK TAKIEJ FUNKCJI
    [   :hr => 3:0.5:5,
        :hc => 1:0.5:2,
        :α => 0.05:0.05:0.15,
        :c => 0.5:0.5:2,
        :v => 0.2:0.1:0.4]
))
"""

# RS reveue-sharing manufacturer hand

go(ParameterSweep(
    "rsm",
    solveRSMaddDecentral, 
    [   :hm => 1:0.5:2,
        :hc => 3:0.5:4,
        :α => 0.1:0.05:0.15,
        :c => 1:0.5:2,
        :v => 0.2:0.1:0.4,
        :r => 0.5:0.5:1
    ]
))


# CRS cost-reveue-sharing manufacturer hand

go(ParameterSweep(
    "crsm",
    solveCRSMaddDecentral, 
    [   :hm => 1:0.5:2,
        :hc => 3:0.5:4,
        :α => 0.1:0.05:0.15,
        :c => 1:0.5:2,
        :v => 0.2:0.1:0.4,
        :r => 0.5:0.5:1
    ]
))


#RS retailer hand
go(ParameterSweep(
    "rsr",
    solveRSRaddDecentral, 
    [   :hr => 3:0.5:4,
        :hc => 0.5:0.5:1.5,
        :α => 0.1:0.05:0.15,
        :c => 1:0.5:2,
        :v => 0.2:0.1:0.4,
        :r => 0.5:0.5:1
    ]
))

# CRS cost-reveue-sharing retailer hand

go(ParameterSweep(
    "crsr",
    solveCRSMaddDecentral, 
    [   :hm => 3:0.5:4,
        :hc => 0.5:0.5:1.5,
        :α => 0.1:0.05:0.15,
        :c => 1:0.5:2,
        :v => 0.2:0.1:0.4,
        :r => 0.5:0.5:1
    ]
))


