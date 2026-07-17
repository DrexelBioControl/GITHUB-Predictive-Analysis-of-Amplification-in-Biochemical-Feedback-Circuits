# ============================================================
# MODEL 2: inactive gate + fuel leak
#
# Biological assumption:
#   Leak occurs through unfolded/inactive gate reacting with fuel:
#       fRSDg + F -> OUT-like product
#
# Model name used in configs:
#   model2_inactive_gate_fuel
# ============================================================

STATE_NAMES = [
    "fRSDg",
    "uRSDg",
    "RSDg",
    "IN",
    "OUT",
    "DRL",
    "ROL",
    "I_RSDg",
    "F",
    "F_RSDg",
]

ROL_INDEX = STATE_NAMES.index("ROL")
output_index = ROL_INDEX


def initial_conditions(d, params):

    return [
        0,          # fRSDg
        0,          # uRSDg
        0,          # RSDg
        0,          # IN
        0,          # OUT
        500,        # DRL
        0,          # ROL
        0,          # I_RSDg
        0,          # F
        0           # F_RSDg
    ]


def rhs(t, x, RSD_temp, IN_temp, F_temp, params):

    # Read parameters from dictionary
    k_txn = params["k_txn"]
    ksd = params["ksd"]
    kfsd = params["kfsd"]
    krev = params["krev"]
    kf_rep = params["kf_rep"]
    kRz = params["kRz"]
    kfld = params["kfld"]
    klk = params["klk"]
    basal_frac = params["basal_frac"]
   
    # Constrain reverse strand displacement between F_RSDg and IN to equal forward strand displacement
    kback = ksd

    # Unpack state variables
    fRSDg  = x[0]
    uRSDg  = x[1]
    RSDg   = x[2]
    IN     = x[3]
    OUT    = x[4]
    DRL    = x[5]
    ROL    = x[6]
    I_RSDg = x[7]
    F      = x[8]
    F_RSDg = x[9]

    # ODEs
    dfRSDg  = k_txn*RSD_temp - kfld*fRSDg - klk*fRSDg*F
    duRSDg  = kfld*fRSDg - kRz*uRSDg
    dRSDg   = kRz*uRSDg - ksd*RSDg*IN + krev*I_RSDg*OUT
    dIN     = k_txn*IN_temp - ksd*RSDg*IN + krev*I_RSDg*OUT + kfsd*I_RSDg*F
    dOUT    = ksd*RSDg*IN - kf_rep*OUT*DRL - krev*I_RSDg*OUT + basal_frac*k_txn*RSD_temp + klk*fRSDg*F
    dDRL    = -kf_rep*OUT*DRL
    dROL    =  kf_rep*OUT*DRL
    dI_RSDg = ksd*RSDg*IN - krev*I_RSDg*OUT - kfsd*I_RSDg*F
    dF      = k_txn*F_temp - kfsd*I_RSDg*F + kback*F_RSDg*IN - klk*fRSDg*F
    dF_RSDg = kfsd*I_RSDg*F - kback*F_RSDg*IN + klk*fRSDg*F

    return [dfRSDg, duRSDg, dRSDg, dIN, dOUT, dDRL, dROL, dI_RSDg, dF, dF_RSDg]
