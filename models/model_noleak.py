# ============================================================
# MODEL: no leak
#
# Biological assumption:
#   no leak at all, not even transciptional leak
#
#
# Model name used in configs:
#   model_noleak
# ============================================================

STATE_NAMES = [
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

    DRL_0 = d["DRL_0"]

    return [
        0,          # uRSDg
        0,          # RSDg
        0,          # IN
        0,          # OUT
        DRL_0,      # DRL
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

    # Unpack state variables
    uRSDg  = x[0]
    RSDg   = x[1]
    IN     = x[2]
    OUT    = x[3]
    DRL    = x[4]
    ROL    = x[5]
    I_RSDg = x[6]
    F      = x[7]
    F_RSDg = x[8]

    # ODEs
    duRSDg = k_txn*RSD_temp- kRz*uRSDg
    dRSDg = kRz*uRSDg - ksd*RSDg*IN + krev*I_RSDg*OUT
    dIN = k_txn*IN_temp - ksd*RSDg*IN + krev*I_RSDg*OUT + kfsd*I_RSDg*F
    dOUT = ksd*RSDg*IN - kf_rep*OUT*DRL - krev*I_RSDg*OUT
    dDRL = - kf_rep*OUT*DRL
    dROL = kf_rep*OUT*DRL
    dI_RSDg = ksd*RSDg*IN - krev*I_RSDg*OUT - kfsd*I_RSDg*F
    dF = k_txn*F_temp - kfsd*I_RSDg*F + ksd*F_RSDg*IN
    dF_RSDg = kfsd*I_RSDg*F - ksd*F_RSDg*IN

    return [duRSDg, dRSDg, dIN, dOUT, dDRL, dROL, dI_RSDg, dF, dF_RSDg]