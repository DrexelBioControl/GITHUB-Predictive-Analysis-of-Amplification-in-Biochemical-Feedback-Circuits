# ============================================================
# MODEL 3: misfolded-gate fuel leak
#
# Biological assumption:
#   A fraction of gate transcripts forms a misfolded gate, mG.
#   The misfolded gate cannot directly interact with the reporter.
#   Fuel reacts with mG to release OUT:
#
#       mG + F -> OUT + waste
#
# Model name used in configs:
#   model3_misfoldedgate_fuel
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
    "mG",
]

ROL_INDEX = STATE_NAMES.index("ROL")
output_index = ROL_INDEX


def initial_conditions(d, params):

    return [
        0,          # uRSDg
        0,          # RSDg
        0,          # IN
        0,          # OUT
        500,        # DRL
        0,          # ROL
        0,          # I_RSDg
        0,          # F
        0,          # F_RSDg
        0,          # mG
    ]


def rhs(t, x, RSD_temp, IN_temp, F_temp, params):

    # Read shared parameters from dictionary
    k_txn = params["k_txn"]
    ksd = params["ksd"]
    kfsd = params["kfsd"]
    krev = params["krev"]
    kf_rep = params["kf_rep"]
    kRz = params["kRz"]
    basal_frac = params["basal_frac"]
   
    # Constrain reverse strand displacement between F_RSDg and IN to equal forward strand displacement
    kback = ksd

    # Model 3 parameters
    misfold_frac = params["misfold_frac"]
    kmG_F = params["kmG_F"]

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
    mG     = x[9]

    # Misfolded-gate fuel reaction
    v_mG_F = kmG_F * mG * F

    # ODEs

    # Only the correctly folded fraction enters the normal gate pathway.
    duRSDg = (
        (1.0 - misfold_frac) * k_txn * RSD_temp
        - kRz * uRSDg
    )

    dRSDg = (
        kRz * uRSDg
        - ksd * RSDg * IN
        + krev * I_RSDg * OUT
    )

    dIN = (
        k_txn * IN_temp
        - ksd * RSDg * IN
        + krev * I_RSDg * OUT
        + kfsd * I_RSDg * F
    )

    dOUT = (
        ksd * RSDg * IN
        - kf_rep * OUT * DRL
        - krev * I_RSDg * OUT
        + basal_frac * k_txn * RSD_temp
        + v_mG_F
    )

    dDRL = (
        -kf_rep * OUT * DRL
    )

    dROL = (
        kf_rep * OUT * DRL
    )

    dI_RSDg = (
        ksd * RSDg * IN
        - krev * I_RSDg * OUT
        - kfsd * I_RSDg * F
    )

    # Fuel is consumed both by the intended regeneration reaction
    # and by the new misfolded-gate leak reaction.
    dF = (
        k_txn * F_temp
        - kfsd * I_RSDg * F
        + kback * F_RSDg * IN
        - v_mG_F
    )

    dF_RSDg = (
        kfsd * I_RSDg * F
        - kback * F_RSDg * IN
    )

    # Misfolded gate is produced from a fraction of gate transcription
    # and consumed when it reacts with fuel.
    dmG = (
        misfold_frac * k_txn * RSD_temp
        - v_mG_F
    )

    return [
        duRSDg,
        dRSDg,
        dIN,
        dOUT,
        dDRL,
        dROL,
        dI_RSDg,
        dF,
        dF_RSDg,
        dmG,
    ]