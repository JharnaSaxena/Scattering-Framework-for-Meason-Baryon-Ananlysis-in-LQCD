# pipeline_adapter.py
# PSQ total momentum vector (integer)
PSQ_TO_D = {
    "PSQ0": (0, 0, 0),
    "PSQ1": (0, 0, 1),
    "PSQ2": (1, 1, 0),
    "PSQ3": (1, 1, 1),
    "PSQ4": (0, 0, 2),
}

def full_irrep_label(psq: str, irrep: str) -> str:
    #convert (psq, irrep) to morningstar label like G1u(0)
    psq_num = psq.replace("PSQ", "")
    return f"{irrep}({psq_num})"
