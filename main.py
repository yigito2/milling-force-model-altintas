from milling_force import milling_force

Mx, My, Mx_ave = milling_force(
    Z=5.08,
    b=9.05,
    ft=0.05,
    V=30,
    MS=1,
    R=9.05,
    N=4,
    beta=30,
    alphar=12,
    tau=613,
    Kte=24,
    Kre=42,
    Kae=2,
    L=150,
)

print("\nReturned Mx_ave =", Mx_ave)