import numpy as np
import matplotlib.pyplot as plt


def sind(x):
    return np.sin(np.radians(x))


def cosd(x):
    return np.cos(np.radians(x))


def tand(x):
    return np.tan(np.radians(x))


def atand(x):
    return np.degrees(np.arctan(x))

def Ft_integral(phi, ft, Ktc, Kte):
    return -(1.0) * (-Ktc * ft * np.cos(phi) + Kte * phi)


def Fx_integral(phi, ft, Ktc, Krc, Kte, Kre):
    return (
        (-ft / 4.0)
        * (-Ktc * np.cos(2 * phi) + Krc * ((2 * phi) - np.sin(2 * phi)))
        - Kte * np.sin(phi)
        + Kre * np.cos(phi)
    )


def Fy_integral(phi, ft, Ktc, Krc, Kte, Kre):
    return (
        (ft / 4.0)
        * (Ktc * ((2 * phi) - np.sin(2 * phi)) + Krc * np.cos(2 * phi))
        - Kte * np.cos(phi)
        - Kre * np.sin(phi)
    )


def Fz_integral(phi, ft, Kac, Kae):
    return ft * Kac * np.cos(phi) - Kae * phi


def milling_force(
    Z,
    b,
    ft,
    V,
    MS,
    R,
    N,
    beta,
    alphar,
    tau,
    Kte,
    Kre,
    Kae,
    L,
):

    betaa = 19.1 + (0.29 * alphar)

    C0 = 1.755 - 0.028 * alphar * np.pi / 180
    C1 = 0.331 - 0.0082 * alphar * np.pi / 180

    r = C0 * (0.02 ** C1)

    alphan = atand(tand(alphar) * np.cos(np.radians(beta)))

    etac = 0

    for eta in np.arange(1, 180, 0.01):

        A = r * cosd(alphan) + np.cos(np.radians(beta)) * tand(betaa)
        B = tand(betaa) * sind(alphan) * sind(beta)
        C = r * sind(alphan) * tand(betaa)
        D = r * tand(betaa) * tand(beta)
        E = sind(beta) * cosd(alphan)

        delta = (
            A * sind(eta)
            - B * cosd(eta)
            - C * sind(eta) * cosd(eta)
            + D * (cosd(eta) ** 2)
            - E
        )

        if abs(delta) < 0.001:
            etac = eta
            break

    print(f"betaa = {betaa:.3f}")
    print(f"r = {r:.3f}")
    print(f"alphan = {alphan:.3f}")

    rt = r * (cosd(etac) / cosd(beta))

    phin = atand(
        rt * cosd(alphar)
        /
        (1 - rt * sind(alphar))
    )

    betan = atand(
        tand(betaa) * cosd(etac)
    )

    print(f"etac = {etac:.3f}")
    print(f"rt = {rt:.3f}")
    print(f"phin = {phin:.3f}")
    print(f"betan = {betan:.3f}")

    C = np.sqrt(
        (cosd(phin + betan - alphan) ** 2)
        + (np.tan(np.radians(etac)) ** 2)
        * (sind(betan) ** 2)
    )

    Ktc = (
        tau
        * (
            cosd(betan - alphan)
            + tand(beta) * tand(etac) * sind(betan)
        )
    ) / (sind(phin) * C)

    Krc = (
        tau * sind(betan - alphan)
    ) / (
        sind(phin) * cosd(beta) * C
    )

    Kac = (
        tau
        * (
            cosd(betan - alphan) * tand(beta)
            - tand(etac) * sind(betan)
        )
    ) / (
        sind(phin) * C
    )

    print(f"Ktc theoretical = {Ktc:.2f}")
    print(f"Krc theoretical = {Krc:.2f}")
    print(f"Kac theoretical = {Kac:.2f}")

    Ktc = 447.88
    Kte = 70.298

    Krc = 374.21
    Kre = 12.603

    Kac = 97.0
    Kae = 0.0

    print("\nUsing calibrated coefficients")

    print(f"Ktc = {Ktc}")
    print(f"Kte = {Kte}")
    print(f"Krc = {Krc}")
    print(f"Kre = {Kre}")
    print(f"Kac = {Kac}")
    print(f"Kae = {Kae}")

    if MS == 1:

        phist = 0

        phiex = np.arccos((R - b) / R)

    elif MS == 2:

        phist = np.pi - np.arccos((R - b) / R)

        phiex = np.pi

    print(f"\nphist = {phist:.5f}")
    print(f"phiex = {phiex:.5f}")

    kbeta = np.tan(np.radians(beta)) / R

    print(f"kbeta = {kbeta:.6f}")

    Ft = np.zeros(360)
    Fx = np.zeros(360)
    Fy = np.zeros(360)
    Fz = np.zeros(360)

    T = np.zeros(360)
    P = np.zeros(360)

    phi1 = np.zeros(360)

    Ftj = np.zeros((360, N))
    Fxj = np.zeros((360, N))
    Fyj = np.zeros((360, N))
    Fzj = np.zeros((360, N))

    phij = np.zeros((360, N))

    Ftjd = np.zeros((360, N))
    Ftju = np.zeros((360, N))

    Fxjd = np.zeros((360, N))
    Fxju = np.zeros((360, N))

    Fyjd = np.zeros((360, N))
    Fyju = np.zeros((360, N))

    Fzjd = np.zeros((360, N))
    Fzju = np.zeros((360, N))

    Condition = np.zeros((360, N, 6))

    print("\nArrays initialized.")
    print("Fx shape:", Fx.shape)
    print("Fxj shape:", Fxj.shape)
    print("Condition shape:", Condition.shape)

    phip = 2 * np.pi / N

    for phi in range(360):

        phi_deg = phi + 1

        phi1[phi] = np.deg2rad(phi_deg)

    print("\nAngular positions generated.")
    print("First angle:", phi1[0])
    print("Last angle :", phi1[-1])

    for phi in range(360):

        for j in range(N):

            phij[phi, j] = phi1[phi] + j * phip

            if phij[phi, j] >= 2 * np.pi:
                phij[phi, j] -= 2 * np.pi
    print("\nTooth angular positions generated.")

    print("phi1[0] =", phi1[0])

    print("First tooth =", phij[0, 0])
    print("Second tooth =", phij[0, 1])
    print("Third tooth =", phij[0, 2])
    print("Fourth tooth =", phij[0, 3])

    for phi in range(360):

        for j in range(N):

            phi_current = phij[phi, j]

            # OUT OF CUT
            if phi_current < phist or phi_current > phiex:
                continue

            Zd = 0.0
            Zu = Z

            phid = phi_current - Zd * kbeta
            phiu = phi_current - Zu * kbeta

            Ftj[phi, j] = (
                -(1 / kbeta)
                * (
                    Ft_integral(phiu, ft, Ktc, Kte)
                    - Ft_integral(phid, ft, Ktc, Kte)
                )
            )

            Fxj[phi, j] = (
                -(1 / kbeta)
                * (
                    Fx_integral(
                        phiu,
                        ft,
                        Ktc,
                        Krc,
                        Kte,
                        Kre,
                    )
                    -
                    Fx_integral(
                        phid,
                        ft,
                        Ktc,
                        Krc,
                        Kte,
                        Kre,
                    )
                )
            )

            Fyj[phi, j] = (
                -(1 / kbeta)
                * (
                    Fy_integral(
                        phiu,
                        ft,
                        Ktc,
                        Krc,
                        Kte,
                        Kre,
                    )
                    -
                    Fy_integral(
                        phid,
                        ft,
                        Ktc,
                        Krc,
                        Kte,
                        Kre,
                    )
                )
            )

            Fzj[phi, j] = (
                -(1 / kbeta)
                * (
                    Fz_integral(
                        phiu,
                        ft,
                        Kac,
                        Kae,
                    )
                    -
                    Fz_integral(
                        phid,
                        ft,
                        Kac,
                        Kae,
                    )
                )
            )

        Ft[phi] = np.sum(Ftj[phi, :])

        Fx[phi] = np.sum(Fxj[phi, :])
        Fy[phi] = np.sum(Fyj[phi, :])
        Fz[phi] = np.sum(Fzj[phi, :])

        T[phi] = R * Ft[phi]
        P[phi] = V * Ft[phi]

    active_count = 0

    for phi in range(360):

        for j in range(N):

            phij0 = phi1[phi] + j * phip

            if phij0 > phiex and phij[phi, j] > (phiex + Z * kbeta):

                Condition[phi, j, 0] = 1
                active_count += 1
    print("\nCondition 1 test")

    print(
        "Number of teeth outside cut =",
        active_count
    )

    print("\nForces calculated")

    print("Fx mean =", np.mean(Fx))
    print("Fy mean =", np.mean(Fy))
    print("Fz mean =", np.mean(Fz))

    z = np.zeros((360, N))

    for j in range(N):

        for phi_d in range(360):

            if MS == 1:

                z[phi_d, j] = (
                    phi1[phi_d] + j * phip
                ) / (np.tan(np.radians(beta)) / (2 * R))

            else:

                z[phi_d, j] = (
                    -np.pi + phi1[phi_d] + j * phip
                ) / (np.tan(np.radians(beta)) / (2 * R))

            if z[phi_d, j] > Z or z[phi_d, j] < 0:
                z[phi_d, j] = 0

    Mx = np.sum(Fyj * (L - z), axis=1) / 1000.0
    My = np.sum(Fxj * (L - z), axis=1) / 1000.0

    Mx_ave = np.mean(Mx)

    print("\nMoment results")
    print("Mean Mx =", Mx_ave)

    plt.figure(figsize=(10, 6))

    plt.plot(phi1, Fx, label="Fx")
    plt.plot(phi1, Fy, label="Fy")
    plt.plot(phi1, Fz, label="Fz")

    plt.xlabel("Rotation (rad)")
    plt.ylabel("Force (N)")
    plt.legend()
    plt.grid(True)

    plt.show()

    plt.figure(figsize=(10, 6))

    plt.plot(phi1, Fx, label="Fx")
    plt.plot(phi1, Fy, label="Fy")
    plt.plot(phi1, Fz, label="Fz")

    plt.xlabel("Rotation (rad)")
    plt.ylabel("Force (N)")
    plt.legend()
    plt.grid(True)

    plt.show()

    z = np.zeros((360, N))

    for j in range(N):

        for phi_d in range(360):

            if MS == 1:

                z[phi_d, j] = (
                    phi1[phi_d]
                    + j * phip
                ) / (np.tan(np.radians(beta)) / (2 * R))

            else:

                z[phi_d, j] = (
                    -np.pi
                    + phi1[phi_d]
                    + j * phip
                ) / (np.tan(np.radians(beta)) / (2 * R))

            if z[phi_d, j] > Z or z[phi_d, j] < 0:
                z[phi_d, j] = 0
        Mx = np.sum(Fyj * (L - z), axis=1) / 1000.0
    My = np.sum(Fxj * (L - z), axis=1) / 1000.0

    Mx_ave = np.mean(Mx)

    print("\nMoment results")
    print("Mean Mx =", Mx_ave)

    return Mx, My, Mx_ave