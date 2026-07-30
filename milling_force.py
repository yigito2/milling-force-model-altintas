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


def calc_segment_force(
    phid,
    phiu,
    kbeta,
    ft,
    Ktc,
    Krc,
    Kac,
    Kte,
    Kre,
    Kae,
):
    Ftjd = -(1 / kbeta) * (
        -Ktc * ft * np.cos(phid)
        + Kte * phid
    )

    Ftju = -(1 / kbeta) * (
        -Ktc * ft * np.cos(phiu)
        + Kte * phiu
    )

    Ft = Ftju - Ftjd

    Fxjd = -(1 / kbeta) * (
        (-ft / 4)
        * (
            -Ktc * np.cos(2 * phid)
            + Krc * ((2 * phid) - np.sin(2 * phid))
        )
        - Kte * np.sin(phid)
        + Kre * np.cos(phid)
    )

    Fxju = -(1 / kbeta) * (
        (-ft / 4)
        * (
            -Ktc * np.cos(2 * phiu)
            + Krc * ((2 * phiu) - np.sin(2 * phiu))
        )
        - Kte * np.sin(phiu)
        + Kre * np.cos(phiu)
    )

    Fx = Fxju - Fxjd

    Fyjd = -(1 / kbeta) * (
        (ft / 4)
        * (
            Ktc * ((2 * phid) - np.sin(2 * phid))
            + Krc * np.cos(2 * phid)
        )
        - Kte * np.cos(phid)
        - Kre * np.sin(phid)
    )

    Fyju = -(1 / kbeta) * (
        (ft / 4)
        * (
            Ktc * ((2 * phiu) - np.sin(2 * phiu))
            + Krc * np.cos(2 * phiu)
        )
        - Kte * np.cos(phiu)
        - Kre * np.sin(phiu)
    )

    Fy = Fyju - Fyjd

    Fzjd = -(1 / kbeta) * (
        ft * Kac * np.cos(phid)
        - Kae * phid
    )

    Fzju = -(1 / kbeta) * (
        ft * Kac * np.cos(phiu)
        - Kae * phiu
    )

    Fz = Fzju - Fzjd

    return Ft, Fx, Fy, Fz


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

    alphan = atand(tand(alphar) * cosd(beta))

    etac = 0.0

    for eta in np.arange(1, 180, 0.01):
        A = r * cosd(alphan) + cosd(beta) * tand(betaa)
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

    rt = r * (cosd(etac) / cosd(beta))

    phin = atand(
        rt * cosd(alphar)
        /
        (1 - rt * sind(alphar))
    )

    betan = atand(
        tand(betaa) * cosd(etac)
    )

    C_force = np.sqrt(
        (cosd(phin + betan - alphan) ** 2)
        + (tand(etac) ** 2)
        * (sind(betan) ** 2)
    )

    Ktc_theoretical = (
        tau
        * (
            cosd(betan - alphan)
            + tand(beta) * tand(etac) * sind(betan)
        )
    ) / (sind(phin) * C_force)

    Krc_theoretical = (
        tau * sind(betan - alphan)
    ) / (
        sind(phin) * cosd(beta) * C_force
    )

    Kac_theoretical = (
        tau
        * (
            cosd(betan - alphan) * tand(beta)
            - tand(etac) * sind(betan)
        )
    ) / (
        sind(phin) * C_force
    )

    print(f"betaa = {betaa:.3f}")
    print(f"r = {r:.3f}")
    print(f"alphan = {alphan:.3f}")
    print(f"etac = {etac:.3f}")
    print(f"rt = {rt:.3f}")
    print(f"phin = {phin:.3f}")
    print(f"betan = {betan:.3f}")

    print(f"Ktc theoretical = {Ktc_theoretical:.2f}")
    print(f"Krc theoretical = {Krc_theoretical:.2f}")
    print(f"Kac theoretical = {Kac_theoretical:.2f}")

    # Calibrated cutting force coefficients from original MATLAB code
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
        phist = 0.0
        phiex = np.arccos((R - b) / R)

    elif MS == 2:
        phist = np.pi - np.arccos((R - b) / R)
        phiex = np.pi

    else:
        raise ValueError("MS must be 1 for up milling or 2 for down milling.")

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

    Condition = np.zeros((360, N, 6))

    phip = 2 * np.pi / N

    for phi in range(360):
        phi_deg = phi + 1
        phi1[phi] = np.deg2rad(phi_deg)

    print("\nAngular positions generated.")
    print("First angle:", phi1[0])
    print("Last angle :", phi1[-1])

    condition_counter = np.zeros(6, dtype=int)

    for phi in range(360):

        for j in range(N):

            phij_raw = phi1[phi] + j * phip
            phij_current = phij_raw
            phij0 = phij_raw

            if phij_current >= 2 * np.pi:
                phij_current -= 2 * np.pi
                phij0 -= 2 * np.pi

            phij[phi, j] = phij_current

            force_is_calculated = False

            # Condition 1: Out
            if (
                phij0 > phiex
                and phij_current > (phiex + Z * kbeta)
            ):
                Condition[phi, j, 0] = 0
                condition_counter[0] += 1

            # Condition 2: In
            if (
                phij_current >= phiex
                and phij_current <= (phiex + Z * kbeta)
            ):
                Zd = (phij_current - phiex) / kbeta
                Zu = Z

                Condition[phi, j, 1] = 1
                condition_counter[1] += 1
                force_is_calculated = True

            # Condition 3: In
            elif (
                phij_current >= phiex
                and phij_current <= (phist + Z * kbeta)
            ):
                Zd = (phij_current - phiex) / kbeta
                Zu = (phij_current - phist) / kbeta

                Condition[phi, j, 2] = 1
                condition_counter[2] += 1
                force_is_calculated = True

            elif (
                phij_current >= phist
                and phij_current <= phiex
            ):

                # Condition 4: In
                if (
                    phij_current >= (phist + Z * kbeta)
                    and phij_current <= (phiex + Z * kbeta)
                ):
                    Zd = 0.0
                    Zu = Z

                    Condition[phi, j, 3] = 1
                    condition_counter[3] += 1
                    force_is_calculated = True

                # Condition 5: In
                elif phij_current <= (phist + Z * kbeta):
                    Zd = 0.0
                    Zu = (phij_current - phist) / kbeta

                    Condition[phi, j, 4] = 1
                    condition_counter[4] += 1
                    force_is_calculated = True

            # Condition 6: Out
            if (
                phij_current < phist
                and phij_current < (phist + Z * kbeta)
            ):
                Ftj[phi, j] = 0.0
                Fxj[phi, j] = 0.0
                Fyj[phi, j] = 0.0
                Fzj[phi, j] = 0.0

                Condition[phi, j, 5] = 1
                condition_counter[5] += 1

                force_is_calculated = False

            if force_is_calculated:
                phid = phij_current - Zd * kbeta
                phiu = phij_current - Zu * kbeta

                Ft_loc, Fx_loc, Fy_loc, Fz_loc = calc_segment_force(
                    phid=phid,
                    phiu=phiu,
                    kbeta=kbeta,
                    ft=ft,
                    Ktc=Ktc,
                    Krc=Krc,
                    Kac=Kac,
                    Kte=Kte,
                    Kre=Kre,
                    Kae=Kae,
                )

                Ftj[phi, j] = Ft_loc
                Fxj[phi, j] = Fx_loc
                Fyj[phi, j] = Fy_loc
                Fzj[phi, j] = Fz_loc

            Ft[phi] += Ftj[phi, j]
            Fx[phi] += Fxj[phi, j]
            Fy[phi] += Fyj[phi, j]
            Fz[phi] += Fzj[phi, j]

            T[phi] = R * Ft[phi]
            P[phi] = V * Ft[phi]

    print("\nCondition counts")
    print("Condition 1:", condition_counter[0])
    print("Condition 2:", condition_counter[1])
    print("Condition 3:", condition_counter[2])
    print("Condition 4:", condition_counter[3])
    print("Condition 5:", condition_counter[4])
    print("Condition 6:", condition_counter[5])

    Fx_ave_st = (
        N * Z * ft / (8 * np.pi)
        * (
            Ktc * np.cos(2 * phist)
            - Krc * (2 * phist - np.sin(2 * phist))
        )
        + N * Z / (2 * np.pi)
        * (
            -Kte * np.sin(phist)
            + Kre * np.cos(phist)
        )
    )

    Fx_ave_ex = (
        N * Z * ft / (8 * np.pi)
        * (
            Ktc * np.cos(2 * phiex)
            - Krc * (2 * phiex - np.sin(2 * phiex))
        )
        + N * Z / (2 * np.pi)
        * (
            -Kte * np.sin(phiex)
            + Kre * np.cos(phiex)
        )
    )

    Fx_ave = Fx_ave_ex - Fx_ave_st

    Fy_ave_st = (
        N * Z * ft / (8 * np.pi)
        * (
            Ktc * (2 * phist - np.sin(2 * phist))
            + Krc * np.cos(2 * phist)
        )
        - N * Z / (2 * np.pi)
        * (
            Kte * np.cos(phist)
            + Kre * np.sin(phist)
        )
    )

    Fy_ave_ex = (
        N * Z * ft / (8 * np.pi)
        * (
            Ktc * (2 * phiex - np.sin(2 * phiex))
            + Krc * np.cos(2 * phiex)
        )
        - N * Z / (2 * np.pi)
        * (
            Kte * np.cos(phiex)
            + Kre * np.sin(phiex)
        )
    )

    Fy_ave = Fy_ave_ex - Fy_ave_st

    Fz_ave_st = (
        N * Z / (2 * np.pi)
        * (
            -Kac * ft * np.cos(phist)
            + Kae * phist
        )
    )

    Fz_ave_ex = (
        N * Z / (2 * np.pi)
        * (
            -Kac * ft * np.cos(phiex)
            + Kae * phiex
        )
    )

    Fz_ave = Fz_ave_ex - Fz_ave_st

    print("\nForces calculated")
    print("Fx mean numerical =", np.mean(Fx))
    print("Fy mean numerical =", np.mean(Fy))
    print("Fz mean numerical =", np.mean(Fz))

    print("\nAverage forces from analytical equations")
    print("Fx_ave analytical =", Fx_ave)
    print("Fy_ave analytical =", Fy_ave)
    print("Fz_ave analytical =", Fz_ave)

    z = np.zeros((360, N))

    for j in range(N):

        for phi_d in range(360):

            if MS == 1:
                z[phi_d, j] = (
                    phi1[phi_d] + j * phip
                ) / (
                    np.tan(np.radians(beta)) / (2 * R)
                )

            elif MS == 2:
                z[phi_d, j] = (
                    -np.pi + phi1[phi_d] + j * phip
                ) / (
                    np.tan(np.radians(beta)) / (2 * R)
                )

            if z[phi_d, j] > Z or z[phi_d, j] < 0:
                z[phi_d, j] = 0.0

    Mx = np.sum(Fyj * (L - z), axis=1) / 1000.0
    My = np.sum(Fxj * (L - z), axis=1) / 1000.0

    Mx_ave = np.mean(Mx)

    print("\nMoment results")
    print("Mean Mx =", Mx_ave)
    print("Mean My =", np.mean(My))

    plt.figure(figsize=(10, 6))
    plt.plot(phi1, Fx, label="Fx")
    plt.plot(phi1, Fy, label="Fy")
    plt.plot(phi1, Fz, label="Fz")
    plt.xlabel("Rotation (rad)")
    plt.ylabel("Cutting Forces (N)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(phi1, Mx, label="Mx")
    plt.plot(phi1, My, label="My")
    plt.xlabel("Rotation (rad)")
    plt.ylabel("Moment (N.m)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return Mx, My, Mx_ave