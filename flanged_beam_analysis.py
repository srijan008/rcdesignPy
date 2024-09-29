import matplotlib.pyplot as plt
import numpy as np

def input_geometry():
    while True:
        try:
            b_f = float(input("Enter flange width in mm: "))
            d_f = float(input("Enter flange depth in mm: "))
            b_w = float(input("Enter web width in mm: "))
            D = float(input("Enter overall depth (D) in mm: "))
            d = float(input("Enter distance from flange to reinforcement (d) in mm: "))
            n_bars = int(input("Enter number of reinforcement bars: "))
            dia_bars = float(input("Enter diameter of reinforcement bars (in mm): "))
            fck = float(input("Enter concrete strength (fck) in MPa: "))
            fy = float(input("Enter steel yield strength (fy) in MPa: "))
            
            if b_f <= 0 or d_f <= 0 or b_w <= 0 or D <= 0 or d <= 0 or n_bars <= 0 or dia_bars <= 0 or fck <= 0 or fy <= 0:
                raise ValueError("All values must be positive and non-zero.")
            break

        except ValueError as e:
            print(f"Invalid input! {e}. Please re-enter the values.")

    return b_f, d_f, b_w, D, d, n_bars, dia_bars, fck, fy

def draw_cross_section(b_f, d_f, b_w, D, d, n_bars, dia_bars, xu=None):
    fig, ax = plt.subplots()
    ax.set_xlim(0, b_f + 20)
    ax.set_ylim(0, D + 50)

    # Draw flange
    rect_flange = plt.Rectangle((10, D - d_f), b_f, d_f, edgecolor='black', facecolor='lightgray')
    ax.add_patch(rect_flange)

    # Draw web
    rect_web = plt.Rectangle((10 + (b_f - b_w) / 2, 0), b_w, D - d_f, edgecolor='black', facecolor='gray')
    ax.add_patch(rect_web)

    for i in range(n_bars):
        x = 10 + (b_f - b_w) / 2 + (b_w / (n_bars + 1)) * (i + 1)
        y = D - d - dia_bars / 2  
        ax.add_patch(plt.Circle((x, y), dia_bars / 2, color='red'))

    if xu is not None:
        plt.axhline(D - xu, color='blue', linestyle='--', label=f'Neutral Axis (xu = {xu:.2f} mm)')
    
    ax.set_aspect('equal')
    plt.gca().invert_yaxis()
    plt.title("Cross Section of Flanged Beam")
    plt.legend()
    plt.show(block=False)


def flexural_analysis(b_f, d_f, b_w, D, d, n_bars, dia_bars, fck, fy):
    Ast = n_bars * np.pi * (dia_bars / 2) ** 2
    
    Xu = 0.87 * fy * Ast / (0.36 * fck * b_f)
    
    if Xu <= d_f:
        xu = Xu
        Mu = 0.36 * fck * b_w * xu * (d - 0.416 * xu)

    else:
        xu1 = 7 / 3 * d_f 
        Cuf = 0.446 * fck * (b_f - b_w) * d_f  
        Tu = 0.87 * fy * Ast 
        xu = (Tu - Cuf) / (0.36 * fck * b_w)
        xu_temp = 100
        yf = d_f
        if xu1 > xu:
           while True:
               yf = 0.15 * xu_temp + 0.65 * d_f
               xu = (0.87 * fy * Ast - 0.446 * fck * (b_f - b_w) * yf) / (0.36 * fck * b_w)
               if abs(xu - xu_temp) < 0.01:
                   break
               xu_temp = xu
        Mu = 0.36 * fck * b_w * xu * (d - 0.416 * xu) + 0.447 * fck * (b_f - b_w) * yf * (d - yf / 2)

    if fy == 250:
        xu_max_by_d = 0.53
    elif fy == 415:
        xu_max_by_d = 0.48
    elif fy == 500:
        xu_max_by_d = 0.46

    if (xu / d) > xu_max_by_d:

        a = 0.0035
        b = -0.0035 * d
        c = 0.87 * fy * Ast / (2 * 10**5)

        discriminant = b**2 - 4 * a * c
        if discriminant >= 0:
            xu1 = (-b + np.sqrt(discriminant)) / (2 * a)
            xu2 = (-b - np.sqrt(discriminant)) / (2 * a)

            xu = max(xu1, xu2) if xu1 > 0 else xu2

            if xu <= d_f:
                Mu = 0.36 * fck * b_w * xu * (d - 0.416 * xu)
                return xu, Mu
            else:
                xu1 = 7 / 3 * d_f
                Cuf = 0.446 * fck * (b_f - b_w) * d_f
                Tu = 0.87 * fy * Ast
                xu = (Tu - Cuf) / (0.36 * fck * b_w)
                xu_temp = 100
                yf = d_f
                if xu1 > xu:
                    while True:
                        yf = 0.15 * xu_temp + 0.65 * d_f
                        xu = (0.87 * fy * Ast - 0.446 * fck * (b_f - b_w) * yf) / (0.36 * fck * b_w)
                        if abs(xu - xu_temp) < 0.01:
                            break
                        xu_temp = xu
                Mu = 0.36 * fck * b_w * xu * (d - 0.416 * xu) + 0.447 * fck * (b_f - b_w) * yf * (d - yf / 2)

    # let's check limitiing moment capacity

    if d_f / d < 0.2:
        Mu_lim = 0.36 * xu_max_by_d* (1 - 0.42*xu_max_by_d) * fck * b_w * d**2 + 0.45 * fck * (b_f - b_w) * d_f * (d - d_f / 2)
    else:
        Mu_lim = 0.36 * xu_max_by_d* (1 - 0.42*xu_max_by_d) * fck * b_w * d**2 + 0.45 * fck * (b_f - b_w) * yf * (d - yf / 2)
    
    return xu,xu_max_by_d * d, Mu_lim, Mu



def main():
    b_f, d_f, b_w, D, d, n_bars, dia_bars, fck, fy = input_geometry()
    draw_cross_section(b_f, d_f, b_w, D, d, n_bars, dia_bars)
    input("Confirm section? (Press Enter to continue)")

    xu,xu_lim, Mu_lim, Mu = flexural_analysis(b_f, d_f, b_w, D, d, n_bars, dia_bars, fck, fy)
    print(f"Neutral Axis Depth (xu): {xu:.2f} mm")
    print(f"Limiting Neutral Axis Depth (xu): {xu_lim:.2f} mm")
    print(f"Actual  Moment Capacity (Mu): {Mu_lim / 1e6:.2f} kNm")
    print(f"Limiting Moment Capacity (Mu): {Mu / 1e6:.2f} kNm")

    draw_cross_section(b_f, d_f, b_w, D, d, n_bars, dia_bars, xu)
    input("Confirm the neutral axis visualization (Press Enter to continue)")

if __name__ == "__main__":
    main()

            # b_f = 1000  # flange width in mm
            # d_f = 100   # flange depth in mm
            # b_w = 250  # web width in mm
            # D = 600    # overall depth in mm
            # d = 520    # distance from flange to reinforcement in mm
            # n_bars = 6 # number of reinforcement bars
            # dia_bars = 28 # diameter of reinforcement bars in mm
            # fck = 20   # concrete strength in MPa
            # fy = 250   # steel yield strength in MPa
        
         