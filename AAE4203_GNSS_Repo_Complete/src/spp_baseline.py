
import csv
import math
import os

# ==========================================
# AAE4203 - Baseline GPS SPP
# Standard Python libraries only
# ==========================================

FOLDER = "outputs/own_inputs"
OUTFILE = "outputs/spp_baseline_diagnostic.csv"

# Initial receiver ECEF estimate (metres)
INITIAL = [-2418000.0, 5385000.0, 2406000.0]

C = 299792458.0
MAX_ITER = 15
TOL = 1e-4


def read_csv(path):
    with open(path, "r", newline="") as f:
        return list(csv.reader(f))


def number(value):
    try:
        x = float(value.strip())
        return x if math.isfinite(x) else None
    except (ValueError, AttributeError):
        return None


def load_matrix(filename):
    rows = read_csv(os.path.join(FOLDER, filename))
    return [[number(v) for v in row] for row in rows]


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def solve_linear(A, b):
    """Solve Ax=b using Gaussian elimination with pivoting."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]

    for i in range(n):
        pivot = max(range(i, n), key=lambda r: abs(M[r][i]))

        if abs(M[pivot][i]) < 1e-12:
            raise ValueError("Singular matrix")

        M[i], M[pivot] = M[pivot], M[i]

        factor = M[i][i]
        M[i] = [v / factor for v in M[i]]

        for j in range(n):
            if j == i:
                continue
            factor = M[j][i]
            M[j] = [
                M[j][k] - factor * M[i][k]
                for k in range(n + 1)
            ]

    return [M[i][-1] for i in range(n)]


def elevation_angle(receiver, satellite):
    """Calculate satellite elevation in local ENU coordinates."""
    x, y, z = receiver
    sx, sy, sz = satellite

    dx, dy, dz = sx-x, sy-y, sz-z
    r = norm([dx, dy, dz])

    lon = math.atan2(y, x)
    p = math.sqrt(x*x + y*y)

    # WGS-84 approximate geodetic latitude
    lat = math.atan2(z, p * (1 - 6.69437999014e-3))

    east = -math.sin(lon)*dx + math.cos(lon)*dy

    north = (
        -math.sin(lat)*math.cos(lon)*dx
        -math.sin(lat)*math.sin(lon)*dy
        +math.cos(lat)*dz
    )

    up = (
        math.cos(lat)*math.cos(lon)*dx
        +math.cos(lat)*math.sin(lon)*dy
        +math.sin(lat)*dz
    )

    elevation = math.asin(max(-1.0, min(1.0, up/r)))
    return math.degrees(elevation)


def solve_epoch(observations, satellite_data,
                clock_data, iono_data, tropo_data,
                initial):

    receiver = initial[:]
    clock_receiver = 0.0

    for iteration in range(MAX_ITER):

        H = []
        residual = []
        weights = []

        for sat in range(len(observations)):

            pr = observations[sat]
            xyz = satellite_data[sat]
            cb = clock_data[sat]
            ion = iono_data[sat]
            trop = tropo_data[sat]

            if any(v is None for v in [pr, cb, ion, trop]):
                continue

            if any(v is None for v in xyz):
                continue

            sx, sy, sz = xyz
            satellite = [sx, sy, sz]

            dx = sx - receiver[0]
            dy = sy - receiver[1]
            dz = sz - receiver[2]

            distance = norm([dx, dy, dz])

            if distance == 0:
                continue

            el = elevation_angle(receiver, satellite)

            if el <= 0:
                continue

            # Corrected pseudorange (Tutorial Eq. 2)
            corrected_pr = pr + cb - ion - trop

            # Observation equation (Tutorial Eq. 3)
            y = corrected_pr - distance - clock_receiver

            ex = dx / distance
            ey = dy / distance
            ez = dz / distance

            H.append([-ex, -ey, -ez, 1.0])
            residual.append(y)

            # Elevation-dependent weighting
            weights.append(math.sin(math.radians(el))**2)

        if len(H) < 4:
            return None

        # Normal equations: (H'WH) dx = H'Wy
        N = [[0.0]*4 for _ in range(4)]
        u = [0.0]*4

        for row, y, w in zip(H, residual, weights):
            for i in range(4):
                u[i] += w * row[i] * y
                for j in range(4):
                    N[i][j] += w * row[i] * row[j]

        try:
            correction = solve_linear(N, u)
        except ValueError:
            return None

        for i in range(3):
            receiver[i] += correction[i]

        clock_receiver += correction[3]

        if norm(correction[:3]) < TOL:
            return receiver, clock_receiver, len(H)

    return receiver, clock_receiver, len(H)


def main():

    print("Loading GNSS CSV files...")

    epochs = read_csv(os.path.join(FOLDER, "epochs.csv"))
    satellites = read_csv(os.path.join(FOLDER, "satellites.csv"))

    if len(epochs) < 2 or len(satellites) < 2:
        raise ValueError("Epoch or satellite file is empty")

    # First row contains column names
    epoch_data = epochs[1:]
    satellite_ids = [r[0].strip() for r in satellites[1:]]

    pseudoranges = load_matrix("pseudoranges_meas.csv")
    positions = load_matrix("satellite_positions.csv")
    clocks = load_matrix("satellite_clock_bias.csv")
    iono = load_matrix("ionospheric_delay.csv")
    tropo = load_matrix("tropospheric_delay.csv")

    n_sat = len(satellite_ids)
    n_epoch = len(epoch_data)

    print("Satellites:", n_sat)
    print("Epochs:", n_epoch)

    # Verify matrix dimensions
    for name, matrix in [
        ("Pseudoranges", pseudoranges),
        ("Satellite positions", positions),
        ("Satellite clocks", clocks),
        ("Ionosphere", iono),
        ("Troposphere", tropo)
    ]:
        if len(matrix) != n_sat:
            raise ValueError(name + ": satellite row mismatch")

    for i in range(n_sat):
        if len(pseudoranges[i]) != n_epoch:
            raise ValueError("Pseudorange epoch mismatch")

        if len(positions[i]) != 3*n_epoch:
            raise ValueError("Satellite position must have 3 columns per epoch")

        for matrix in [clocks, iono, tropo]:
            if len(matrix[i]) != n_epoch:
                raise ValueError("Correction matrix epoch mismatch")

    print("CSV dimensions verified.")
    print("Starting GPS-only SPP...")

    results = []
    initial = INITIAL[:]

    for k in range(n_epoch):

        obs = [pseudoranges[i][k] for i in range(n_sat)]

        xyz = [
            positions[i][3*k:3*k+3]
            for i in range(n_sat)
        ]

        cb = [clocks[i][k] for i in range(n_sat)]
        ion = [iono[i][k] for i in range(n_sat)]
        trop = [tropo[i][k] for i in range(n_sat)]

        solution = solve_epoch(
            obs, xyz, cb, ion, trop, initial
        )

        if solution is None:
            continue

        pos, receiver_clock, used = solution

        initial = pos[:]

        gps_week = int(float(epoch_data[k][0]))
        gps_tow = float(epoch_data[k][1])

        results.append([
            gps_week,
            gps_tow,
            pos[0],
            pos[1],
            pos[2],
            receiver_clock,
            used
        ])

        if (k+1) % 50 == 0:
            print("Processed:", k+1, "/", n_epoch)

    os.makedirs("outputs", exist_ok=True)

    with open(OUTFILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gps_week", "gps_tow_s",
            "x_ecef_m", "y_ecef_m", "z_ecef_m",
            "receiver_clock_bias_m", "num_satellites"
        ])
        writer.writerows(results)

    print("\n========== SPP COMPLETE ==========")
    print("Total epochs:", n_epoch)
    print("Successful solutions:", len(results))
    print("Failed epochs:", n_epoch-len(results))
    print("Saved:", OUTFILE)

    if results:
        print("\nFirst 5 solutions:")
        for row in results[:5]:
            print(row)


if __name__ == "__main__":
    main()

