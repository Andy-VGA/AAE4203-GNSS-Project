#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AAE4203 - GPS-only Single Point Positioning (SPP)
Weighted least squares + Gauss-Newton iteration.

Expected input directory:
    outputs/own_inputs/

Required files:
    epochs.csv
    satellites.csv
    pseudoranges_meas.csv
    satellite_positions.csv
    satellite_clock_bias.csv
    ionospheric_delay.csv
    tropospheric_delay.csv
    group_delay.csv
    elevation_deg.csv

The matrix files are satellite-major: one row per satellite, one column per
epoch; satellite_positions.csv has 3 columns per epoch (X,Y,Z).

All range-like corrections are assumed to be in metres, consistent with the
MATLAB observation model used in this project:
    P_corrected = P + satellite_clock - ionosphere - troposphere - TGD

Every epoch is retained. Failed epochs are written with NaN coordinates and
status="FAILED".
"""

import csv
import math
import os

FOLDER = "outputs/own_inputs"
OUTFILE = "outputs/spp_baseline_diagnostic.csv"
INITIAL = [-2418000.0, 5385000.0, 2406000.0]
MAX_ITER = 15
TOL = 1e-4
MIN_ELEVATION_DEG = 15.0
MIN_SATELLITES = 4
WEIGHT_FLOOR = 1e-6


def read_csv(path):
    with open(path, "r", newline="") as f:
        return list(csv.reader(f))


def number(value):
    try:
        x = float(str(value).strip())
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def load_matrix(filename):
    rows = read_csv(os.path.join(FOLDER, filename))
    return [[number(v) for v in row] for row in rows]


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def solve_linear(A, b):
    """Solve Ax=b by Gaussian elimination with partial pivoting."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for i in range(n):
        pivot = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[pivot][i]) < 1e-12:
            raise ValueError("Singular normal matrix")
        M[i], M[pivot] = M[pivot], M[i]
        pivot_value = M[i][i]
        M[i] = [v / pivot_value for v in M[i]]
        for j in range(n):
            if j == i:
                continue
            factor = M[j][i]
            M[j] = [M[j][k] - factor * M[i][k] for k in range(n + 1)]
    return [M[i][-1] for i in range(n)]


def solve_epoch(observations, satellite_data, clock_data, iono_data,
                tropo_data, tgd_data, elevation_data, initial):
    """Estimate receiver ECEF XYZ and receiver clock bias (metres)."""
    receiver = initial[:]
    receiver_clock = 0.0
    used_last = 0

    for iteration in range(1, MAX_ITER + 1):
        H, residual, weights = [], [], []
        for sat in range(len(observations)):
            pr, xyz = observations[sat], satellite_data[sat]
            cb, ion = clock_data[sat], iono_data[sat]
            trop, tgd, el = tropo_data[sat], tgd_data[sat], elevation_data[sat]
            if any(v is None for v in [pr, cb, ion, trop, tgd, el]):
                continue
            if any(v is None for v in xyz) or el < MIN_ELEVATION_DEG:
                continue

            sx, sy, sz = xyz
            dx, dy, dz = sx - receiver[0], sy - receiver[1], sz - receiver[2]
            distance = norm([dx, dy, dz])
            if not math.isfinite(distance) or distance <= 0:
                continue

            corrected_pr = pr + cb - ion - trop - tgd
            y = corrected_pr - distance - receiver_clock
            ex, ey, ez = dx / distance, dy / distance, dz / distance
            H.append([-ex, -ey, -ez, 1.0])
            residual.append(y)
            weights.append(max(math.sin(math.radians(el)) ** 2, WEIGHT_FLOOR))

        used_last = len(H)
        if used_last < MIN_SATELLITES:
            return None

        N = [[0.0] * 4 for _ in range(4)]
        u = [0.0] * 4
        for row, v, w in zip(H, residual, weights):
            for i in range(4):
                u[i] += w * row[i] * v
                for j in range(4):
                    N[i][j] += w * row[i] * row[j]

        try:
            correction = solve_linear(N, u)
        except ValueError:
            return None

        for i in range(3):
            receiver[i] += correction[i]
        receiver_clock += correction[3]

        if norm(correction[:3]) < TOL:
            return {"position": receiver, "clock": receiver_clock,
                    "used": used_last, "iterations": iteration,
                    "converged": True}

    return {"position": receiver, "clock": receiver_clock, "used": used_last,
            "iterations": MAX_ITER, "converged": False}


def validate_dimensions(n_sat, n_epoch, matrices):
    for name, matrix in matrices.items():
        if len(matrix) != n_sat:
            raise ValueError(f"{name}: expected {n_sat} satellite rows, got {len(matrix)}")
    for i in range(n_sat):
        for name in ["pseudoranges", "satellite_clocks", "ionosphere",
                     "troposphere", "group_delay", "elevation"]:
            if len(matrices[name][i]) != n_epoch:
                raise ValueError(f"{name}: satellite row {i + 1} has "
                                 f"{len(matrices[name][i])} epochs; expected {n_epoch}")
        if len(matrices["satellite_positions"][i]) != 3 * n_epoch:
            raise ValueError(f"satellite_positions: row {i + 1} must contain "
                             f"{3 * n_epoch} values; got "
                             f"{len(matrices['satellite_positions'][i])}")


def main():
    print("Loading GNSS CSV files...")
    epochs = read_csv(os.path.join(FOLDER, "epochs.csv"))
    satellites = read_csv(os.path.join(FOLDER, "satellites.csv"))
    if len(epochs) < 2 or len(satellites) < 2:
        raise ValueError("Epoch or satellite file is empty")

    epoch_data = epochs[1:]
    satellite_ids = [row[0].strip() for row in satellites[1:] if row]
    matrices = {
        "pseudoranges": load_matrix("pseudoranges_meas.csv"),
        "satellite_positions": load_matrix("satellite_positions.csv"),
        "satellite_clocks": load_matrix("satellite_clock_bias.csv"),
        "ionosphere": load_matrix("ionospheric_delay.csv"),
        "troposphere": load_matrix("tropospheric_delay.csv"),
        "group_delay": load_matrix("group_delay.csv"),
        "elevation": load_matrix("elevation_deg.csv"),
    }

    n_sat, n_epoch = len(satellite_ids), len(epoch_data)
    print("Satellites:", n_sat)
    print("Epochs:", n_epoch)
    validate_dimensions(n_sat, n_epoch, matrices)
    print("CSV dimensions verified, including TGD and elevation.")

    results = []
    solved_count = converged_count = failed_count = nonconverged_count = 0
    initial = INITIAL[:]

    for k in range(n_epoch):
        obs = [matrices["pseudoranges"][i][k] for i in range(n_sat)]
        xyz = [matrices["satellite_positions"][i][3*k:3*k+3] for i in range(n_sat)]
        cb = [matrices["satellite_clocks"][i][k] for i in range(n_sat)]
        ion = [matrices["ionosphere"][i][k] for i in range(n_sat)]
        trop = [matrices["troposphere"][i][k] for i in range(n_sat)]
        tgd = [matrices["group_delay"][i][k] for i in range(n_sat)]
        elevation = [matrices["elevation"][i][k] for i in range(n_sat)]
        solution = solve_epoch(obs, xyz, cb, ion, trop, tgd, elevation, initial)
        gps_week, gps_tow = int(float(epoch_data[k][0])), float(epoch_data[k][1])

        if solution is None:
            failed_count += 1
            results.append([gps_week, gps_tow, "NaN", "NaN", "NaN", "NaN", 0, 0, "FAILED"])
        else:
            pos, clock = solution["position"], solution["clock"]
            used, iterations = solution["used"], solution["iterations"]
            converged = solution["converged"]
            initial = pos[:]
            solved_count += 1
            if converged:
                converged_count += 1
                status = "CONVERGED"
            else:
                nonconverged_count += 1
                status = "MAX_ITER"
            results.append([gps_week, gps_tow, pos[0], pos[1], pos[2],
                            clock, used, iterations, status])

        if (k + 1) % 50 == 0 or (k + 1) == n_epoch:
            print(f"Processed: {k + 1} / {n_epoch}")

    os.makedirs("outputs", exist_ok=True)
    with open(OUTFILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gps_week", "gps_tow_s", "x_ecef_m", "y_ecef_m",
                         "z_ecef_m", "receiver_clock_bias_m", "num_satellites",
                         "iterations", "status"])
        writer.writerows(results)

    print("\n========== SPP COMPLETE ==========")
    print("Total epochs:", n_epoch)
    print("Solved epochs:", solved_count)
    print("Converged epochs:", converged_count)
    print("Non-converged epochs:", nonconverged_count)
    print("Failed epochs:", failed_count)
    print("Output rows (excluding header):", len(results))
    print("Saved:", OUTFILE)
    if results:
        print("\nFirst 5 output rows:")
        for row in results[:5]:
            print(row)


if __name__ == "__main__":
    main()
