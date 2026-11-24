"""The NumPy twin of TwoBeamThickness.fit (DiffractoGeometry, Sources/DiffractoGeometry/TwoBeamThickness.swift), step for
step, to read what the Swift fit should answer on the committed fixtures and to try a change before porting it.

Run (any venv with NumPy):
    python3 -I -B tools/pins/two-beam-thickness/fit.py [Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-disc.json]
Prints each case's fitted thickness and its error against the simulated one, the extinction distance, the residual, the
fringe line and a second candidate. It writes nothing; the Swift tests hold their own bounds (3 %, two candidates).

The method (Kelly, Jostsons, Blake and Napier, phys. stat. sol. (a) 31, 771 (1975); Allen, Philos. Mag. A 43, 325 (1981)):
the two-beam intensity of the diffracted beam I(s) = sin^2(pi t s_eff) / (xi s_eff)^2, s_eff^2 = (s - c)^2 + 1/xi^2,
times a scale plus a constant background, fitted by least squares over the whole profile: a grid over thickness (ratio
1.002) and 25 extinction distances in the window 0.7 to 1.15 of the two-beam value, scale and background solved
linearly at every point, then Levenberg-Marquardt on (t, xi, c) from the best local minima of the grid. The fringe line
goes through the measured minima nearest the ones the fit predicts.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
WINDOW = (0.7, 1.15)
XI_STEPS = 25
T_RATIO = 1.002
SECOND_MARGIN = 0.05   # a second thickness fits "as well" within 5 % of the best rms ...
SECOND_APART = 0.05    # ... and lies more than 5 % away


def shape(s, t, xi, c):
    se = np.sqrt((s - c) ** 2 + 1 / xi ** 2)
    return np.sin(np.pi * t * se) ** 2 / (xi * se) ** 2


def linear(p, y):
    """Scale and background by least squares (centred); a scale <= 0 is no fit: the constant ȳ."""
    pm, ym = p.mean(), y.mean()
    sxx = float(np.sum((p - pm) ** 2))
    a = float(np.sum((p - pm) * (y - ym))) / sxx if sxx > 0 else 0.0
    if not a > 0:
        return 0.0, ym
    return a, ym - a * pm


def ssr(s, y, t, xi, c):
    p = shape(s, t, xi, c)
    a, b = linear(p, y)
    r = a * p + b - y
    return float(np.sum(r * r))


def lm(s, y, x0, lo, hi):
    """Levenberg-Marquardt with central-difference Jacobian on scaled parameters, clamped to [lo, hi]."""
    x = np.array(x0, float)
    scale = np.array([x0[0], x0[1], s[-1] - s[0]])

    def res(q):
        p = shape(s, q[0], q[1], q[2])
        a, b = linear(p, y)
        return a * p + b - y

    r = res(x)
    f = float(r @ r)
    lam = 1e-3
    for _ in range(200):
        jac = np.empty((len(s), 3))
        for k in range(3):
            h = 1e-7 * scale[k]
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            jac[:, k] = (res(xp) - res(xm)) / (2 * h) * scale[k]
        jtj = jac.T @ jac
        jtr = jac.T @ r
        improved = False
        while lam < 1e12:
            a = jtj + lam * np.diag(np.diag(jtj))
            try:
                step = -np.linalg.solve(a, jtr)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            xn = np.clip(x + step * scale, lo, hi)
            rn = res(xn)
            fn = float(rn @ rn)
            if fn < f:
                done = np.all(np.abs(xn - x) <= 1e-12 * scale)
                x, r, f = xn, rn, fn
                lam = max(lam / 10, 1e-12)
                improved = True
                break
            lam *= 10
        if not improved or done:
            break
    return x, f


def minima(s, y):
    """Local minima of the samples, parabola-interpolated."""
    out = []
    ds = s[1] - s[0]
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] <= y[i + 1]:
            den = y[i - 1] - 2 * y[i] + y[i + 1]
            out.append(s[i] + (0.5 * (y[i - 1] - y[i + 1]) / den * ds if den > 0 else 0.0))
    return np.array(out)


def predicted(t, xi, n):
    v = n * n / (t * t) - 1 / (xi * xi)
    return math.sqrt(v) if v >= 0 else None


def fringe_line(orders, devs):
    x = 1 / np.asarray(orders, float) ** 2
    yy = (np.asarray(devs) / np.asarray(orders, float)) ** 2
    xm, ym = x.mean(), yy.mean()
    slope = float(np.sum((x - xm) * (yy - ym)) / np.sum((x - xm) ** 2))
    icpt = ym - slope * xm
    return 1 / math.sqrt(icpt), (1 / math.sqrt(-slope) if slope < 0 else float("nan"))


def fit(s, y, xi0, trange=(100.0, 5000.0)):
    s = np.asarray(s, float)
    y = np.asarray(y, float)
    w = y - y.min()
    c0 = float(np.sum(s * w) / np.sum(w))
    xis = xi0 * WINDOW[0] * (WINDOW[1] / WINDOW[0]) ** (np.arange(XI_STEPS) / (XI_STEPS - 1))
    nt = int(math.floor(math.log(trange[1] / trange[0]) / math.log(T_RATIO))) + 1
    ts = trange[0] * T_RATIO ** np.arange(nt)
    r = np.full(nt, np.inf)
    bx = np.zeros(nt)
    for i, t in enumerate(ts):
        for xi in xis:
            v = ssr(s, y, t, xi, c0)
            if v < r[i]:
                r[i], bx[i] = v, xi
    loc = [i for i in range(1, nt - 1) if r[i] < r[i - 1] and r[i] <= r[i + 1]]
    loc.sort(key=lambda i: r[i])
    lo = [trange[0] * 0.5, xi0 * WINDOW[0], s[0]]
    hi = [trange[1] * 2, xi0 * WINDOW[1], s[-1]]
    cands = []
    for i in loc[:4]:
        if r[i] > 4 * r[loc[0]] and i != loc[0]:
            break
        x, f = lm(s, y, [ts[i], bx[i], c0], lo, hi)
        cands.append((f, x))
    cands.sort(key=lambda q: q[0])
    best = cands[0]
    rms = lambda f: math.sqrt(f / len(s)) / float(np.max(np.abs(y)))  # noqa: E731
    second = None
    for f, x in cands[1:]:
        if abs(x[0] / best[1][0] - 1) > SECOND_APART and math.sqrt(f) <= (1 + SECOND_MARGIN) * math.sqrt(best[0]):
            second = (x, rms(f))
            break
    t, xi, c = best[1]
    n0 = int(math.floor(t / xi)) + 1
    left = right = 0
    n = n0
    while True:
        sp = predicted(t, xi, n)
        if c + sp > s[-1] and c - sp < s[0]:
            break
        right += c + sp <= s[-1]
        left += c - sp >= s[0]
        n += 1
    # fringe line
    mins = minima(s, y)
    orders, devs = [], []
    n = n0
    while len(mins) > 0:
        sp = predicted(t, xi, n)
        sn = predicted(t, xi, n + 1)
        if c + sp > s[-1] and c - sp < s[0]:
            break
        tol = 0.35 * (sn - sp)
        pair = []
        for sgn in (-1, 1):
            k = int(np.argmin(np.abs(mins - (c + sgn * sp))))
            if abs(mins[k] - (c + sgn * sp)) < tol:
                pair.append(sgn * (mins[k] - c))
        if len(pair) == 2:
            orders.append(n)
            devs.append(0.5 * (pair[0] + pair[1]))
        n += 1
    line = fringe_line(orders, devs) if len(orders) >= 2 else None
    return dict(t=t, xi=xi, c=c, rms=rms(best[0]), left=left, right=right, line=line, orders=orders, second=second,
                cands=[(x[0], rms(f)) for f, x in cands])


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Tests/DiffractoGeometryTests/Fixtures/two-beam-thickness-disc.json"
    doc = json.loads(path.read_text())
    xi0 = float(sys.argv[2]) if len(sys.argv) > 2 else 967.0335034036704  # the kit's two-beam value for Si 220, 200 kV
    for case in doc["cases"]:
        q = fit(case["deviation"], case["intensity"], xi0)
        true = case["thickness"]
        line = f"{q['line'][0] / 10:.2f} nm ({100 * (q['line'][0] / true - 1):+.2f} %)" if q["line"] else "none"
        sec = f"{q['second'][0][0] / 10:.1f} nm (rms {q['second'][1]:.4f})" if q["second"] else "none"
        print(f"{case['name']}: {q['t'] / 10:.2f} nm ({100 * (q['t'] / true - 1):+.2f} %), xi {q['xi'] / 10:.1f} nm, "
              f"rms {q['rms']:.4f}, minima {q['left']}/{q['right']}, fringe line {line}, second {sec}; "
              f"candidates {[(round(a / 10, 1), round(b, 4)) for a, b in q['cands']]}")


if __name__ == "__main__":
    main()
