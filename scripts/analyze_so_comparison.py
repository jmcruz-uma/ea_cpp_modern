#!/usr/bin/env python3
"""
analyze_so_comparison.py
Análisis estadístico de la comparativa monobjetivo.

Uso:
    python3 scripts/analyze_so_comparison.py OUTDIR [--csv] [--plot]

Argumentos:
    OUTDIR   directorio con los JSONs de salida de run_so_comparison.sh
    --csv    exporta resultados a OUTDIR/comparison_results.csv
    --plot   genera gráficas (requiere matplotlib)

Los JSONs esperados en OUTDIR:
    rastrigin-20-1-stats.json         (ea_cpp_original)
    modern-rastrigin-20-1-stats.json  (ea_cpp_modern)
"""

import json
import sys
import statistics
import argparse
from pathlib import Path


def load_runs(path):
    with open(path) as f:
        return json.load(f)


def final_best(run):
    return run["rundata"][0]["idata"]["best"][-1]


def best_at(run, checkpoint):
    evals = run["rundata"][0]["idata"]["evals"]
    best  = run["rundata"][0]["idata"]["best"]
    idx = max((i for i, e in enumerate(evals) if e <= checkpoint), default=None)
    return best[idx] if idx is not None else None


def summary_stats(values, label, width=14):
    s = sorted(values)
    n = len(s)
    q1 = s[n // 4]
    q3 = s[3 * n // 4]
    hits = sum(1 for x in s if x < 1e-6)
    print(f"  {label:<20} "
          f"med={statistics.median(s):>{width}.6f}  "
          f"mean={statistics.mean(s):>{width}.6f}  "
          f"Q1={q1:>{width}.6f}  "
          f"Q3={q3:>{width}.6f}  "
          f"max={max(s):>{width}.6f}  "
          f"hits≈0={hits}/{n}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="Directorio con los JSONs de salida")
    parser.add_argument("--csv",  action="store_true", help="Exportar CSV")
    parser.add_argument("--plot", action="store_true", help="Generar gráficas")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    orig_path   = outdir / "rastrigin-20-1-stats.json"
    modern_path = outdir / "modern-rastrigin-20-1-stats.json"

    for p in [orig_path, modern_path]:
        if not p.exists():
            print(f"ERROR: no encontrado {p}")
            sys.exit(1)

    orig   = load_runs(orig_path)
    modern = load_runs(modern_path)
    n_orig   = len(orig)
    n_modern = len(modern)

    print("=" * 72)
    print("  Comparativa monobjetivo: ea_cpp_original vs ea_cpp_modern")
    print("=" * 72)
    print(f"  Runs: original={n_orig}  moderno={n_modern}")
    print()

    # ── Tiempos ─────────────────────────────────────────────────────────────
    t_o = [r["time"] for r in orig]
    t_m = [r["time"] for r in modern]
    print("── Tiempo por run (segundos) ────────────────────────────────────────")
    summary_stats(t_o, "ea_cpp_original")
    summary_stats(t_m, "ea_cpp_modern")
    speedup = statistics.mean(t_o) / statistics.mean(t_m)
    speedup_med = statistics.median(t_o) / statistics.median(t_m)
    print(f"\n  Speedup (media)   : {speedup:.3f}×")
    print(f"  Speedup (mediana) : {speedup_med:.3f}×")
    print()

    # ── Fitness final ────────────────────────────────────────────────────────
    b_o = [final_best(r) for r in orig]
    b_m = [final_best(r) for r in modern]
    print("── Fitness final (minimización, 0 = óptimo global) ──────────────────")
    summary_stats(b_o, "ea_cpp_original")
    summary_stats(b_m, "ea_cpp_modern")
    print()

    # ── Convergencia en hitos ────────────────────────────────────────────────
    # Inferir maxevals del JSON
    max_evals = orig[0]["rundata"][0]["idata"]["evals"][-1]
    checkpoints = []
    for scale in [100, 1_000, 10_000, 100_000, 500_000, 1_000_000, 10_000_000]:
        if scale <= max_evals:
            checkpoints.append(scale)
    if max_evals not in checkpoints:
        checkpoints.append(max_evals)

    print("── Fitness medio por hito de evaluaciones ───────────────────────────")
    print(f"  {'Evals':>12}  {'Original':>14}  {'Moderno':>14}  {'Ratio M/O':>10}")
    print("  " + "-" * 54)

    rows = []
    for cp in checkpoints:
        vo = [v for r in orig   if (v := best_at(r, cp)) is not None]
        vm = [v for r in modern if (v := best_at(r, cp)) is not None]
        if not vo or not vm:
            continue
        mo, mm = statistics.mean(vo), statistics.mean(vm)
        if mo > 1e-12:
            ratio_s = f"{mm/mo:.3f}"
        elif mm < 1e-10:
            ratio_s = "→0"
        else:
            ratio_s = "∞"
        print(f"  {cp:>12,}  {mo:>14.6f}  {mm:>14.6f}  {ratio_s:>10}")
        rows.append((cp, mo, mm))
    print()

    # ── Test de Wilcoxon (si scipy disponible) ────────────────────────────────
    try:
        from scipy import stats as sp
        if len(b_o) == len(b_m) and len(b_o) >= 10:
            stat, pval = sp.wilcoxon(b_o, b_m, alternative="two-sided")
            print("── Test de Wilcoxon (fitness final, bilateral) ──────────────────────")
            print(f"  W={stat:.1f}  p={pval:.4f}  "
                  f"{'sin diferencia significativa (p>0.05)' if pval > 0.05 else 'DIFERENCIA SIGNIFICATIVA (p≤0.05)'}")
            print()
    except ImportError:
        print("  (scipy no disponible — instala con: pip install scipy)")
        print()

    # ── CSV ──────────────────────────────────────────────────────────────────
    if args.csv:
        import csv
        csv_path = outdir / "comparison_results.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["system", "run", "seed", "time_s", "final_best"])
            for i, r in enumerate(orig):
                w.writerow(["original", i, r.get("seed", i+1), r["time"], final_best(r)])
            for i, r in enumerate(modern):
                w.writerow(["modern", i, r.get("seed", i+1), r["time"], final_best(r)])
        print(f"  CSV exportado → {csv_path}")

    # ── Gráficas ─────────────────────────────────────────────────────────────
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle("Comparativa monobjetivo: ea_cpp_original vs ea_cpp_modern\n"
                         "Rastrigin d=20, (µ=100,λ=99)-ES, BLX-α + Gaussian, 1M evals")

            # Boxplot tiempos
            ax = axes[0]
            ax.boxplot([t_o, t_m], labels=["original", "modern"], widths=0.5)
            ax.set_ylabel("Tiempo por run (s)")
            ax.set_title(f"Tiempo  (speedup media: {speedup:.2f}×)")
            ax.grid(axis="y", alpha=0.4)

            # Curvas de convergencia (mediana por hito)
            ax = axes[1]
            cp_arr = np.array([r[0] for r in rows])
            med_o  = []
            med_m  = []
            for cp, _, _ in rows:
                vo = [v for r in orig   if (v := best_at(r, cp)) is not None]
                vm = [v for r in modern if (v := best_at(r, cp)) is not None]
                med_o.append(statistics.median(vo))
                med_m.append(statistics.median(vm))
            ax.semilogy(cp_arr, med_o, "o-", label="original",  color="steelblue")
            ax.semilogy(cp_arr, med_m, "s-", label="modern",    color="tomato")
            ax.set_xlabel("Evaluaciones")
            ax.set_ylabel("Fitness mediano (log)")
            ax.set_title("Convergencia (mediana, 30 runs)")
            ax.legend()
            ax.grid(alpha=0.4)

            plt.tight_layout()
            plot_path = outdir / "comparison_plot.png"
            plt.savefig(plot_path, dpi=150)
            print(f"  Gráfica guardada → {plot_path}")
            plt.show()
        except ImportError:
            print("  (matplotlib no disponible — instala con: pip install matplotlib)")


if __name__ == "__main__":
    main()
