#!/usr/bin/env python3
"""
analyze_mo_comparison.py
Análisis estadístico de la comparativa multiobjetivo.

Uso:
    python3 scripts/analyze_mo_comparison.py OUTDIR [--csv] [--plot]

Argumentos:
    OUTDIR   directorio con results_eacpp.csv y results_jmetal.csv
    --csv    exporta tabla resumen a OUTDIR/mo_summary.csv
    --plot   genera boxplots (requiere matplotlib)

Los CSV esperados tienen el formato:
    system,problem,run,time_ms,igd
"""

import sys
import csv
import argparse
import statistics
from pathlib import Path
from collections import defaultdict


def load_csv(path):
    """Carga CSV y devuelve dict: problem → list of {run, time_ms, igd}."""
    data = defaultdict(list)
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            data[row["problem"]].append({
                "run":     int(row["run"]),
                "time_ms": float(row["time_ms"]),
                "igd":     float(row["igd"]),
            })
    return data


def stats(values):
    s = sorted(values)
    n = len(s)
    return {
        "n":      n,
        "mean":   statistics.mean(s),
        "median": statistics.median(s),
        "stdev":  statistics.stdev(s) if n > 1 else 0.0,
        "min":    s[0],
        "max":    s[-1],
        "q1":     s[n // 4],
        "q3":     s[3 * n // 4],
    }


def print_block(label, cpp_vals, jmetal_vals, unit="", fmt=".3f"):
    w = 10
    print(f"  {'':20}  {'ea-cpp':>{w}}  {'jMetal':>{w}}")
    print(f"  {'─'*44}")
    for key in ("median", "mean", "stdev", "min", "max"):
        c = getattr(cpp_vals[key] if isinstance(cpp_vals, dict) else cpp_vals, '__float__', lambda: cpp_vals[key])()
        j = jmetal_vals[key]
        c = cpp_vals[key]
        line = f"  {key:<20}  {c:>{w}{fmt}}  {j:>{w}{fmt}}"
        if unit:
            line += f"  {unit}"
        print(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir")
    parser.add_argument("--csv",  action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    cpp_path    = outdir / "results_eacpp.csv"
    jmetal_path = outdir / "results_jmetal.csv"

    if not cpp_path.exists():
        sys.exit(f"ERROR: no encontrado {cpp_path}")

    cpp_data = load_csv(cpp_path)
    jmetal_data = load_csv(jmetal_path) if jmetal_path.exists() else {}

    problems = sorted(set(list(cpp_data.keys()) + list(jmetal_data.keys())))

    print("=" * 68)
    print("  Comparativa multiobjetivo: ea_cpp_modern vs jMetal 7.4")
    print(f"  Problemas: {', '.join(problems)}")
    print("=" * 68)

    summary_rows = []

    for prob in problems:
        print(f"\n── {prob} {'─'*(60-len(prob))}")
        cpp_runs    = cpp_data.get(prob, [])
        jmetal_runs = jmetal_data.get(prob, [])

        cpp_times  = [r["time_ms"] for r in cpp_runs]
        jmetal_times = [r["time_ms"] for r in jmetal_runs]
        cpp_igds   = [r["igd"] for r in cpp_runs]
        jmetal_igds = [r["igd"] for r in jmetal_runs]

        if cpp_times:
            cs = stats(cpp_times)
            print(f"  Tiempo (ms)  ea-cpp   : "
                  f"med={cs['median']:8.1f}  mean={cs['mean']:8.1f}  "
                  f"σ={cs['stdev']:6.1f}  [{cs['min']:.1f}, {cs['max']:.1f}]")
        if jmetal_times:
            js = stats(jmetal_times)
            print(f"  Tiempo (ms)  jMetal   : "
                  f"med={js['median']:8.1f}  mean={js['mean']:8.1f}  "
                  f"σ={js['stdev']:6.1f}  [{js['min']:.1f}, {js['max']:.1f}]")
        if cpp_times and jmetal_times:
            speedup_mean   = js['mean']   / cs['mean']
            speedup_median = js['median'] / cs['median']
            print(f"  Speedup      (media)  : {speedup_mean:.3f}×")
            print(f"  Speedup      (mediana): {speedup_median:.3f}×")
        else:
            speedup_mean = speedup_median = float("nan")

        print()
        if cpp_igds:
            ci = stats(cpp_igds)
            print(f"  IGD          ea-cpp   : "
                  f"med={ci['median']:.6f}  mean={ci['mean']:.6f}  "
                  f"σ={ci['stdev']:.6f}  [{ci['min']:.6f}, {ci['max']:.6f}]")
        if jmetal_igds:
            ji = stats(jmetal_igds)
            print(f"  IGD          jMetal   : "
                  f"med={ji['median']:.6f}  mean={ji['mean']:.6f}  "
                  f"σ={ji['stdev']:.6f}  [{ji['min']:.6f}, {ji['max']:.6f}]")

        # Test de Wilcoxon sobre IGD
        if cpp_igds and jmetal_igds and len(cpp_igds) == len(jmetal_igds):
            try:
                from scipy import stats as sp
                diffs = [c - j for c, j in zip(cpp_igds, jmetal_igds)]
                if all(d == 0.0 for d in diffs):
                    print("  Wilcoxon IGD : todas diferencias = 0")
                else:
                    stat, pval = sp.wilcoxon(cpp_igds, jmetal_igds,
                                             alternative="two-sided", zero_method="zsplit")
                    sig = "SIGNIFICATIVA (p≤0.05)" if pval <= 0.05 else "no significativa (p>0.05)"
                    print(f"  Wilcoxon IGD : W={stat:.1f}  p={pval:.4f}  → {sig}")
            except ImportError:
                pass

        summary_rows.append({
            "problem":        prob,
            "n_cpp":          len(cpp_runs),
            "n_jmetal":       len(jmetal_runs),
            "cpp_time_med":   statistics.median(cpp_times)   if cpp_times   else float("nan"),
            "jmetal_time_med":statistics.median(jmetal_times) if jmetal_times else float("nan"),
            "speedup_median": speedup_median,
            "speedup_mean":   speedup_mean,
            "cpp_igd_med":    statistics.median(cpp_igds)    if cpp_igds    else float("nan"),
            "jmetal_igd_med": statistics.median(jmetal_igds) if jmetal_igds else float("nan"),
        })

    # ── Resumen agregado de speedup ──────────────────────────────────────────
    valid_speedups = [r["speedup_median"] for r in summary_rows
                      if not __import__("math").isnan(r["speedup_median"])]
    if valid_speedups:
        print(f"\n{'═'*68}")
        print(f"  Speedup mediano promedio sobre {len(valid_speedups)} problemas: "
              f"{statistics.mean(valid_speedups):.3f}×")
        print(f"{'═'*68}")

    # ── CSV ─────────────────────────────────────────────────────────────────
    if args.csv:
        csv_out = outdir / "mo_summary.csv"
        with open(csv_out, "w", newline="") as f:
            fields = list(summary_rows[0].keys()) if summary_rows else []
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(summary_rows)
        print(f"\n  Resumen exportado → {csv_out}")

    # ── Gráficas ─────────────────────────────────────────────────────────────
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            n_prob = len(problems)
            fig, axes = plt.subplots(2, n_prob, figsize=(4 * n_prob, 8))
            fig.suptitle("NSGA-II ea_cpp_modern vs jMetal 7.4\n"
                         "pop=100, max_evals=25000, SBX+PM (jMetal protocol)")

            for col, prob in enumerate(problems):
                cpp_runs    = cpp_data.get(prob, [])
                jmetal_runs = jmetal_data.get(prob, [])

                cpp_times   = [r["time_ms"] for r in cpp_runs]
                jmetal_times = [r["time_ms"] for r in jmetal_runs]
                cpp_igds    = [r["igd"]     for r in cpp_runs]
                jmetal_igds  = [r["igd"]     for r in jmetal_runs]

                # Boxplot tiempos
                ax = axes[0, col]
                data_t = [d for d in [cpp_times, jmetal_times] if d]
                labels_t = [l for l, d in [("ea-cpp", cpp_times), ("jMetal", jmetal_times)] if d]
                ax.boxplot(data_t, labels=labels_t, widths=0.5)
                ax.set_title(f"{prob} — Tiempo (ms)")
                ax.set_ylabel("ms" if col == 0 else "")
                ax.grid(axis="y", alpha=0.4)

                # Boxplot IGD
                ax = axes[1, col]
                data_i = [d for d in [cpp_igds, jmetal_igds] if d]
                labels_i = [l for l, d in [("ea-cpp", cpp_igds), ("jMetal", jmetal_igds)] if d]
                ax.boxplot(data_i, labels=labels_i, widths=0.5)
                ax.set_title(f"{prob} — IGD")
                ax.set_ylabel("IGD" if col == 0 else "")
                ax.grid(axis="y", alpha=0.4)

            plt.tight_layout()
            plot_path = outdir / "mo_comparison_plot.png"
            plt.savefig(plot_path, dpi=150)
            print(f"  Gráfica guardada → {plot_path}")
            plt.show()
        except ImportError:
            print("  (matplotlib no disponible — pip install matplotlib)")


if __name__ == "__main__":
    main()
