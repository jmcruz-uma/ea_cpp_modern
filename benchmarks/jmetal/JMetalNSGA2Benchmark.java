import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.operator.crossover.impl.SBXCrossover;
import org.uma.jmetal.operator.mutation.impl.PolynomialMutation;
import org.uma.jmetal.problem.doubleproblem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.zdt.ZDT1;
import org.uma.jmetal.problem.multiobjective.zdt.ZDT2;
import org.uma.jmetal.problem.multiobjective.zdt.ZDT3;
import org.uma.jmetal.problem.multiobjective.zdt.ZDT4;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;

/**
 * jMetal 7.4 NSGA-II benchmark runner — ZDT1-ZDT4.
 *
 * Parámetros (deben coincidir con nsga2_runner.cpp):
 *   Population : 100
 *   Max evals  : 25000
 *   Crossover  : SBX,  pc=0.9,  eta_c=20
 *   Mutation   : PM,   pm=1/D,  eta_m=20
 *   ZDT4       : dim=10  (resto: dim=30)
 *
 * Salida CSV: system,problem,run,time_ms,igd
 *
 * Uso:
 *   javac -cp <jmetal-jars> JMetalNSGA2Benchmark.java
 *   java  -server -Xms512m -Xmx2g -cp <jmetal-jars>:. JMetalNSGA2Benchmark results_jmetal.csv [NUM_RUNS [BASE_SEED]]
 */
public class JMetalNSGA2Benchmark {

    static final int    POP_SIZE   = 100;
    static final int    MAX_EVALS  = 25000;
    static final double CX_PROB    = 0.9;
    static final double CX_ETA     = 20.0;
    static final double MUT_ETA    = 20.0;
    static final int    WARMUP     = 3;   // runs para JIT (no se graban)
    static final int    REF_PTS    = 500;
    static final int    NUM_RUNS   = 30;
    static final long   BASE_SEED  = 42L;

    record ProblemSpec(String name, DoubleProblem problem, int dim) {}

    public static void main(String[] args) throws Exception {
        String csvPath  = (args.length > 0) ? args[0] : "results_jmetal.csv";
        int    numRuns  = (args.length > 1) ? Integer.parseInt(args[1])  : NUM_RUNS;
        long   baseSeed = (args.length > 2) ? Long.parseLong(args[2])    : BASE_SEED;

        System.err.println("=== jMetal 7.4 NSGA-II MO Benchmark ===");
        System.err.printf("pop=%d  max_evals=%d  runs=%d  seed=%d%n",
                          POP_SIZE, MAX_EVALS, numRuns, baseSeed);

        var specs = List.of(
            new ProblemSpec("ZDT1", new ZDT1(30), 30),
            new ProblemSpec("ZDT2", new ZDT2(30), 30),
            new ProblemSpec("ZDT3", new ZDT3(30), 30),
            new ProblemSpec("ZDT4", new ZDT4(10), 10)
        );

        try (var pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println("system,problem,run,time_ms,igd");

            for (var spec : specs) {
                System.err.println("  " + spec.name() + " (dim=" + spec.dim() + ")...");

                // Warmup — permite al JIT compilar los hot paths
                System.err.println("    warmup (" + WARMUP + " runs)...");
                for (int w = 0; w < WARMUP; w++) {
                    org.uma.jmetal.util.pseudorandom.JMetalRandom.getInstance()
                        .setSeed(baseSeed + w);
                    runOnce(spec.problem());
                }

                // Runs medidos
                for (int r = 0; r < numRuns; r++) {
                    org.uma.jmetal.util.pseudorandom.JMetalRandom.getInstance()
                        .setSeed(baseSeed + WARMUP + r);
                    long t0     = System.nanoTime();
                    var  result = runOnce(spec.problem());
                    long t1     = System.nanoTime();

                    double timeMs = (t1 - t0) / 1_000_000.0;
                    double igd    = computeIGD(spec.name(), result);

                    pw.printf("jmetal,%s,%d,%.3f,%.6f%n",
                              spec.name(), r + 1, timeMs, igd);
                    pw.flush();
                }
            }
        }
        System.err.println("Done. Resultados en: " + csvPath);
    }

    static List<DoubleSolution> runOnce(DoubleProblem problem) {
        double pm = 1.0 / problem.numberOfVariables();
        var algorithm = new NSGAIIBuilder<>(
                problem,
                new SBXCrossover(CX_PROB, CX_ETA),
                new PolynomialMutation(pm, MUT_ETA),
                POP_SIZE)
            .setMaxEvaluations(MAX_EVALS)
            .setSolutionListEvaluator(new SequentialSolutionListEvaluator<>())
            .build();
        algorithm.run();
        return algorithm.result();
    }

    // IGD analítico con REF_PTS puntos (mismo criterio que nsga2_runner.cpp)
    static double computeIGD(String name, List<DoubleSolution> front) {
        double[][] ref = makeRefFront(name);
        double sum = 0.0;
        for (double[] r : ref) {
            double minDist = Double.MAX_VALUE;
            for (var sol : front) {
                double d = 0.0;
                for (int i = 0; i < r.length; i++) {
                    double diff = r[i] - sol.objectives()[i];
                    d += diff * diff;
                }
                minDist = Math.min(minDist, Math.sqrt(d));
            }
            sum += minDist;
        }
        return sum / ref.length;
    }

    static double[][] makeRefFront(String name) {
        double[][] pf = new double[REF_PTS][2];
        if (name.equals("ZDT1") || name.equals("ZDT4")) {
            for (int i = 0; i < REF_PTS; i++) {
                double f1 = (double) i / (REF_PTS - 1);
                pf[i][0] = f1;
                pf[i][1] = 1.0 - Math.sqrt(f1);
            }
        } else if (name.equals("ZDT2")) {
            for (int i = 0; i < REF_PTS; i++) {
                double f1 = (double) i / (REF_PTS - 1);
                pf[i][0] = f1;
                pf[i][1] = 1.0 - f1 * f1;
            }
        } else if (name.equals("ZDT3")) {
            double[][] segs = {
                {0.0,         0.0830015349},
                {0.182228780, 0.257762881},
                {0.409525530, 0.453861524},
                {0.617876091, 0.652585810},
                {0.823295506, 0.851832239}
            };
            int pps = REF_PTS / segs.length;
            int idx = 0;
            for (double[] seg : segs) {
                double lo = seg[0], hi = seg[1];
                for (int j = 0; j < pps && idx < REF_PTS; j++) {
                    double f1 = lo + (hi - lo) * j / (pps - 1);
                    pf[idx][0] = f1;
                    pf[idx][1] = 1.0 - Math.sqrt(f1) - f1 * Math.sin(10 * Math.PI * f1);
                    idx++;
                }
            }
        }
        return pf;
    }
}
