import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.operator.crossover.impl.SBXCrossover;
import org.uma.jmetal.operator.mutation.impl.PolynomialMutation;
import org.uma.jmetal.problem.doubleproblem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.uf.UF1;
import org.uma.jmetal.problem.multiobjective.uf.UF2;
import org.uma.jmetal.problem.multiobjective.uf.UF3;
import org.uma.jmetal.problem.multiobjective.uf.UF4;
import org.uma.jmetal.problem.multiobjective.uf.UF5;
import org.uma.jmetal.problem.multiobjective.uf.UF6;
import org.uma.jmetal.problem.multiobjective.uf.UF7;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;

/**
 * jMetal 7.4 NSGA-II benchmark runner — UF1-UF7 (CEC2009).
 *
 * Parámetros (deben coincidir con uf_runner.cpp):
 *   Population : 100
 *   Max evals  : 25000
 *   Crossover  : SBX, pc=0.9, eta_c=20
 *   Mutation   : PM,  pm=1/dim, eta_m=20
 *   Variables  : 30 (todos UF1-7)
 *
 * Salida CSV: system,problem,run,time_ms,igd
 *
 * Uso:
 *   javac -cp <jmetal-jars> JMetalUFBenchmark.java
 *   java  -server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch \
 *         -XX:+DisableExplicitGC -Djava.awt.headless=true \
 *         -cp <jmetal-jars>:. JMetalUFBenchmark results_jmetal_uf.csv [NUM_RUNS [BASE_SEED]]
 */
public class JMetalUFBenchmark {

    static final int    POP_SIZE  = 100;
    static final int    MAX_EVALS = 25000;
    static final double CX_PROB   = 0.9;
    static final double CX_ETA    = 20.0;
    static final double MUT_ETA   = 20.0;
    static final int    WARMUP    = 3;
    static final int    NUM_RUNS  = 30;
    static final long   BASE_SEED = 42L;
    static final int    REF_PTS   = 500;

    record ProblemSpec(String name, DoubleProblem problem) {}

    public static void main(String[] args) throws Exception {
        String csvPath  = (args.length > 0) ? args[0] : "results_jmetal_uf.csv";
        int    numRuns  = (args.length > 1) ? Integer.parseInt(args[1]) : NUM_RUNS;
        long   baseSeed = (args.length > 2) ? Long.parseLong(args[2])   : BASE_SEED;

        System.err.println("=== jMetal 7.4 NSGA-II UF Benchmark ===");
        System.err.printf("pop=%d  max_evals=%d  runs=%d  seed=%d%n",
                          POP_SIZE, MAX_EVALS, numRuns, baseSeed);

        var specs = List.of(
            new ProblemSpec("UF1", new UF1(30)),
            new ProblemSpec("UF2", new UF2(30)),
            new ProblemSpec("UF3", new UF3(30)),
            new ProblemSpec("UF4", new UF4(30)),
            new ProblemSpec("UF5", new UF5(30, 10, 0.1)),
            new ProblemSpec("UF6", new UF6(30, 2, 0.1)),
            new ProblemSpec("UF7", new UF7(30))
        );

        try (var pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println("system,problem,run,time_ms,igd");

            for (var spec : specs) {
                System.err.println("  " + spec.name() + " (dim=30)...");

                System.err.println("    warmup (" + WARMUP + " runs)...");
                for (int w = 0; w < WARMUP; w++) {
                    org.uma.jmetal.util.pseudorandom.JMetalRandom.getInstance()
                        .setSeed(baseSeed + w);
                    runOnce(spec.problem());
                }

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

    // IGD con mismo criterio que uf_runner.cpp
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
        return switch (name) {
            case "UF1", "UF2", "UF3" -> {
                // f2 = 1 - sqrt(f1), f1 in [0,1]
                double[][] pf = new double[REF_PTS][2];
                for (int i = 0; i < REF_PTS; i++) {
                    double f1 = (double) i / (REF_PTS - 1);
                    pf[i][0] = f1;
                    pf[i][1] = 1.0 - Math.sqrt(f1);
                }
                yield pf;
            }
            case "UF4" -> {
                // f2 = 1 - f1^2, f1 in [0,1]
                double[][] pf = new double[REF_PTS][2];
                for (int i = 0; i < REF_PTS; i++) {
                    double f1 = (double) i / (REF_PTS - 1);
                    pf[i][0] = f1;
                    pf[i][1] = 1.0 - f1 * f1;
                }
                yield pf;
            }
            case "UF5" -> {
                // 21 puntos discretos en f1+f2=1 con f1 = k/20
                double[][] pf = new double[21][2];
                for (int k = 0; k <= 20; k++) {
                    pf[k][0] = k / 20.0;
                    pf[k][1] = 1.0 - k / 20.0;
                }
                yield pf;
            }
            case "UF6" -> {
                // Dos segmentos de f1+f2=1: f1 in [1/4,1/2] y [3/4,1], 250 pts c/u
                int half = REF_PTS / 2;
                double[][] pf = new double[REF_PTS][2];
                for (int i = 0; i < half; i++) {
                    double f1 = 0.25 + 0.25 * i / (half - 1);
                    pf[i][0] = f1;
                    pf[i][1] = 1.0 - f1;
                }
                for (int i = 0; i < half; i++) {
                    double f1 = 0.75 + 0.25 * i / (half - 1);
                    pf[half + i][0] = f1;
                    pf[half + i][1] = 1.0 - f1;
                }
                yield pf;
            }
            case "UF7" -> {
                // f1+f2=1, f1 in [0,1]
                double[][] pf = new double[REF_PTS][2];
                for (int i = 0; i < REF_PTS; i++) {
                    double f1 = (double) i / (REF_PTS - 1);
                    pf[i][0] = f1;
                    pf[i][1] = 1.0 - f1;
                }
                yield pf;
            }
            default -> new double[0][0];
        };
    }
}
