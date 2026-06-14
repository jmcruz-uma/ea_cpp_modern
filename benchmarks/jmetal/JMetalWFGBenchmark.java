import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.operator.crossover.impl.SBXCrossover;
import org.uma.jmetal.operator.mutation.impl.PolynomialMutation;
import org.uma.jmetal.problem.doubleproblem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.wfg.WFG1;
import org.uma.jmetal.problem.multiobjective.wfg.WFG2;
import org.uma.jmetal.problem.multiobjective.wfg.WFG3;
import org.uma.jmetal.problem.multiobjective.wfg.WFG4;
import org.uma.jmetal.problem.multiobjective.wfg.WFG5;
import org.uma.jmetal.problem.multiobjective.wfg.WFG6;
import org.uma.jmetal.problem.multiobjective.wfg.WFG7;
import org.uma.jmetal.problem.multiobjective.wfg.WFG8;
import org.uma.jmetal.problem.multiobjective.wfg.WFG9;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * jMetal 7.4 NSGA-II benchmark runner — WFG1-9.
 *
 * Parámetros (deben coincidir con wfg_runner.cpp):
 *   Population : 100
 *   Max evals  : 25000
 *   Crossover  : SBX, pc=0.9, eta_c=20
 *   Mutation   : PM,  pm=1/dim, eta_m=20
 *   k=2, l=4, M=2, dim=6  (DefaultWFGSettings)
 *
 * Salida CSV: system,problem,run,time_ms,igd
 *
 * Uso:
 *   javac -cp <jmetal-jars> JMetalWFGBenchmark.java
 *   java  -server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch \
 *         -XX:+DisableExplicitGC -Djava.awt.headless=true \
 *         -cp <jmetal-jars>:. JMetalWFGBenchmark results_jmetal_wfg.csv [NUM_RUNS [BASE_SEED [REF_DIR]]]
 */
public class JMetalWFGBenchmark {

    static final int    POP_SIZE  = 100;
    static final int    MAX_EVALS = 25000;
    static final double CX_PROB   = 0.9;
    static final double CX_ETA    = 20.0;
    static final double MUT_ETA   = 20.0;
    static final int    WARMUP    = 3;
    static final int    NUM_RUNS  = 30;
    static final long   BASE_SEED = 42L;
    static final String DEFAULT_REF_DIR =
        "/home/alumno/jMetal/resources/referenceFrontsCSV";

    record ProblemSpec(String name, DoubleProblem problem) {}

    public static void main(String[] args) throws Exception {
        String csvPath  = (args.length > 0) ? args[0] : "results_jmetal_wfg.csv";
        int    numRuns  = (args.length > 1) ? Integer.parseInt(args[1]) : NUM_RUNS;
        long   baseSeed = (args.length > 2) ? Long.parseLong(args[2])   : BASE_SEED;
        String refDir   = (args.length > 3) ? args[3]                   : DEFAULT_REF_DIR;

        System.err.println("=== jMetal 7.4 NSGA-II WFG Benchmark ===");
        System.err.printf("pop=%d  max_evals=%d  runs=%d  seed=%d%n",
                          POP_SIZE, MAX_EVALS, numRuns, baseSeed);
        System.err.println("ref_dir=" + refDir);

        // All problems use DefaultWFGSettings: k=2, l=4, M=2
        var specs = List.of(
            new ProblemSpec("WFG1", new WFG1()),
            new ProblemSpec("WFG2", new WFG2()),
            new ProblemSpec("WFG3", new WFG3()),
            new ProblemSpec("WFG4", new WFG4()),
            new ProblemSpec("WFG5", new WFG5()),
            new ProblemSpec("WFG6", new WFG6()),
            new ProblemSpec("WFG7", new WFG7()),
            new ProblemSpec("WFG8", new WFG8()),
            new ProblemSpec("WFG9", new WFG9())
        );

        try (var pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println("system,problem,run,time_ms,igd");

            for (var spec : specs) {
                System.err.println("  " + spec.name() + " (dim=" +
                                   spec.problem().numberOfVariables() + ")...");

                double[][] ref = loadRefFront(refDir + "/" + spec.name() + ".2D.csv");
                System.err.println("    ref front: " + ref.length + " pts");

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
                    double igd    = computeIGD(ref, result);

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

    static double[][] loadRefFront(String path) throws Exception {
        var lines = new ArrayList<double[]>();
        try (var br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                String[] parts = line.split(",");
                lines.add(new double[]{Double.parseDouble(parts[0]),
                                       Double.parseDouble(parts[1])});
            }
        }
        return lines.toArray(new double[0][]);
    }

    static double computeIGD(double[][] ref, List<DoubleSolution> front) {
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
}
