import org.uma.jmetal.algorithm.multiobjective.nsgaiii.NSGAIIIBuilder;
import org.uma.jmetal.operator.crossover.impl.SBXCrossover;
import org.uma.jmetal.operator.mutation.impl.PolynomialMutation;
import org.uma.jmetal.operator.selection.impl.BinaryTournamentSelection;
import org.uma.jmetal.problem.doubleproblem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ1;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ2;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ3;
import org.uma.jmetal.problem.multiobjective.dtlz.DTLZ4;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;

/**
 * jMetal 7.4 NSGA-III benchmark runner — DTLZ1-4 (3 objectives).
 *
 * Parámetros (deben coincidir con nsga3_runner.cpp):
 *   n_obj=3, divisions=12  → 91 ref points → pop_size=92
 *   maxIterations=272       → 271×92+92 = 25024 ≈ 25000 evaluaciones totales
 *   SBX(p=0.9, η_c=20), PolMut(p=1/dim, η_m=20)
 *   DTLZ1: 7 variables, DTLZ2/3/4: 12 variables
 *
 * Salida CSV: system,problem,run,time_ms,igd
 *
 * Uso:
 *   javac -cp <jmetal-jars> JMetalNSGA3Benchmark.java
 *   java -server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch \
 *        -XX:+DisableExplicitGC -Djava.awt.headless=true \
 *        -cp <jmetal-jars>:. JMetalNSGA3Benchmark results_jmetal.csv [NUM_RUNS [BASE_SEED]]
 */
public class JMetalNSGA3Benchmark {

    static final int    N_OBJ         = 3;
    static final int    DIVISIONS     = 12;   // 91 ref points → pop_size=92
    static final int    MAX_ITERATIONS = 272;  // 271 EA loops × 92 = 24932 + 92 init ≈ 25024 evals
    static final double SBX_PROB      = 0.9;
    static final double SBX_ETA       = 20.0;
    static final double MUT_ETA       = 20.0;
    static final int    WARMUP        = 3;
    static final int    REF_DIV       = 30;   // Das-Dennis divisions for reference front
    static final int    NUM_RUNS      = 30;
    static final long   BASE_SEED     = 42L;

    record ProblemSpec(String name, DoubleProblem problem, int dim) {}

    public static void main(String[] args) throws Exception {
        String csvPath   = (args.length > 0) ? args[0]                    : "results_jmetal.csv";
        int    numRuns   = (args.length > 1) ? Integer.parseInt(args[1])  : NUM_RUNS;
        long   baseSeed  = (args.length > 2) ? Long.parseLong(args[2])    : BASE_SEED;

        System.err.println("=== jMetal 7.4 NSGA-III Benchmark ===");
        System.err.printf(
            "n_obj=%d  divisions=%d  max_iterations=%d  SBX(p=%.1f,η=%.0f)  PolMut(η=%.0f)  runs=%d  seed=%d%n",
            N_OBJ, DIVISIONS, MAX_ITERATIONS, SBX_PROB, SBX_ETA, MUT_ETA, numRuns, baseSeed);

        var specs = List.of(
            new ProblemSpec("DTLZ1", new DTLZ1(7,  N_OBJ), 7),
            new ProblemSpec("DTLZ2", new DTLZ2(12, N_OBJ), 12),
            new ProblemSpec("DTLZ3", new DTLZ3(12, N_OBJ), 12),
            new ProblemSpec("DTLZ4", new DTLZ4(12, N_OBJ), 12)
        );

        double[][] ref1   = makeDTLZ1Ref();
        double[][] refSph = makeSphereRef();

        try (var pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println("system,problem,run,time_ms,igd");

            for (var spec : specs) {
                System.err.println("  " + spec.name() + " (dim=" + spec.dim() + ")...");

                double[][] ref = spec.name().equals("DTLZ1") ? ref1 : refSph;

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
                    double igd    = computeIGD(result, ref);

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
        var crossover = new SBXCrossover(SBX_PROB, SBX_ETA);
        var mutation  = new PolynomialMutation(pm, MUT_ETA);
        var selection = new BinaryTournamentSelection<DoubleSolution>(
                new RankingAndCrowdingDistanceComparator<>());

        var algorithm = new NSGAIIIBuilder<>(problem)
                .setCrossoverOperator(crossover)
                .setMutationOperator(mutation)
                .setSelectionOperator(selection)
                .setMaxIterations(MAX_ITERATIONS)
                .setNumberOfDivisions(DIVISIONS)
                .build();

        algorithm.run();
        return algorithm.result(); // non-dominated front
    }

    // IGD (Euclidean, p=1): mean distance from each ref point to closest solution
    static double computeIGD(List<DoubleSolution> front, double[][] ref) {
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

    // DTLZ1 reference front (3 obj): simplex sum=0.5
    // Das-Dennis with REF_DIV divisions → C(REF_DIV+2, 2) points
    static double[][] makeDTLZ1Ref() {
        var pts = generateDasDennis(N_OBJ, REF_DIV);
        double[][] ref = new double[pts.size()][N_OBJ];
        for (int i = 0; i < pts.size(); i++)
            for (int j = 0; j < N_OBJ; j++)
                ref[i][j] = pts.get(i)[j] * 0.5;
        return ref;
    }

    // DTLZ2/3/4 reference front (3 obj): unit quarter sphere
    // Das-Dennis points normalized to unit sphere
    static double[][] makeSphereRef() {
        var pts = generateDasDennis(N_OBJ, REF_DIV);
        double[][] ref = new double[pts.size()][N_OBJ];
        for (int i = 0; i < pts.size(); i++) {
            double norm = 0.0;
            for (int j = 0; j < N_OBJ; j++)
                norm += pts.get(i)[j] * pts.get(i)[j];
            norm = Math.sqrt(norm);
            if (norm < 1e-14) norm = 1.0;
            for (int j = 0; j < N_OBJ; j++)
                ref[i][j] = pts.get(i)[j] / norm;
        }
        return ref;
    }

    // Das-Dennis reference point generation (recursive)
    static java.util.List<double[]> generateDasDennis(int nObj, int divisions) {
        var result = new java.util.ArrayList<double[]>();
        double[] point = new double[nObj];
        generateRecursive(result, point, nObj, divisions, divisions, 0);
        return result;
    }

    static void generateRecursive(java.util.List<double[]> result, double[] point,
                                   int nObj, int left, int total, int element) {
        if (element == nObj - 1) {
            point[element] = (double) left / total;
            result.add(point.clone());
        } else {
            for (int i = 0; i <= left; i++) {
                point[element] = (double) i / total;
                generateRecursive(result, point, nObj, left - i, total, element + 1);
            }
        }
    }
}
