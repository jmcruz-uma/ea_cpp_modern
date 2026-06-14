import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.nsgaiii.NSGAIIIBuilder;
import org.uma.jmetal.operator.crossover.impl.SBXCrossover;
import org.uma.jmetal.operator.mutation.impl.PolynomialMutation;
import org.uma.jmetal.operator.selection.impl.BinaryTournamentSelection;
import org.uma.jmetal.problem.doubleproblem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F1;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F2;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F3;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F4;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F5;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F6;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F7;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F8;
import org.uma.jmetal.problem.multiobjective.lz09.LZ09F9;
import org.uma.jmetal.solution.doublesolution.DoubleSolution;
import org.uma.jmetal.util.comparator.RankingAndCrowdingDistanceComparator;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * jMetal 7.4 benchmark runner — LZ09F1-F9 (Li & Zhang 2009).
 *
 * Parámetros (deben coincidir con lz09_runner.cpp):
 *   NSGA-II  (F1-F5, F7-F9): pop=100, max_evals=25000, SBX(pc=0.9, η=20), PM(η=20, p=1/dim)
 *   NSGA-III (F6, 3-obj)   : divisions=12 → 91 ref pts → pop=92, max_iterations=272
 *
 * Salida CSV: system,problem,run,time_ms,igd
 *
 * Uso:
 *   javac -cp <jmetal-jars> JMetalLZ09Benchmark.java
 *   java  -server -Xms512m -Xmx2g -XX:+UseG1GC -XX:+AlwaysPreTouch \
 *         -XX:+DisableExplicitGC -Djava.awt.headless=true \
 *         -cp <jmetal-jars>:. JMetalLZ09Benchmark results_jmetal_lz09.csv [NUM_RUNS [BASE_SEED]]
 */
public class JMetalLZ09Benchmark {

    // NSGA-II (F1-F5, F7-F9)
    static final int    POP_SIZE     = 100;
    static final int    MAX_EVALS    = 25000;
    // NSGA-III (F6)
    static final int    DIVISIONS    = 12;   // 91 ref pts → pop_size=92
    static final int    MAX_ITER_F6  = 272;  // 272 × 92 ≈ 25024 evaluaciones
    // Operadores (ambos algoritmos)
    static final double CX_PROB      = 0.9;
    static final double CX_ETA       = 20.0;
    static final double MUT_ETA      = 20.0;
    // Benchmark
    static final int    WARMUP       = 3;
    static final int    NUM_RUNS     = 30;
    static final long   BASE_SEED    = 42L;
    static final int    REF_PTS      = 500;
    static final int    REF_DIV      = 30;   // Das-Dennis divisions para frente 3-obj

    enum Algo { NSGA2, NSGA3 }

    record ProblemSpec(String name, DoubleProblem problem, Algo algo, double[][] refFront) {}

    public static void main(String[] args) throws Exception {
        String csvPath  = (args.length > 0) ? args[0]                   : "results_jmetal_lz09.csv";
        int    numRuns  = (args.length > 1) ? Integer.parseInt(args[1]) : NUM_RUNS;
        long   baseSeed = (args.length > 2) ? Long.parseLong(args[2])   : BASE_SEED;

        System.err.println("=== jMetal 7.4 LZ09 Benchmark ===");
        System.err.printf("NSGA-II:  pop=%d  max_evals=%d  runs=%d  seed=%d%n",
                          POP_SIZE, MAX_EVALS, numRuns, baseSeed);
        System.err.printf("NSGA-III: divisions=%d  max_iter=%d  (LZ09F6)%n",
                          DIVISIONS, MAX_ITER_F6);

        double[][] ref21  = makeRef21();    // ptype=21: f2=1-sqrt(f1)
        double[][] ref22  = makeRef22();    // ptype=22: f2=1-f1²
        double[][] refSph = makeSphereRef(); // ptype=31: cuarto de esfera

        var specs = List.of(
            new ProblemSpec("LZ09F1", new LZ09F1(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F2", new LZ09F2(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F3", new LZ09F3(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F4", new LZ09F4(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F5", new LZ09F5(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F6", new LZ09F6(), Algo.NSGA3, refSph),
            new ProblemSpec("LZ09F7", new LZ09F7(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F8", new LZ09F8(), Algo.NSGA2, ref21),
            new ProblemSpec("LZ09F9", new LZ09F9(), Algo.NSGA2, ref22)
        );

        try (var pw = new PrintWriter(new FileWriter(csvPath))) {
            pw.println("system,problem,run,time_ms,igd");

            for (var spec : specs) {
                System.err.printf("  %s (dim=%d, %s)...%n",
                                  spec.name(),
                                  spec.problem().numberOfVariables(),
                                  spec.algo() == Algo.NSGA3 ? "3-obj, NSGA-III" : "2-obj, NSGA-II");

                System.err.println("    warmup (" + WARMUP + " runs)...");
                for (int w = 0; w < WARMUP; w++) {
                    org.uma.jmetal.util.pseudorandom.JMetalRandom.getInstance()
                        .setSeed(baseSeed + w);
                    runOnce(spec);
                }

                for (int r = 0; r < numRuns; r++) {
                    org.uma.jmetal.util.pseudorandom.JMetalRandom.getInstance()
                        .setSeed(baseSeed + WARMUP + r);

                    long t0     = System.nanoTime();
                    var  result = runOnce(spec);
                    long t1     = System.nanoTime();

                    double timeMs = (t1 - t0) / 1_000_000.0;
                    double igd    = computeIGD(result, spec.refFront());

                    pw.printf("jmetal,%s,%d,%.3f,%.6f%n",
                              spec.name(), r + 1, timeMs, igd);
                    pw.flush();
                }
            }
        }
        System.err.println("Done. Resultados en: " + csvPath);
    }

    static List<DoubleSolution> runOnce(ProblemSpec spec) {
        return switch (spec.algo()) {
            case NSGA2 -> runOnceNSGA2(spec.problem());
            case NSGA3 -> runOnceNSGA3(spec.problem());
        };
    }

    static List<DoubleSolution> runOnceNSGA2(DoubleProblem problem) {
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

    static List<DoubleSolution> runOnceNSGA3(DoubleProblem problem) {
        double pm = 1.0 / problem.numberOfVariables();
        var crossover = new SBXCrossover(CX_PROB, CX_ETA);
        var mutation  = new PolynomialMutation(pm, MUT_ETA);
        var selection = new BinaryTournamentSelection<DoubleSolution>(
                new RankingAndCrowdingDistanceComparator<>());
        var algorithm = new NSGAIIIBuilder<>(problem)
                .setCrossoverOperator(crossover)
                .setMutationOperator(mutation)
                .setSelectionOperator(selection)
                .setMaxIterations(MAX_ITER_F6)
                .setNumberOfDivisions(DIVISIONS)
                .build();
        algorithm.run();
        return algorithm.result();
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

    // ── Frentes de referencia ─────────────────────────────────────────────────

    // ptype=21: f2 = 1 - sqrt(f1), f1 in [0,1]   (LZ09F1-F5, F7-F8)
    static double[][] makeRef21() {
        double[][] pf = new double[REF_PTS][2];
        for (int i = 0; i < REF_PTS; i++) {
            double f1 = (double) i / (REF_PTS - 1);
            pf[i][0] = f1;
            pf[i][1] = 1.0 - Math.sqrt(f1);
        }
        return pf;
    }

    // ptype=22: f2 = 1 - f1², f1 in [0,1]   (LZ09F9)
    static double[][] makeRef22() {
        double[][] pf = new double[REF_PTS][2];
        for (int i = 0; i < REF_PTS; i++) {
            double f1 = (double) i / (REF_PTS - 1);
            pf[i][0] = f1;
            pf[i][1] = 1.0 - f1 * f1;
        }
        return pf;
    }

    // ptype=31: cuarto de esfera f1²+f2²+f3²=1, fi≥0   (LZ09F6)
    // Das-Dennis con REF_DIV divisiones, normalizado a esfera unitaria.
    static double[][] makeSphereRef() {
        var pts = generateDasDennis(3, REF_DIV);
        double[][] ref = new double[pts.size()][3];
        for (int i = 0; i < pts.size(); i++) {
            double norm = 0.0;
            for (int j = 0; j < 3; j++) norm += pts.get(i)[j] * pts.get(i)[j];
            norm = Math.sqrt(norm);
            if (norm < 1e-14) norm = 1.0;
            for (int j = 0; j < 3; j++) ref[i][j] = pts.get(i)[j] / norm;
        }
        return ref;
    }

    // Das-Dennis reference point generation (recursive simplex lattice)
    static List<double[]> generateDasDennis(int nObj, int divisions) {
        var result = new ArrayList<double[]>();
        double[] point = new double[nObj];
        generateRecursive(result, point, nObj, divisions, divisions, 0);
        return result;
    }

    static void generateRecursive(List<double[]> result, double[] point,
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
