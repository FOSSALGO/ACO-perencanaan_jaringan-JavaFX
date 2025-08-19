import java.util.*;

public class AntColonyOptimization {

    private int numAnts;       // jumlah semut
    private int numNodes;      // jumlah node
    private int maxIterations; // jumlah iterasi maksimum
    private double alpha;      // pengaruh feromon
    private double beta;       // pengaruh heuristik (1/jarak)
    private double evaporation; // laju penguapan feromon
    private double Q;           // konstanta penambah feromon

    private double[][] graph;   // matriks adjacency
    private double[][] pheromone; // matriks feromon
    private Random random;

    private int source;
    private int destination;

    public AntColonyOptimization(double[][] graph, int source, int destination,
                                 int numAnts, int maxIterations,
                                 double alpha, double beta,
                                 double evaporation, double Q) {
        this.graph = graph;
        this.numNodes = graph.length;
        this.source = source;
        this.destination = destination;
        this.numAnts = numAnts;
        this.maxIterations = maxIterations;
        this.alpha = alpha;
        this.beta = beta;
        this.evaporation = evaporation;
        this.Q = Q;
        this.pheromone = new double[numNodes][numNodes];
        this.random = new Random();

        // Inisialisasi feromon
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                pheromone[i][j] = 0.1; // nilai awal kecil
            }
        }
    }

    public Result solve() {
        List<Integer> bestPath = null;
        double bestCost = Double.MAX_VALUE;

        for (int iter = 0; iter < maxIterations; iter++) {
            List<List<Integer>> allPaths = new ArrayList<>();
            List<Double> allCosts = new ArrayList<>();

            // Bangun solusi untuk setiap semut
            for (int k = 0; k < numAnts; k++) {
                List<Integer> path = buildPath();
                double cost = calculateCost(path);
                allPaths.add(path);
                allCosts.add(cost);

                if (cost < bestCost) {
                    bestCost = cost;
                    bestPath = new ArrayList<>(path);
                }
            }

            // Update feromon
            evaporatePheromone();
            for (int k = 0; k < numAnts; k++) {
                depositPheromone(allPaths.get(k), allCosts.get(k));
            }

            System.out.println("Iterasi " + (iter+1) + " → Jalur: " + allPaths.get(0) + " | Cost: " + allCosts.get(0));
        }

        return new Result(bestPath, bestCost);
    }

    private List<Integer> buildPath() {
        List<Integer> path = new ArrayList<>();
        boolean[] visited = new boolean[numNodes];
        int currentNode = source;
        path.add(currentNode);
        visited[currentNode] = true;

        while (currentNode != destination) {
            int nextNode = selectNextNode(currentNode, visited);
            if (nextNode == -1) break; // tidak ada jalur
            path.add(nextNode);
            visited[nextNode] = true;
            currentNode = nextNode;
        }

        return path;
    }

    private int selectNextNode(int currentNode, boolean[] visited) {
        double[] probabilities = new double[numNodes];
        double sum = 0.0;

        for (int j = 0; j < numNodes; j++) {
            if (!visited[j] && graph[currentNode][j] > 0) {
                probabilities[j] = Math.pow(pheromone[currentNode][j], alpha)
                        * Math.pow(1.0 / graph[currentNode][j], beta);
                sum += probabilities[j];
            }
        }

        if (sum == 0) return -1; // tidak ada jalur

        // Roulette wheel selection
        double r = random.nextDouble() * sum;
        double total = 0;
        for (int j = 0; j < numNodes; j++) {
            total += probabilities[j];
            if (total >= r) {
                return j;
            }
        }

        return -1;
    }

    private double calculateCost(List<Integer> path) {
        double cost = 0.0;
        for (int i = 0; i < path.size() - 1; i++) {
            int u = path.get(i);
            int v = path.get(i + 1);
            cost += graph[u][v];
        }
        return cost;
    }

    private void evaporatePheromone() {
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                pheromone[i][j] *= (1 - evaporation);
            }
        }
    }

    private void depositPheromone(List<Integer> path, double cost) {
        for (int i = 0; i < path.size() - 1; i++) {
            int u = path.get(i);
            int v = path.get(i + 1);
            pheromone[u][v] += Q / cost;
            pheromone[v][u] += Q / cost; // jika graf tak berarah
        }
    }

    // Class hasil
    public static class Result {
        public List<Integer> bestPath;
        public double bestCost;

        public Result(List<Integer> path, double cost) {
            this.bestPath = path;
            this.bestCost = cost;
        }
    }

    // --- Main Program ---
    public static void main(String[] args) {
        // Matriks adjacency (contoh: 5 node, nilai = latensi)
        double[][] graph = {
                {0, 10, 0, 30, 100},
                {10, 0, 50, 0, 0},
                {0, 50, 0, 20, 10},
                {30, 0, 20, 0, 60},
                {100, 0, 10, 60, 0}
        };

        AntColonyOptimization aco = new AntColonyOptimization(
                graph,
                0, // source (node 0)
                4, // destination (node 4)
                10, // jumlah semut
                50, // jumlah iterasi
                1.0, // alpha
                5.0, // beta
                0.5, // evaporation
                100  // Q
        );

        Result result = aco.solve();
        System.out.println("\n=== Hasil Akhir ===");
        System.out.println("Jalur terbaik: " + result.bestPath);
        System.out.println("Biaya (latensi total): " + result.bestCost);
    }
}
