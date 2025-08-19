import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Point2D;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.*;
import javafx.scene.input.MouseButton;
import javafx.scene.layout.*;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.stage.Stage;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Ant Colony Optimization (ACO) Visualizer for Internet Routing
 * Fitur:
 *  - Canvas: zoom (scroll), pan (drag kanan), drag vertex (drag kiri), tambah vertex (klik kiri),
 *    buat/hapus edge (Shift + klik dua vertex).
 *  - Panel parameter ACO + tombol Run.
 *  - Highlight jalur terbaik.
 *
 * Kompilasi (Java 11+):
 *   javac ACOJavaFXApp.java
 * Jalankan:
 *   java ACOJavaFXApp
 */
public class ACOJavaFXApp extends Application {

    // ====== Model graf sederhana ======
    static class Node {
        int id;
        double x, y;
        Node(int id, double x, double y) { this.id = id; this.x = x; this.y = y; }
    }

    static class Edge {
        int a, b; // node ids
        double weight; // Euclidean distance (auto)
        Edge(int a, int b, double weight) { this.a = a; this.b = b; this.weight = weight; }
        boolean involves(int id){ return a==id || b==id; }
        int other(int id){ return a==id ? b : a; }
    }

    static class Graph {
        List<Node> nodes = new ArrayList<>();
        List<Edge> edges = new ArrayList<>();

        Node addNode(double x, double y){
            Node n = new Node(nodes.size(), x, y);
            nodes.add(n);
            return n;
        }

        void toggleEdge(int a, int b){
            if (a==b) return;
            // ensure a<b for uniqueness
            int u = Math.min(a,b), v = Math.max(a,b);
            for (Iterator<Edge> it = edges.iterator(); it.hasNext();){
                Edge e = it.next();
                int eu = Math.min(e.a,e.b), ev = Math.max(e.a,e.b);
                if (eu==u && ev==v){
                    it.remove(); // remove existing edge
                    return;
                }
            }
            Node na = nodes.get(u), nb = nodes.get(v);
            double w = Math.hypot(na.x - nb.x, na.y - nb.y);
            edges.add(new Edge(u, v, w));
        }

        double[][] toAdjacency(){
            int n = nodes.size();
            double[][] A = new double[n][n];
            for (Edge e: edges){
                A[e.a][e.b] = e.weight;
                A[e.b][e.a] = e.weight; // undirected
            }
            return A;
        }

        void recomputeWeightsEuclidean(){
            for (Edge e: edges){
                Node na = nodes.get(e.a), nb = nodes.get(e.b);
                e.weight = Math.hypot(na.x - nb.x, na.y - nb.y);
            }
        }

        void clear(){
            nodes.clear();
            edges.clear();
        }
    }

    // ====== ACO core ======
    static class ACO {
        private final double[][] graph; // adjacency, 0 = no edge
        private final int n;
        private final int ants, iters;
        private final double alpha, beta, rho, Q;
        private final int source, dest;
        private final Random rng = new Random();

        private final double[][] tau;

        ACO(double[][] graph, int source, int dest, int ants, int iters,
            double alpha, double beta, double rho, double Q) {
            this.graph = graph;
            this.n = graph.length;
            this.source = source;
            this.dest = dest;
            this.ants = ants;
            this.iters = iters;
            this.alpha = alpha;
            this.beta = beta;
            this.rho = rho;
            this.Q = Q;
            this.tau = new double[n][n];
            for (int i=0;i<n;i++) Arrays.fill(tau[i], 0.1);
        }

        Result solve(){
            List<Integer> bestPath = null;
            double bestCost = Double.MAX_VALUE;

            for (int t=0;t<iters;t++){
                List<List<Integer>> paths = new ArrayList<>();
                List<Double> costs = new ArrayList<>();

                for (int k=0;k<ants;k++){
                    List<Integer> path = constructPath();
                    double cost = pathCost(path);
                    paths.add(path);
                    costs.add(cost);
                    if (cost < bestCost){
                        bestCost = cost;
                        bestPath = new ArrayList<>(path);
                    }
                }
                evaporate();
                for (int k=0;k<ants;k++){
                    deposit(paths.get(k), costs.get(k));
                }
            }
            return new Result(bestPath, bestCost);
        }

        private List<Integer> constructPath(){
            boolean[] visited = new boolean[n];
            List<Integer> path = new ArrayList<>();
            int cur = source;
            path.add(cur);
            visited[cur] = true;

            int safety = 0;
            while (cur != dest && safety++ < n*2){
                int next = selectNext(cur, visited);
                if (next == -1) break;
                path.add(next);
                visited[next] = true;
                cur = next;
            }
            return path;
        }

        private int selectNext(int i, boolean[] visited){
            double[] p = new double[n];
            double sum = 0;
            for (int j=0;j<n;j++){
                if (!visited[j] && graph[i][j] > 0){
                    double eta = 1.0 / graph[i][j];
                    p[j] = Math.pow(tau[i][j], alpha) * Math.pow(eta, beta);
                    sum += p[j];
                }
            }
            if (sum == 0) return -1;
            double r = rng.nextDouble() * sum;
            double acc = 0;
            for (int j=0;j<n;j++){
                if (p[j] == 0) continue;
                acc += p[j];
                if (acc >= r) return j;
            }
            return -1;
        }

        private double pathCost(List<Integer> path){
            if (path.size() < 2 || path.get(path.size()-1) != dest) return Double.MAX_VALUE/4;
            double c = 0;
            for (int k=0;k<path.size()-1;k++){
                int u = path.get(k), v = path.get(k+1);
                c += graph[u][v];
            }
            return c;
        }

        private void evaporate(){
            for (int i=0;i<n;i++){
                for (int j=0;j<n;j++){
                    tau[i][j] *= (1 - rho);
                }
            }
        }

        private void deposit(List<Integer> path, double cost){
            if (cost == Double.MAX_VALUE/4) return;
            double add = Q / cost;
            for (int k=0;k<path.size()-1;k++){
                int u = path.get(k), v = path.get(k+1);
                tau[u][v] += add;
                tau[v][u] += add; // undirected
            }
        }

        static class Result {
            final List<Integer> bestPath;
            final double bestCost;
            Result(List<Integer> bestPath, double bestCost){
                this.bestPath = bestPath;
                this.bestCost = bestCost;
            }
        }
    }

    // ====== State UI / Kanvas ======
    private final Graph graph = new Graph();

    private double scale = 1.0;
    private double offsetX = 0;
    private double offsetY = 0;

    private Node draggingNode = null;
    private double lastMouseX, lastMouseY;

    private Integer edgeFirstSelection = null; // for Shift+click edge creation

    private List<Integer> bestPath = Collections.emptyList();
    private double bestCost = Double.NaN;

    // UI controls
    private TextField antsField, itersField, alphaField, betaField, rhoField, qField, srcField, dstField;
    private Label statusLabel;

    private Canvas canvas;

    @Override
    public void start(Stage stage) {
        // Left panel (controls)
        VBox control = buildControlPanel();

        // Canvas
        canvas = new Canvas(1000, 700);
        StackPane canvasHolder = new StackPane(canvas);
        canvasHolder.setStyle("-fx-background-color: #1e1f22;");
        canvas.widthProperty().addListener((obs, o, n) -> draw());
        canvas.heightProperty().addListener((obs, o, n) -> draw());

        initCanvasInteractions();

        BorderPane root = new BorderPane();
        root.setLeft(control);
        root.setCenter(canvasHolder);

        Scene scene = new Scene(root, 1300, 760);
        stage.setTitle("ACO Internet Routing — JavaFX Visualizer");
        stage.setScene(scene);
        stage.show();

        draw();
    }

    private VBox buildControlPanel(){
        antsField = new TextField("20");
        itersField = new TextField("80");
        alphaField = new TextField("1.0");
        betaField = new TextField("5.0");
        rhoField = new TextField("0.5");
        qField = new TextField("100");
        srcField = new TextField("0");
        dstField = new TextField("1");

        GridPane grid = new GridPane();
        grid.setHgap(8);
        grid.setVgap(8);
        int r=0;
        grid.add(new Label("Ants (m)"), 0, r); grid.add(antsField, 1, r++);
        grid.add(new Label("Iterations (T)"), 0, r); grid.add(itersField, 1, r++);
        grid.add(new Label("Alpha (α)"), 0, r); grid.add(alphaField, 1, r++);
        grid.add(new Label("Beta (β)"), 0, r); grid.add(betaField, 1, r++);
        grid.add(new Label("Evaporation (ρ)"), 0, r); grid.add(rhoField, 1, r++);
        grid.add(new Label("Q"), 0, r); grid.add(qField, 1, r++);
        grid.add(new Label("Source"), 0, r); grid.add(srcField, 1, r++);
        grid.add(new Label("Destination"), 0, r); grid.add(dstField, 1, r++);

        Button runBtn = new Button("Run ACO");
        runBtn.setMaxWidth(Double.MAX_VALUE);
        runBtn.setOnAction(e -> runACO());

        Button clearBtn = new Button("Clear");
        clearBtn.setMaxWidth(Double.MAX_VALUE);
        clearBtn.setOnAction(e -> {
            graph.clear();
            bestPath = Collections.emptyList();
            bestCost = Double.NaN;
            draw();
            status("Cleared.");
        });

        Button autoBtn = new Button("Auto Layout");
        autoBtn.setMaxWidth(Double.MAX_VALUE);
        autoBtn.setOnAction(e -> {
            autoLayout();
            draw();
            status("Auto layout applied.");
        });

        statusLabel = new Label("Tip: Klik kiri tambah vertex. Shift+klik dua vertex untuk edge.");
        statusLabel.setWrapText(true);

        VBox left = new VBox(10, titled("ACO Parameters", grid),
                                  runBtn, autoBtn, clearBtn,
                                  titled("Info", statusLabel));
        left.setPadding(new Insets(12));
        left.setPrefWidth(300);
        return left;
    }

    private TitledPane titled(String title, javafx.scene.Node node){
        TitledPane tp = new TitledPane(title, node);
        tp.setExpanded(true);
        return tp;
    }

    private void initCanvasInteractions(){
        GraphicsContext g = canvas.getGraphicsContext2D();

        canvas.setOnMousePressed(e -> {
            lastMouseX = e.getX();
            lastMouseY = e.getY();
            Point2D world = screenToWorld(e.getX(), e.getY());

            if (e.getButton() == MouseButton.SECONDARY){
                // start panning
                return;
            }

            // check if clicked a node
            Node hit = findNodeAt(world.getX(), world.getY());
            if (hit != null){
                if (e.isShiftDown()){
                    // Edge creation: if this is the first selection, store; else toggle edge
                    if (edgeFirstSelection == null){
                        edgeFirstSelection = hit.id;
                        status("Edge: pilih vertex kedua (Shift+klik).");
                    } else {
                        graph.toggleEdge(edgeFirstSelection, hit.id);
                        edgeFirstSelection = null;
                        graph.recomputeWeightsEuclidean();
                        draw();
                        status("Edge diubah.");
                    }
                } else {
                    draggingNode = hit; // start dragging node
                }
            } else {
                if (!e.isShiftDown() && e.getButton()==MouseButton.PRIMARY){
                    Node n = graph.addNode(world.getX(), world.getY());
                    draw();
                    status("Tambah vertex id=" + n.id + " di (" + (int)world.getX() + "," + (int)world.getY() + ")");
                }
            }
        });

        canvas.setOnMouseDragged(e -> {
            Point2D world = screenToWorld(e.getX(), e.getY());
            if (e.getButton() == MouseButton.SECONDARY){
                // panning
                offsetX += (e.getX() - lastMouseX);
                offsetY += (e.getY() - lastMouseY);
                lastMouseX = e.getX();
                lastMouseY = e.getY();
                draw();
                return;
            }

            if (draggingNode != null){
                draggingNode.x = world.getX();
                draggingNode.y = world.getY();
                graph.recomputeWeightsEuclidean();
                draw();
            }
        });

        canvas.setOnMouseReleased(e -> draggingNode = null);

        canvas.setOnScroll(e -> {
            double factor = (e.getDeltaY() > 0) ? 1.1 : 0.9;
            Point2D before = screenToWorld(e.getX(), e.getY());
            scale *= factor;
            // keep mouse point stable
            Point2D after = screenToWorld(e.getX(), e.getY());
            offsetX += (after.getX() - before.getX()) * scale;
            offsetY += (after.getY() - before.getY()) * scale;
            draw();
        });

        canvas.setOnMouseMoved(e -> {
            Point2D w = screenToWorld(e.getX(), e.getY());
            g.setFont(Font.font(12));
            // optional overlay tooltip handled in draw()
        });
    }

    private void runACO(){
        try {
            int ants = Integer.parseInt(antsField.getText().trim());
            int iters = Integer.parseInt(itersField.getText().trim());
            double alpha = Double.parseDouble(alphaField.getText().trim());
            double beta = Double.parseDouble(betaField.getText().trim());
            double rho = Double.parseDouble(rhoField.getText().trim());
            double Q = Double.parseDouble(qField.getText().trim());
            int src = Integer.parseInt(srcField.getText().trim());
            int dst = Integer.parseInt(dstField.getText().trim());

            if (graph.nodes.isEmpty()){
                status("Graf kosong. Tambahkan vertex dahulu.");
                return;
            }
            if (src < 0 || src >= graph.nodes.size() || dst < 0 || dst >= graph.nodes.size()){
                status("Source/Destination di luar rentang node.");
                return;
            }

            double[][] A = graph.toAdjacency();
            ACO aco = new ACO(A, src, dst, ants, iters, alpha, beta, rho, Q);
            ACO.Result res = aco.solve();

            bestPath = (res.bestPath == null) ? Collections.emptyList() : res.bestPath;
            bestCost = res.bestCost;

            if (bestPath.isEmpty() || bestCost >= Double.MAX_VALUE/8){
                status("Tidak ditemukan rute dari " + src + " ke " + dst + ".");
            } else {
                status("Best path: " + bestPath + " | cost ≈ " + String.format(Locale.US, "%.2f", bestCost));
            }
            draw();
        } catch (NumberFormatException ex){
            status("Parameter tidak valid: " + ex.getMessage());
        } catch (Exception ex){
            status("Error: " + ex.getMessage());
            ex.printStackTrace();
        }
    }

    private void autoLayout(){
        int n = graph.nodes.size();
        if (n == 0) return;
        double cx = 0, cy = 0;
        for (Node nd: graph.nodes){ cx += nd.x; cy += nd.y; }
        cx /= n; cy /= n;
        double R = 150 + 10*n;
        for (int i=0;i<n;i++){
            double ang = 2*Math.PI * i / n;
            graph.nodes.get(i).x = cx + R*Math.cos(ang);
            graph.nodes.get(i).y = cy + R*Math.sin(ang);
        }
        graph.recomputeWeightsEuclidean();
    }

    // ====== Rendering ======
    private void draw(){
        GraphicsContext g = canvas.getGraphicsContext2D();
        double w = canvas.getWidth(), h = canvas.getHeight();

        // background
        g.setFill(Color.web("#1e1f22"));
        g.fillRect(0,0,w,h);

        // grid (world space)
        drawGrid(g);

        // edges
        g.setLineWidth(2.0);
        g.setStroke(Color.web("#5e6875"));
        for (Edge e: graph.edges){
            Node a = graph.nodes.get(e.a), b = graph.nodes.get(e.b);
            Point2D sa = worldToScreen(a.x, a.y);
            Point2D sb = worldToScreen(b.x, b.y);
            g.setStroke(Color.web("#5e6875"));
            g.strokeLine(sa.getX(), sa.getY(), sb.getX(), sb.getY());

            // weight label mid
            double mx = (sa.getX() + sb.getX())/2.0;
            double my = (sa.getY() + sb.getY())/2.0;
            g.setFill(Color.web("#aab2bf"));
            g.setFont(Font.font(11));
            g.fillText(String.format(Locale.US, "%.1f", e.weight), mx+4, my-4);
        }

        // highlight best path
        if (bestPath != null && bestPath.size() >= 2 && bestCost < Double.MAX_VALUE/8){
            g.setLineWidth(4.0);
            g.setStroke(Color.web("#4cc9f0"));
            for (int i=0;i<bestPath.size()-1;i++){
                Node a = graph.nodes.get(bestPath.get(i));
                Node b = graph.nodes.get(bestPath.get(i+1));
                Point2D sa = worldToScreen(a.x, a.y);
                Point2D sb = worldToScreen(b.x, b.y);
                g.strokeLine(sa.getX(), sa.getY(), sb.getX(), sb.getY());
            }
        }

        // nodes
        double R = 10;
        for (Node nd: graph.nodes){
            Point2D s = worldToScreen(nd.x, nd.y);
            // fill
            g.setFill(Color.web("#ffd166"));
            g.fillOval(s.getX()-R, s.getY()-R, 2*R, 2*R);
            g.setStroke(Color.web("#2d2f33"));
            g.setLineWidth(2);
            g.strokeOval(s.getX()-R, s.getY()-R, 2*R, 2*R);
            // label
            g.setFill(Color.WHITE);
            g.setFont(Font.font(13));
            g.fillText(String.valueOf(nd.id), s.getX()+R+3, s.getY()-R-3);
        }

        // HUD
        g.setFill(Color.web("#ffffffaa"));
        g.setFont(Font.font(12));
        g.fillText(String.format(Locale.US, "Nodes: %d  Edges: %d  Scale: %.2f",
                graph.nodes.size(), graph.edges.size(), scale), 12, h - 14);
    }

    private void drawGrid(GraphicsContext g){
        // world rect visible
        Point2D topLeft = screenToWorld(0,0);
        Point2D botRight = screenToWorld(canvas.getWidth(), canvas.getHeight());
        double minX = Math.min(topLeft.getX(), botRight.getX());
        double maxX = Math.max(topLeft.getX(), botRight.getX());
        double minY = Math.min(topLeft.getY(), botRight.getY());
        double maxY = Math.max(topLeft.getY(), botRight.getY());

        double step = niceGridStep(50/scale); // roughly 50 px
        double startX = Math.floor(minX/step)*step;
        double startY = Math.floor(minY/step)*step;

        g.setStroke(Color.web("#2b2f36"));
        g.setLineWidth(1.0);
        for (double x = startX; x <= maxX; x += step){
            Point2D s1 = worldToScreen(x, minY);
            Point2D s2 = worldToScreen(x, maxY);
            g.strokeLine(s1.getX(), s1.getY(), s2.getX(), s2.getY());
        }
        for (double y = startY; y <= maxY; y += step){
            Point2D s1 = worldToScreen(minX, y);
            Point2D s2 = worldToScreen(maxX, y);
            g.strokeLine(s1.getX(), s1.getY(), s2.getX(), s2.getY());
        }
    }

    private double niceGridStep(double target){
        double[] steps = {1,2,5};
        double scale = Math.pow(10, Math.floor(Math.log10(target)));
        for (double s: steps){
            double cand = s*scale;
            if (cand >= target) return cand;
        }
        return 10*scale;
    }

    // ====== Utils koordinat ======
    private Point2D worldToScreen(double x, double y){
        return new Point2D(x * scale + offsetX, y * scale + offsetY);
    }
    private Point2D screenToWorld(double sx, double sy){
        return new Point2D((sx - offsetX) / scale, (sy - offsetY) / scale);
    }

    private Node findNodeAt(double x, double y){
        double r = 12 / scale; // hit radius in world units
        for (int i=graph.nodes.size()-1; i>=0; i--){
            Node nd = graph.nodes.get(i);
            if (Math.hypot(nd.x - x, nd.y - y) <= r) return nd;
        }
        return null;
    }

    private void status(String msg){
        statusLabel.setText(msg);
    }

    public static void main(String[] args){
        launch(args);
    }
}
