/*
 * Copyright 2017 HugeGraph Authors
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

package com.baidu.hugegraph.example;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.tinkerpop.gremlin.process.traversal.P;
import org.apache.tinkerpop.gremlin.process.traversal.Path;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversal;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversalSource;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.__;
import org.apache.tinkerpop.gremlin.structure.Edge;
import org.apache.tinkerpop.gremlin.structure.T;
import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.slf4j.Logger;

import com.baidu.hugegraph.HugeGraph;
import com.baidu.hugegraph.backend.id.IdGenerator;
import com.baidu.hugegraph.schema.EdgeLabel;
import com.baidu.hugegraph.schema.SchemaManager;
import com.baidu.hugegraph.schema.VertexLabel;
import com.baidu.hugegraph.structure.HugeEdge;
import com.baidu.hugegraph.structure.HugeVertex;
import com.baidu.hugegraph.type.define.Directions;
import com.baidu.hugegraph.util.E;
import com.baidu.hugegraph.util.Log;
import com.google.common.collect.ImmutableSet;

public class Example2 {

    private static final Logger LOG = Log.logger(Example2.class);

    public static void main(String[] args) throws InterruptedException {
        LOG.info("Example2 start!");

        HugeGraph graph = ExampleUtil.loadGraph();

        Example2.load(graph);

        int maxDepth = 25;
//        Map<Object, Double> rank = Example2.personRank(graph, IdGenerator.of("A"), 0.5, 3);
        Map<Object, Double> rank = Example2.neighborRank2(graph, IdGenerator.of("A"), "rate", 0.5, maxDepth);
        // 根据值大小排序
        StringBuilder result = new StringBuilder();
        result.append("| ").append(maxDepth).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("A"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("B"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("C"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("a"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("b"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("c"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append(((int) (rank.getOrDefault(IdGenerator.of("d"), 0.0) * 1000)) / 1000.0).append(" | ");
        result.append("1.0").append(" |");
        System.out.println(result);

        System.out.println(sortByValue(rank));
//        traversal(graph);
        graph.close();

        HugeGraph.shutdown(30L);
    }

    public static Map<Object, Double> personRank(final HugeGraph graph,
                                                 Object vertexId,
                                                 double alpha,
                                                 int iterCount) {
        long beg = System.currentTimeMillis();

        GraphTraversalSource g = graph.traversal();

        Map<Object, Double> rank = makeEmptyRank(g);
        // 起点选择概率为1,其他顶点为0
        rank.put(vertexId, 1.0);

        // 开始迭代
        for (int i = 0; i < iterCount; i++) {
            Map<Object, Double> tmp = makeEmptyRank(g);
            // 遍历每一个顶点
            GraphTraversal<Vertex, Vertex> vertices = g.V();
            while (vertices.hasNext()) {
                Vertex vertex = vertices.next();
                // 遍历每个顶点连接的顶点
                List<Vertex> boths = g.V(vertex.id()).out().toList();
                for (Vertex bothVertex : boths) {
                    double weight = tmp.get(bothVertex.id());
                    weight += alpha * rank.get(vertex.id()) / boths.size();
                    tmp.put(bothVertex.id(), weight);
                }
            }
            double rootWeight = tmp.get(vertexId);
            rootWeight += 1 - alpha;
            tmp.put(vertexId, rootWeight);

            rank = tmp;
        }

        long end = System.currentTimeMillis();
        System.out.println("耗时：" + (end - beg) / 1000.0F);
        System.out.println("总和：" + rank.values().stream().reduce((a, b) -> a + b).get());
        return rank;
    }

    public static Map<Object, Double> neighborRank2(HugeGraph graph,
                                                   Object vertexId,
                                                   String label,
                                                   double alpha,
                                                   int maxDepth) {
        E.checkArgument(maxDepth >= 1, "Must > 1");
        long beg = System.currentTimeMillis();

        GraphTraversalSource g = graph.traversal();
        // TODO: check exist
        HugeVertex vertex = (HugeVertex) g.V(vertexId).next();
        VertexLabel vertexLabel = vertex.schemaLabel();
        HugeEdge edge = (HugeEdge) g.V(vertexId).bothE(label).next();
        EdgeLabel edgeLabel = edge.schemaLabel();

//        Map<String, Directions> labelDirs = new HashMap<>();
//        if (edgeLabel.sourceLabel().equals(vertexLabel.id())) {
//            labelDirs.put(graph.vertexLabel(edgeLabel.sourceLabel()).name(), Directions.OUT);
//            labelDirs.put(graph.vertexLabel(edgeLabel.targetLabel()).name(), Directions.IN);
//        } else {
//            labelDirs.put(graph.vertexLabel(edgeLabel.sourceLabel()).name(), Directions.IN);
//            labelDirs.put(graph.vertexLabel(edgeLabel.targetLabel()).name(), Directions.OUT);
//        }

        Map<Object, Double> rank = new HashMap<>();
        rank.put(vertex.id(), 1.0);
        rank = collectWeight2(graph, vertexId, label, 1, rank, alpha, maxDepth);

        long end = System.currentTimeMillis();
        System.out.println("耗时：" + (end - beg) / 1000.0F);
        System.out.println("总和：" + rank.values().stream().reduce((a, b) -> a + b).get());
        return rank;
    }

    /**
     * 我自己的算法 both 版
     * @param graph
     * @param depth
     * @param rank
     * @param alpha
     * @param maxDepth
     * @return
     */
    private static Map<Object, Double> collectWeight2(HugeGraph graph,
                                                      Object source,
                                                      String edgeLabel,
//                                                      Map<String, Directions> labelDirs,
                                                      int depth,
                                                      Map<Object, Double> rank,
                                                      double alpha,
                                                      int maxDepth) {
        assert depth <= maxDepth;
        GraphTraversalSource g = graph.traversal();

        // 收集当前节点和其指定方向和label的邻居节点
        Set<Vertex> vertices = new HashSet<>();
        for (Object vid : rank.keySet()) {
            vertices.add(g.V(vid).next());
            vertices.addAll(g.V(vid).both(edgeLabel).toSet());
        }

        Map<Object, Double> tmpRank = new HashMap<>();
        for (Vertex vertex : vertices) {
            Object vertexId = vertex.id();
            // 自己对自己的贡献
            double nextWeight = tmpRank.getOrDefault(vertexId, 0.0);
//            double nextWeight = selfContribution(graph, vertexId, edgeLabel, rank, alpha);

            // 它的所有邻居顶点对自己的贡献
            GraphTraversal<Vertex, Vertex> neighborVertices = g.V(vertexId).both(edgeLabel);
            while (neighborVertices.hasNext()) {
                Vertex neighborVertex = neighborVertices.next();
                Double neighborVertexWeight = rank.get(neighborVertex.id());
                // 如果当前顶点还不存在权重，则跳过
                if (neighborVertexWeight == null) {
                    continue;
                }
                long degree = g.V(neighborVertex.id()).both().count().next();
//                Directions dir = labelDirs.get(neighborVertex.label());
//                if (dir == Directions.OUT) {
//                    // 计算顶点的出度
//                    degree = g.V(neighborVertex.id()).out().count().next();
//                } else {
//                    assert dir == Directions.IN;
//                    // 计算顶点的入度
//                    degree = g.V(neighborVertex.id()).in().count().next();
//                }
                assert degree > 0;
                nextWeight += neighborVertexWeight * alpha / degree;
            }
            // 更新顶点的权重
            tmpRank.put(vertexId, nextWeight);
        }
        double sourceRank = tmpRank.get(source);
        sourceRank += (1 - alpha);
        tmpRank.put(source, sourceRank);

        // 如果已经达到指定步数，则返回，否则继续
        if (++depth > maxDepth) {
            return tmpRank;
        }
        return collectWeight2(graph, source, edgeLabel, depth, tmpRank, alpha, maxDepth);
    }

    private static double selfContribution(HugeGraph graph, Object vertexId, String edgeLabel,
                                           Map<Object, Double> rank, double alpha) {
        GraphTraversalSource g = graph.traversal();

        // 它自己对自己的贡献
        Double lastWeight = rank.get(vertexId);
        if (lastWeight == null) {
            lastWeight = 0.0;
        }
        double nextWeight = lastWeight * (1 - alpha);
        long outCount = g.V(vertexId).both(edgeLabel).count().next();
        if (outCount == 0) {
            nextWeight = lastWeight;
        }
        return nextWeight;
    }


    public static Map<Object, Double> neighborRank3(HugeGraph graph,
                                                    Object vertexId,
                                                    String label,
                                                    double alpha,
                                                    int maxDepth) {
        E.checkArgument(maxDepth >= 1, "Must > 1");
        long beg = System.currentTimeMillis();

        GraphTraversalSource g = graph.traversal();
        // TODO: check exist
        HugeVertex vertex = (HugeVertex) g.V(vertexId).next();

        Map<Object, Double> rank = new HashMap<>();
        rank.put(vertex.id(), 1.0);
        rank = collectWeight3(graph, ImmutableSet.of(vertexId), label, 1, rank, alpha, maxDepth);

        long end = System.currentTimeMillis();
        System.out.println("耗时：" + (end - beg) / 1000.0F);
        System.out.println("总和：" + rank.values().stream().reduce((a, b) -> a + b).get());
        return rank;
    }

    private static Map<Object, Double> collectWeight3(HugeGraph graph,
                                                      Set<Object> vertices,
                                                      String edgeLabel,
                                                      int depth,
                                                      Map<Object, Double> rank,
                                                      double alpha,
                                                      int maxDepth) {
        assert depth <= maxDepth;
        GraphTraversalSource g = graph.traversal();

        Set<Object> tmpVertices = new HashSet<>();
        Map<Object, Double> tmpRank = new HashMap<>();
        // 外层做贡献的顶点
        for (Object vertexId : vertices) {
            long degree = g.V(vertexId).bothE(edgeLabel).count().next();
            assert degree > 0;

            // 这个应该肯定存在的
            double vertexWeight = rank.get(vertexId);
            double spreadWeight = vertexWeight * alpha / degree;

            GraphTraversal<Vertex, Vertex> neighborVertices = g.V(vertexId).both(edgeLabel);
            while (neighborVertices.hasNext()) {
                Vertex neighborVertex = neighborVertices.next();
                Double neighborVertexWeight = tmpRank.get(neighborVertex.id());
                // 如果当前顶点还不存在权重，则赋上初始值
                if (neighborVertexWeight == null) {
                    neighborVertexWeight = 0.0;
                }
                // 目标顶点加上被贡献的值（不应该原地更新）
                neighborVertexWeight += spreadWeight;
                tmpVertices.add(neighborVertex.id());
                tmpRank.put(neighborVertex.id(), neighborVertexWeight);
            }
        }

        // 源顶点直接乘以(1 - alpha)即可
        for (Object vertexId : vertices) {
            double vertexWeight = rank.get(vertexId);
            vertexWeight *= (1 - alpha);
            rank.put(vertexId, vertexWeight);
        }
        // 目标顶点加上原有的值
        for (Map.Entry<Object, Double> entry : tmpRank.entrySet()) {
            Double oldWeight = rank.get(entry.getKey());
            if (oldWeight == null) {
                oldWeight = 0.0;
            }
            double increment = entry.getValue();
            rank.put(entry.getKey(), oldWeight + increment);
        }
        // 如果已经达到指定步数，则返回，否则继续
        if (++depth > maxDepth) {
            return rank;
        }
        return collectWeight3(graph, tmpVertices, edgeLabel, depth, rank, alpha, maxDepth);
    }

    private static Map<Object, Double> makeEmptyRank(GraphTraversalSource g) {
        Map<Object, Double> rank = new HashMap<>();
        GraphTraversal<Vertex, Vertex> vertices = g.V();
        // 将所有的点的权重都赋值为0
        while (vertices.hasNext()) {
            Vertex vertex = vertices.next();
            rank.put(vertex.id(), 0.0);
        }
        return rank;
    }

    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Map.Entry<K, V>> list = new ArrayList<>(map.entrySet());
        list.sort(Collections.reverseOrder(Map.Entry.comparingByValue()));

        Map<K, V> result = new LinkedHashMap<>();
        for (Map.Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }

        return result;
    }

    public static void traversal(final HugeGraph graph) {

        GraphTraversalSource g = graph.traversal();

        GraphTraversal<Vertex, Vertex> vertexs = g.V();
        System.out.println(">>>> query all vertices: size=" +
                           vertexs.toList().size());

        List<Edge> edges = g.E().toList();
        System.out.println(">>>> query all edges: size=" +
                           edges.size());

        List<Object> names = g.V().inE("knows").limit(2)
                              .outV().values("name").toList();
        System.out.println(">>>> query vertex(with props) of edges: " + names);
        assert names.size() == 2 : names.size();

        names = g.V().as("a")
                 .out("knows")
                 .and()
                 .out("created").in("created").as("a").values("name")
                 .toList();
        System.out.println(">>>> query with AND: " + names);
        assert names.size() == 1 : names.size();

        List<Vertex> vertex = g.V().has("age", 29).toList();
        System.out.println(">>>> age = 29: " + vertex);
        assert vertex.size() == 1 &&
               vertex.get(0).value("name").equals("marko");

        vertex = g.V()
                  .has("age", 29)
                  .has("city", "Beijing")
                  .toList();
        System.out.println(">>>> age = 29 and city is Beijing: " + vertex);
        assert vertex.size() == 1 &&
               vertex.get(0).value("name").equals("marko");

        edges = g.E().has("weight", P.lt(1.0)).toList();
        System.out.println(">>>> edges with weight < 1.0: " + edges);
        assert edges.size() == 4;

        String person = graph.schema().getVertexLabel("person").id().asString();
        String software = graph.schema().getVertexLabel("software").id()
                          .asString();
        String markoId = String.format("%s:%s", person, "marko");
        String joshId = String.format("%s:%s", person, "josh");
        String lopId = String.format("%s:%s", software, "lop");

        vertex = g.V(joshId)
                  .bothE("created")
                  .has("weight", P.lt(1.0))
                  .otherV()
                  .toList();
        System.out.println(">>>> josh's both edges with weight < 1.0: " +
                           vertex);
        assert vertex.size() == 1 && vertex.get(0).value("name").equals("lop");

        List<Path> paths = g.V(markoId).out().out().path().by("name").toList();
        System.out.println(">>>> test out path: " + paths);
        assert paths.size() == 2;
        assert paths.get(0).get(0).equals("marko");
        assert paths.get(0).get(1).equals("josh");
        assert paths.get(0).get(2).equals("lop");
        assert paths.get(1).get(0).equals("marko");
        assert paths.get(1).get(1).equals("josh");
        assert paths.get(1).get(2).equals("ripple");

        paths = shortestPath(graph, markoId, lopId, 5);
        System.out.println(">>>> test shortest path: " + paths.get(0));
        assert paths.get(0).get(0).equals("marko");
        assert paths.get(0).get(1).equals("lop");
    }

    public static List<Path> shortestPath(final HugeGraph graph,
                                          Object from, Object to,
                                          int maxDepth) {
        GraphTraversal<Vertex, Path> t = graph.traversal()
                .V(from)
                .repeat(__.out().simplePath())
                .until(__.hasId(to).or().loops().is(P.gt(maxDepth)))
                .hasId(to)
                .path().by("name")
                .limit(1);
        return t.toList();
    }

    public static void load(final HugeGraph graph) {
        SchemaManager schema = graph.schema();

        schema.propertyKey("name").asText().ifNotExist().create();
        schema.propertyKey("weight").asDouble().ifNotExist().create();

        schema.vertexLabel("person")
              .useCustomizeStringId()
              .ifNotExist()
              .create();

        schema.vertexLabel("movie")
              .useCustomizeStringId()
              .ifNotExist()
              .create();

        schema.edgeLabel("rate")
              .sourceLabel("person")
              .targetLabel("movie")
              .properties("weight")
              .ifNotExist()
              .create();

        Vertex A = graph.addVertex(T.label, "person", T.id, "A");
        Vertex B = graph.addVertex(T.label, "person", T.id, "B");
        Vertex C = graph.addVertex(T.label, "person", T.id, "C");

        Vertex a = graph.addVertex(T.label, "movie", T.id, "a");
        Vertex b = graph.addVertex(T.label, "movie", T.id, "b");
        Vertex c = graph.addVertex(T.label, "movie", T.id, "c");
        Vertex d = graph.addVertex(T.label, "movie", T.id, "d");

        A.addEdge("rate", a, "weight", 1.0);
        A.addEdge("rate", c, "weight", 1.0);

        B.addEdge("rate", a, "weight", 1.0);
        B.addEdge("rate", b, "weight", 1.0);
        B.addEdge("rate", c, "weight", 1.0);
        B.addEdge("rate", d, "weight", 1.0);

        C.addEdge("rate", c, "weight", 1.0);
        C.addEdge("rate", d, "weight", 1.0);

        graph.tx().commit();
    }

    public static void load2(final HugeGraph graph) {
        SchemaManager schema = graph.schema();

        schema.propertyKey("name").asText().ifNotExist().create();
        schema.propertyKey("weight").asDouble().ifNotExist().create();

        schema.vertexLabel("person")
              .useCustomizeStringId()
              .ifNotExist()
              .create();

        schema.edgeLabel("follow")
              .sourceLabel("person")
              .targetLabel("person")
              .properties("weight")
              .ifNotExist()
              .create();

        Vertex A = graph.addVertex(T.label, "person", T.id, "A");
        Vertex B = graph.addVertex(T.label, "person", T.id, "B");
        Vertex C = graph.addVertex(T.label, "person", T.id, "C");
        Vertex D = graph.addVertex(T.label, "person", T.id, "D");
        Vertex F = graph.addVertex(T.label, "person", T.id, "F");
        Vertex G = graph.addVertex(T.label, "person", T.id, "G");
        Vertex V = graph.addVertex(T.label, "person", T.id, "V");
        Vertex S = graph.addVertex(T.label, "person", T.id, "S");

        A.addEdge("follow", C, "weight", 1.0);
        A.addEdge("follow", F, "weight", 1.0);

        B.addEdge("follow", C, "weight", 1.0);
        B.addEdge("follow", V, "weight", 1.0);

        C.addEdge("follow", V, "weight", 1.0);

        D.addEdge("follow", V, "weight", 1.0);

        F.addEdge("follow", D, "weight", 1.0);
        F.addEdge("follow", G, "weight", 1.0);
        F.addEdge("follow", V, "weight", 1.0);

        G.addEdge("follow", A, "weight", 1.0);

        V.addEdge("follow", S, "weight", 1.0);

        graph.tx().commit();
    }
}
