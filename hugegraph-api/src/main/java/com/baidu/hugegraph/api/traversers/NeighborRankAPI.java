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

package com.baidu.hugegraph.api.traversers;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.inject.Singleton;
import javax.ws.rs.DefaultValue;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.QueryParam;
import javax.ws.rs.core.Context;

import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversal;
import org.apache.tinkerpop.gremlin.process.traversal.dsl.graph.GraphTraversalSource;
import org.apache.tinkerpop.gremlin.structure.Vertex;
import org.slf4j.Logger;

import com.baidu.hugegraph.HugeGraph;
import com.baidu.hugegraph.api.API;
import com.baidu.hugegraph.api.graph.VertexAPI;
import com.baidu.hugegraph.backend.id.Id;
import com.baidu.hugegraph.core.GraphManager;
import com.baidu.hugegraph.server.RestServer;
import com.baidu.hugegraph.util.CollectionUtil;
import com.baidu.hugegraph.util.E;
import com.baidu.hugegraph.util.JsonUtil;
import com.baidu.hugegraph.util.Log;
import com.codahale.metrics.annotation.Timed;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;

@Path("graphs/{graph}/traversers/neighborrank")
@Singleton
public class NeighborRankAPI extends API {

    private static final Logger LOG = Log.logger(RestServer.class);

    @GET
    @Timed
    @Produces(APPLICATION_JSON_WITH_CHARSET)
    public String neighborRank(@Context GraphManager manager,
                               @PathParam("graph") String graph,
                               @QueryParam("source") String source,
                               @QueryParam("label") String edgeLabel,
                               @QueryParam("alpha") double alpha,
                               @QueryParam("max_depth") int maxDepth,
                               @QueryParam("sort_result") @DefaultValue("true")
                               boolean sortResult) {
        LOG.debug("Graph [{}] get neighbor rank from '{}' with " +
                  "edge label '{}', alpha '{}', max depth '{}'",
                  graph, source, edgeLabel, alpha, maxDepth);

        E.checkNotNull(source, "source vertex id");
        E.checkNotNull(edgeLabel, "edge label");
        E.checkArgument(alpha >= 0.0 && alpha <= 1.0,
                        "The alpha must between [0, 1], but got '%s'", alpha);
        E.checkArgument(maxDepth >= 1,
                        "The max depth must >= 1, but got '%s'", maxDepth);

        Id sourceId = VertexAPI.checkAndParseVertexId(source);
        HugeGraph g = graph(manager, graph);

        long degree = g.traversal().V(source).bothE(edgeLabel).count().next();
        if (degree < 1) {
            return JsonUtil.toJson(ImmutableMap.of(sourceId, "1.0"));
        }

        Set<Object> seedVertices = ImmutableSet.of(sourceId);
        Map<Object, Double> rank = new HashMap<>();
        rank.put(sourceId, 1.0);
        rank = neighborRank(g, seedVertices, rank, edgeLabel, alpha, 1, maxDepth);
        if (sortResult) {
            rank = CollectionUtil.sortByValue(rank, false);
        }
        return JsonUtil.toJson(rank);
    }

    private static Map<Object, Double> neighborRank(HugeGraph graph,
                                                    Set<Object> seedVertices,
                                                    Map<Object, Double> rank,
                                                    String edgeLabel,
                                                    double alpha,
                                                    int depth,
                                                    int maxDepth) {
        assert depth <= maxDepth;
        GraphTraversalSource g = graph.traversal();

        Set<Object> neighborVertices = new HashSet<>();
        Map<Object, Double> rankIncrs = new HashMap<>();
        for (Object vertexId : seedVertices) {
            long degree = g.V(vertexId).bothE(edgeLabel).count().next();
            assert degree > 0;
            // Must be exist
            double originRank = rank.get(vertexId);
            double spreadRank = originRank * alpha / degree;

            GraphTraversal<Vertex, Vertex> neighbors = g.V(vertexId)
                                                        .both(edgeLabel);
            // Collect all neighbors increment
            while (neighbors.hasNext()) {
                Vertex neighbor = neighbors.next();
                /*
                 * Assign an initial value when firstly update
                 * neighbor vertex rank
                 */
                double neighborRankIncr = rankIncrs.getOrDefault(neighbor.id(),
                                                                 0.0);
                neighborRankIncr += spreadRank;
                neighborVertices.add(neighbor.id());
                rankIncrs.put(neighbor.id(), neighborRankIncr);
            }
        }
        // Update the rank of the neighbor vertices
        for (Map.Entry<Object, Double> entry : rankIncrs.entrySet()) {
            double oldRank = rank.getOrDefault(entry.getKey(), 0.0);
            double incr = entry.getValue();
            rank.put(entry.getKey(), oldRank + incr);
        }
        // The rank of each seed vertex is attenuated to original * (1 - alpha)
        for (Object vertexId : seedVertices) {
            double oldRank = rank.get(vertexId) * (1 - alpha);
            rank.put(vertexId, oldRank);
        }

        if (++depth > maxDepth) {
            return rank;
        }
        return neighborRank(graph, neighborVertices, rank, edgeLabel,
                            alpha, depth, maxDepth);
    }
}