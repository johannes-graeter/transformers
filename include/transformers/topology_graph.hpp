// Copyright 2019. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@kit.edu)
//  Johannes Beck (johannes.beck@kit.edu)

#pragma once
#include <iostream>
#include <variant>
#include <boost/bimap.hpp>
//#include <boost/variant.hpp>
#include <variant>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/range/irange.hpp>

#include "topology_path.hpp"
#include "internal/tf_traits.hpp"

namespace ct {

/**
 * @brief The FoundGoalException struct
 * Exception that is thrown when graph traversal found the goal.
 * This is the method to terminate f.e. breadth first search descrobed by boost.
 */
struct FoundGoalException : std::exception {
    const char* what() const noexcept override {
        return "Found goal in graph search. Catch this!";
    }
};


/**
 * @brief The EndVisitor class
 * Visitor for graph traversals, to end when goal is found.
 */
template <typename EventType, typename Vertex>
class EndVisitor : public boost::default_bfs_visitor {
public:
    using event_filter = EventType; // NOLINT since boost needs event_filter

public:
    explicit EndVisitor(Vertex goal) : goal_(goal) {
    }
    template <typename OtherVertex, typename Graph>
    void operator()(OtherVertex u, const Graph& /*g*/) const {
        if (u == goal_) {
            throw FoundGoalException();
        }
    }

private:
    Vertex goal_;
};


/**
 * @class TopologyGraph
 * @tparam Ts... Templates that must be dynamic coordinate systems types (TfTraits are defined and dynamic).
 * @todo What to do with static types?
 *
 * Defines and connects topologies for coordinate transforms
 */
template <typename VarType>
class TopologyGraph {
    using Node = size_t;
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

public:
    template <typename T0, typename T1>
    void addEdge(const T0& cs0, const T1& cs1) {
        Node node0 = addNode(cs0);
        Node node1 = addNode(cs1);

        boost::add_edge(node0, node1, graph_);
    }


    /**
     * @brief get Get topology from graph. Use breadth first search since edges have same weight.
     * @param to coord systems to transform to.
     * @param from coord system from which shall be transfered.
     * @param time time instance to search at.
     * @return Path from - to at time t as Topology.
     */
    template <typename ToCoord, typename FromCoord>
    TopologyPath<VarType> get(const ToCoord& to, const FromCoord& from, Time time = Time(0)) const {


        // This throws if node where not added yet.
        Node toNode = nodeList_.left.at(VarType(TfTrait<ToCoord>::toDynamic(to)));
        Node fromNode = nodeList_.left.at(VarType(TfTrait<FromCoord>::toDynamic(from)));

        // Convert to boost types.
        using Vertex = Graph::vertex_descriptor;
        Vertex start = boost::vertex(fromNode, graph_);
        Vertex goal = boost::vertex(toNode, graph_);

        // Instantiate visitor that throws exception when goal reached (yes this is what boost suggests...)
        EndVisitor<boost::on_discover_vertex, Vertex> endVisitor(goal);
        // Array to store predecessor (parent) of each vertex. This will be used as a decorator (actually, its iterator
        // will be).
        std::vector<Vertex> p(boost::num_vertices(graph_));
        p[start] = start;

        // Make combined visitor (visitors must implement operator() and have event_filter type).
        auto vis = boost::visitor(boost::make_bfs_visitor(
            std::make_pair(boost::record_predecessors(&p[0], boost::on_tree_edge()), endVisitor)));

        // Do the search.
        try {
            boost::breadth_first_search(graph_, start, vis);
        } catch (const FoundGoalException& /*e*/) {
        }

        // Walk back on predecessor_map and add to topology.
        TopologyPath<VarType> out;
        Vertex current = goal;
        while (current != start) {
            auto it = nodeList_.right.find(current);
            out.add(it->second, time);
            current = p[current];
        }
        auto it = nodeList_.right.find(start);
        out.add(it->second, time);

        return out;
    }

    void print() const {
        boost::print_graph(graph_);
    }


private:
    template <typename OtherT>
    Node addNode(const OtherT& cs) {
        // Cast cooordinate system to dynamic, so that we do not have to give all static coords in the Ts.
        // Dynamic type must be one of the Ts.
        VarType node = TfTrait<OtherT>::toDynamic(cs);

        // Add to node list.
        Node out;
        auto it = nodeList_.left.find(node);
        if (it == nodeList_.left.end()) { // Does not exist.
            Node ind = nodeList_.size();
            out = ind;
            nodeList_.insert({node, ind});
        } else { // Does exist.
            out = it->second;
        }
        return out;
    }

    ///@brief Nodes are converted to ints(=enumeration) for comparison, store them here.
    boost::bimap<VarType, size_t> nodeList_;

    ///@brief Adjacency list to save topology graph.
    Graph graph_;
};
} // namespace ct
