// Copyright 2018. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@kit.edu)
//  and others

#pragma once
#include <functional>
#include <map>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "save_coordinate_transform.hpp"
#include "topology_path.hpp"

namespace ct {

/**
 *@brief Traits for casting trasnform types.
 * This is necessary if you want to use f.e. Jet types (ceres-solver)
 */
template <typename T>
struct CastTrait;

template <>
struct CastTrait<Eigen::Isometry3d> {
    template <typename ScalarType>
    static Eigen::Transform<ScalarType, 3, Eigen::Isometry> apply(const Eigen::Isometry3d& transform) {
        return transform.cast<ScalarType>;
    }
};

/**
 * @brief The Storage struct
 * Container for transforms.
 * TransformType must defince .inverse()
 * @todo solve inverse with traits.
 */
template <typename TransformType>
struct Storage {
    /**
     * @brief get, get transform by index. If index is negative, inverse is returned.
     * @param index
     * @return Transform
     */
    TransformType get(int index) const {
        if (index > 0) {
            return data_.at(index - 1);
        } else if (index < 0) {
            return data_.at(std::abs(index) - 1).inverse();
        } else {
            throw std::runtime_error("Transformationgraph::Storage is 1 indicated.");
        }
    }

    void pushBack(const TransformType& t) {
        data_.push_back(t);
    }

    size_t size() const {
        return data_.size();
    }

    void clear() {
        data_.clear();
    }

private:
    std::vector<TransformType, Eigen::aligned_allocator<TransformType>> data_; ///> Data must implement "inverse()"
};

/**
 * @brief The DummyDynamic struct
 * This class is used for Dynamic transforms, if the TransformationGraph is purely static.
 */
template <typename StaticTransform>
struct DummyDynamic {
    template <typename A, typename B>
    StaticTransform operator()(const A& /*a*/, const B& /*b*/) const {
        return StaticTransform{};
    }
    DummyDynamic inverse() const {
        return DummyDynamic{};
    }
};

/**
 * @class TransformationGraph
 * @param VarType, types of coordinate systems. Usually an std::variant<...>.
 * @param StaticTransformType, type of temporaly static transforms (f.e. Eigen::Isometry3d)
 * @param DynamicTransformType, type of temporaly dynmic transforms (f.e. std::function<Eigen::Isometry3d(Time,Time)>).
 * Must implement operator()(Time, Time).
 * Graph that stores transformations in order to generate a transformation from a topology.
 */
template <typename VarType,
          typename StaticTransformType,
          typename DynamicTransformType = DummyDynamic<StaticTransformType>>
class TransformationGraph {
public: // Public classes/enums/types etc...
    // Specify Eigen Alignment, should be obsolete with c++17
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using CoordType = VarType;
    using StaticStorage = Storage<StaticTransformType>;
    using DynamicStorage = std::vector<DynamicTransformType, Eigen::aligned_allocator<DynamicTransformType>>;
    using Index = int;
    ///@brief Graph to save indices of static and dynamic trasnforms (as adjacency list)
    using Graph = std::map<CoordType, std::map<CoordType, Index>>;

public: // Public methods.
    /**
     * @brief add; add dynamic transform to storage.
     * @param transform; dynamic transform with identical coordinate system.
     */
    template <typename Coord>
    void add(Transform<Coord, Coord, DynamicTransformType> transform) {
        // Convert to dynamic.
        CoordType toCoord = TfTrait<Coord>::toDynamic(transform.getToCoord());
        CoordType fromCoord = TfTrait<Coord>::toDynamic(transform.getFromCoord());
        assert(TfTrait<decltype(toCoord)>::isSame(toCoord, fromCoord));

        // Save index, in graph (1 indicated).
        // Negative index means transform will be inverted.
        // When coordinate system is same, it must be a temporal mutable transform.
        // Add data to storage.
        dynamicTransformations_.push_back(transform.getData());
        graph_[fromCoord][toCoord] = dynamicTransformations_.size();
        graph_[toCoord][fromCoord] = dynamicTransformations_.size();
    }

    /**
     * @brief add; add static transform to storage.
     * @param transform; static transform to be added.
     */
    template <typename ToCoord, typename FromCoord>
    void add(Transform<ToCoord, FromCoord, StaticTransformType> transform) {
        // Convert to dynamic.
        ///@todo How to do that for compile time static coords?
        CoordType toCoord = TfTrait<ToCoord>::toDynamic(transform.getToCoord());
        CoordType fromCoord = TfTrait<FromCoord>::toDynamic(transform.getFromCoord());
        assert(!TfTrait<decltype(toCoord)>::isSame(toCoord, fromCoord));

        // Save index, in graph (1 indicated).
        // Negative index means transform will be inverted.
        // When coordinate system is same, it must be a temporal mutable transform.
        // Add data to storage.
        staticTransformations_.pushBack(transform.getData());
        graph_[fromCoord][toCoord] = -staticTransformations_.size();
        graph_[toCoord][fromCoord] = staticTransformations_.size(); // Order is important if coord system is same.
    }

    /**
     * @brief get, get transform for topology with all static coord systems.
     * @todo how to integrate time? How to integrate variable coords.
     * @param topology
     * @return accumulated transform.
     */
    Transform<VarType, VarType, StaticTransformType> get(const TopologyPath<VarType>& topology) {
        // Get only valid for topolgies with size >1.
        // Alternative: return identity if size 1, define trait for identity of transform.
        if (topology.size() < 2) {
            throw std::runtime_error("Error in TransformGraph::get : Topology needs to have minimum size 2");
        }

        // Iterate pairwise throught the topology, get transforms and accumulate them.
        // First element of topology is the "to" coordinate system.
        auto toIt = topology.getData().cbegin();
        auto fromIt = std::next(topology.getData().cbegin());

        // Initilialize.
        auto ind = graph_.at(toIt->first).at(fromIt->first);
        StaticTransformType out = getTransform(ind, *toIt, *fromIt);
        std::advance(fromIt, 1);
        std::advance(toIt, 1);

        // Accumulate rest.
        for (; fromIt != topology.getData().cend(); fromIt++, toIt++) {
            auto ind = graph_[toIt->first][fromIt->first];
            out = out * getTransform(ind, *toIt, *fromIt);
        }

        // Get output coordinate systems.
        VarType endCoord = topology.getData().front().first;
        VarType startCoord = topology.getData().back().first;

        return Transform<VarType, VarType, StaticTransformType>(out, endCoord, startCoord);
    }

private:
    // ///////// Methods
    /**
     * @brief get_transform, get transform etiher from static of dynamic transformations.
     * @param ind, index of the trasnformation, can also be negative.
     * @param to, coord system and time to transfer to.
     * @param from, coord system and time to transfer from.
     * @return transform, always static.
     */
    StaticTransformType getTransform(const Index& ind,
                                     typename TopologyPath<VarType>::TimedCoord to,
                                     typename TopologyPath<VarType>::TimedCoord from) const {
        StaticTransformType out;
        bool isDynamic = TfTrait<decltype(to.first)>::isSame(to.first, from.first);
        if (isDynamic) {
            // Index is 1 indicated to be consistent with static container.
            ///@todo wrap dynamic transforms in storage?
            auto t = dynamicTransformations_.at(ind - 1);
            out = t(to.second, from.second);
        } else {
            out = staticTransformations_.get(ind);
        }
        return out;
    }

    // ////////// Attributes
    ///@brief Graph to store indices for each pair of coordinate systems.
    /// Negative indices mean inverse transform.
    Graph graph_;
    StaticStorage staticTransformations_;   ///> Storage for temporaly static transforms.
    DynamicStorage dynamicTransformations_; ///> Storage for temporaly dynamic transforms.
};
} // namespace ct
