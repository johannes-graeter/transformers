// Copyright 2018. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@kit.edu)
//  and others

#pragma once

#include <chrono>
#include <functional>
#include <iostream>
#include <utility>
#include <variant>
#include <vector>

namespace ct {

/**
 * @brief Streamer class to print std::variant<Ts...>.
 *  This is fancy c++17 stuff, called "class template argument deduction guide".
 */
template <typename T>
struct Streamer {
    const T& val;
};
template <typename T>
Streamer(T)->Streamer<T>;

template <class T>
std::ostream& operator<<(std::ostream& os, Streamer<T> s) {
    os << s.val;
    return os;
}

template <typename... Ts>
std::ostream& operator<<(std::ostream& os, Streamer<std::variant<Ts...>> sv) {
    std::visit([&os](const auto& v) { os << Streamer{v}; }, sv.val);
    return os;
}

///@brief Define time type.
using Time = std::chrono::nanoseconds;

/**
 * @brief The NotConnectedException struct
 * Exception which is thrown if non connected topologies shall be concatenated.
 */
struct NotConnectedException : public std::exception {
    const char* what() const noexcept override {
        return "Topologies are unconnected.";
    }
};

/**
 * @brief The TopologyPath class
 * Stores a sequence of nodes and time instances.
 * This is used to get the path of a TopologyGraph and use it for transformations.
 * Especially usefull, if paths shall be connected at different time instances.
 */
template <typename VarType>
class TopologyPath {
public:
    using TimedCoord = std::pair<VarType, Time>; ///> Datatype that is stored in topology.

public:
    ///@brief Default construct.
    TopologyPath() = default;

    /**
     * @brief TopologyPath Construct from vector of nodes.
     * @param data Data of nodes of topology path.
     */
    explicit TopologyPath(std::vector<TimedCoord> data) : data_(std::move(data)) {
    }

    /**
     * @brief operator* Concatenation operator for topology path.
     *  Concatenated topology path is optimized so that unnecessary transformations are eliminated.
     * @param t Topology to be appended.
     * @return concatenated topology
     */
    TopologyPath<VarType> operator*(const TopologyPath<VarType>& t) {
        TopologyPath<VarType> out(this->getData()); // Do I need that copy?

        if (out.getData().back().first == t.getData().front().first) {
            // Advance iterators until redundant info is gone.
            auto startIt = t.getData().cbegin();
            auto removeIt = out.getData().crbegin();

            // Compare the ends of the topology that should be connected.
            // If they are teh same advance iterator so that node is present only ones.
            if (out.getData().back().second == t.getData().front().second) {
                std::advance(removeIt, 1);
            }

            // Go through topology and skip identical coords with iterator.
            // In that way those nodes will not be copied to output topology path.
            while (removeIt->second == std::next(startIt)->second && removeIt->first == std::next(startIt)->first &&
                   removeIt != out.getData().crend() && std::next(startIt) != t.getData().cend()) {
                std::advance(removeIt, 1);
                std::advance(startIt, 1);
            }

            int dist = std::distance(removeIt, out.getData().crend());
            out.getData().resize(dist + std::distance(startIt, t.getData().cend()));

            std::copy(startIt, t.getData().cend(), std::next(out.getData().begin(), dist));
        } else {
            throw NotConnectedException();
        }

        return out;
    }

    /**
     * @brief inverse Invert topology in place.
     * @return inverted topology
     */
    TopologyPath<VarType>& inverse() {
        std::reverse(getData().begin(), getData().end());
        return *this;
    }

    /**
     * @brief size, convenience function for getting size of the topology.
     */
    size_t size() const {
        return getData().size();
    }

    /**
     * @brief operator== Equality operator, topologies must be in the same order.
     * @param other
     * @return
     */
    bool operator==(const TopologyPath<VarType>& other) const {
        bool out = true;

        auto dataIter = getData().cbegin();
        auto otherIter = other.getData().cbegin();
        if (getData().size() != other.getData().size()) {
            return false;
        }
        for (; dataIter != getData().cend() && otherIter != other.getData().cend(); ++dataIter, ++otherIter) {
            if (dataIter->first != otherIter->first || dataIter->second != otherIter->second) {
                out = false;
                break;
            }
        }
        return out;
    }

    /**
     * @brief reserve Reserve data of topology.
     * @param n
     */
    void reserve(size_t n) {
        data_.reserve(n);
    }

    /**
     * @brief add Add data to topology
     * @param node coordinate system, that shall be added.
     * @param time time at whcih node shall be added.
     */
    void add(const VarType& node, const Time& time) {
        data_.push_back(std::make_pair(node, time));
    }

    /**
     * @brief getData Constant getter for data.
     * @return data as vector of coord system and time.
     */
    const std::vector<TimedCoord>& getData() const {
        return data_;
    }

    /**
     * @brief getData Non constant getter for editing topology.
     * @return data as vector of coord system and time.
     */
    std::vector<TimedCoord>& getData() {
        return data_;
    }

    /**
     * @brief print Print topology to cout.
     * @todo turn that into operator<<
     */
    void print() const {
        for (const auto& el : getData()) {
            std::cout << "node " << Streamer{el.first} << " at " << el.second.count() << " ns; ";
        }
        std::cout << std::endl;
    }

private:
    std::vector<TimedCoord> data_; ///> Topology data with coordinate systems and time.
};
} // namespace ct
