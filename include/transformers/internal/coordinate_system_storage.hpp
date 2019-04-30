// Copyright 2019. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@posteo.de)
//  Johannes Beck (johannes.beck@kit.edu)

#pragma once

/**
 * @brief The CoordinateTag enum
 */
enum class CoordinateTag { From, To, Other };

/**
 * @brief The CoordinateSystemStorage struct
 * This class hides the coordinate system management from the Transform and Point class.
 * If the coord system is dynamic its value is saved on this class, otherwise it is optimizied by empty base class
 * optimization.
 * @tparam IsStatic is it a compiletime static coordinate system?
 * @tparam Tag defines if it is a From, To or Other coordinate system in order to be able to distinguish them.
 * @tparam Coord coordinate system type.
 */
template <bool IsStatic, CoordinateTag Tag, typename Coord>
struct CoordinateSystemStorage;

/**
 * @brief The CoordinateSystemStorage<false, Tag, Coord> struct
 * Specialization of a CoordinateSystemStorage for dynamic systems.
 */
template <CoordinateTag Tag, typename Coord>
struct CoordinateSystemStorage<false, Tag, Coord> {
public:
    explicit CoordinateSystemStorage(Coord c) : coord_(std::move(c)) {
    }
    const Coord& get() const {
        return coord_;
    }

private:
    Coord coord_;
};


/**
 * @brief The CoordinateSystemStorage<true, Tag, Coord> struct
 * Specialization of a CoordinateSystemStorage for static systems.
 */
template <CoordinateTag Tag, typename Coord>
struct CoordinateSystemStorage<true, Tag, Coord> {

    /**
     * @brief CoordinateSystemStorage
     * Default constructor, used purely static transform construction.
     */
    CoordinateSystemStorage() {
    }

    /**
     * @brief CoordinateSystemStorage
     * Same interface for constructor is needed, than for dynamic systems.
     */
    explicit CoordinateSystemStorage(const Coord& /*c*/) {
    }
    Coord get() const {
        return Coord();
    }
};
