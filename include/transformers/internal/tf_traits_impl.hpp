// Copyright 2019. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@posteo.de)
//  Johannes Beck (johannes.beck@kit.edu)

#pragma once
#include <cstddef>
#include <type_traits>
#include <typeinfo>

/**
 * @brief The StaticTf struct
 * Class that static coordinate systems inherit from.
 */
struct StaticTf {
protected:
    // Define default constructor as protected to prevent instantiation.
    StaticTf() = default;
};

using ConvertedDynamicType = decltype(typeid(StaticTf).hash_code());

/**
 * @brief The TfTrait struct
 * Default trait.
 * Coordinate system is static if it inherits from StaticTf.
 */
template <typename T>
struct TfStaticTrait {
    static constexpr bool IsStatic = std::is_base_of<StaticTf, T>::value;
};

/**
 * @brief StaticIsSame<T, OtherT>
 * Implement static is same method in order to account for static coord and converted static coord.
 */
template <bool IsTStatic, bool IsOtherTStatic, typename T, typename OtherT>
struct IsSameImpl;

template <typename T, typename OtherT>
struct IsSameImpl<true, true, T, OtherT> {
    static bool apply(const T& /*a*/, const OtherT& /*b*/) {
        return std::is_same<T, OtherT>::value;
    }
};

/**
 * @brief StaticIsSame<T, DefaultDynamicType>
 * Specialization for converted coord.
 */
template <typename T>
struct IsSameImpl<true, false, T, ConvertedDynamicType> {
    static bool apply(const T& /*a*/, const ConvertedDynamicType& b) {
        return typeid(T).hash_code() == b;
    }
};

template <typename OtherT>
struct IsSameImpl<false, true, ConvertedDynamicType, OtherT> {
    static bool apply(const ConvertedDynamicType& a, const OtherT& /*b*/) {
        return typeid(OtherT).hash_code() == a;
    }
};

/**
 * @brief StaticIsSame<T, DefaultDynamicType>
 * Specialization for dynamic coord.
 */
template <typename T, typename OtherT>
struct IsSameImpl<false, false, T, OtherT> {
    static bool apply(const T& a, const OtherT& b) {
        return a == b;
    }
};

/**
 * @brief The IsSame struct
 *  Test if coordinate systems are the same.
 *  Can  be compile time static coordinate systems (input given by templates) or dynamic (input given as vairbales).
 *  This is used in TfTrait inorder to make interface clearer for special types.
 *  Inputs: Static: templates, dynamic: parameters.
 * @return bool (true if is same)
 * @todo How to add better error warnings if trait is not implemented because it doesn't make sense?
 */
template <typename T, typename OtherT>
struct IsSame : public IsSameImpl<TfStaticTrait<T>::IsStatic, TfStaticTrait<OtherT>::IsStatic, T, OtherT> {};


/**
 * @brief TfTraitImpl<IsStatic, T>
 * Interface declaration of TfTrait.
 * Must implement:
 * static constexpr bool IsStatic
 * static ConvertedDynamicType toDynamic(const T&)
 * static bool isSame(const T& a, const OtherT& b)
 *
 * Could implement:
 * static bool isSame(const OtherT& b), convenience method for static types
 *
 * @tparam IsStatic is true if coordinate systems are compile time static.
 * @tparam T is type of the coordinate system
 *
 */
template <bool IsStatic, typename T>
struct TfTraitImpl;

/**
 * @brief The TfTraitImpl<true, T> struct
 * Specialization for static coordinate systems.
 */
template <typename T>
struct TfTraitImpl<true, T> {
    static constexpr bool IsStatic = true;

    /**
     * @brief toDynamic, convert static coords to dynamic coords, by hashing the typeid,
     * input defined to have the same interface as dynamic coords.
     * @return hash of typieid (should be size_t)
     */
    static ConvertedDynamicType toDynamic(const T& /*a*/) {
        return typeid(T).hash_code();
    }

    /**
     * @brief isSame use StaticIsSame struct to test for coords equality
     * for convenience, interface with unsued input is also defined.
     * @param b coordinate system to compare to.
     * @return bool
     */
    template <typename OtherT>
    static bool isSame(const T& /*a*/, const OtherT& b) {
        return IsSame<T, OtherT>::apply(T{}, b);
    }
    template <typename OtherT>
    static bool isSame(const OtherT& b) {
        return IsSame<T, OtherT>::apply(T{}, b);
    }
};

/**
 * @brief The TfTraitImpl<false, T> struct
 * Specialization for dynamic coordinate systems
 */
template <typename T>
struct TfTraitImpl<false, T> {
    static constexpr bool IsStatic = false;

    /**
     * @brief toDynamic, dummy conversion returning same object, defined for convenience.
     * @param a
     * @return same object.
     */
    static const T& toDynamic(const T& a) {
        return a;
    }

    /**
     * @brief isSame test dynamic coords for equality.
     * @param a first dynamic coordinate system.
     * @param b second dynamic coordinate system.
     * @return bool
     */
    template <typename OtherT>
    static bool isSame(const T& a, const OtherT& b) {
        return IsSame<T, OtherT>::apply(a, b);
    }
};
