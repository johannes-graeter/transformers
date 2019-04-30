// Copyright 2018. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
// authors:
//  Johannes Graeter (johannes.graeter@posteo.de)
//  Johannes Beck (johannes.beck@kit.edu)

#pragma once
#include <string>
#include <type_traits>
#include <internal/tf_traits_impl.hpp>


/**
 * @brief Macro for convenience to define static coordinate system
 *   usage: DEFINE_STATIC_COORDINATE_SYSTEM(C1)
 */
#define DEFINE_STATIC_COORDINATE_SYSTEM(name) \
    struct name : public StaticTf {};

/**
 * @brief General purpose trait. Define isStatic in TfStaticTrait.
 */
template <typename T>
struct TfTrait : public TfTraitImpl<TfStaticTrait<T>::IsStatic, T> {};

// ////////// Traits for default dynamic types.

/**
 * @brief The TfTrait<int> struct
 * Specialization for dynamic coordinate system type int.
 */
template <>
struct TfTrait<int> : public TfTraitImpl<false, int> {};

/**
 * @brief The TfTrait<std::string> struct
 * Specialization for dynamic coordinate system type string.
 */
template <>
struct TfTrait<std::string> : public TfTraitImpl<false, std::string> {};
