// Copyright 2019. All rights reserved.
// Institute of Measurement and Control Systems
// Karlsruhe Institute of Technology, Germany
//
// authors:
//  Johannes Graeter (johannes.graeter@posteo.de)
//  Johannes Beck (johannes.beck@kit.edu)

#pragma once
#include <exception>
#include <iostream>
#include <type_traits>
#include <Eigen/Dense>
#include <internal/coordinate_system_storage.hpp>
#include <internal/tf_traits.hpp>

namespace ct {

// ///////////////// Template declarations.
template <typename ToCoord, typename FromCoord, typename DataType>
class Transform;

template <typename C, typename DataType>
class Point;

// ///////////////// Namespace for all non-trait conversions.
namespace conversion {
/**
 * @brief toDynamic, wraps up traits to cast a transform (static, partially static or dynamic) to dynamic.
 * @todo Auto only returns by value. If toDynamic is called for already purely dynamic transforms, the copy is
 * unnecessary.
 * @param t, transform to cast (static, partially static or dynamic).
 * @return casted dynamic transform
 */
template <typename To, typename From, typename DataType>
static auto toDynamic(const Transform<To, From, DataType>& t) {
    auto coordTo = TfTrait<To>::toDynamic(t.getToCoord());
    auto coordFrom = TfTrait<From>::toDynamic(t.getFromCoord());
    return Transform<decltype(coordTo), decltype(coordFrom), DataType>(t.getData(), coordTo, coordFrom);
}

/**
 * @brief toDynamic, wraps up traits to cast a Point (static dynamic) to dynamic.
 * @todo Auto only returns by value. If toDynamic is called for already purely dynamic points, the copy is
 * unnecessary.
 * @param t, Point to cast (static or dynamic).
 * @return casted dynamic Point
 */
template <typename T, typename DataType>
static auto toDynamic(const Point<T, DataType>& p) {
    auto coord = TfTrait<T>::toDynamic(p.getCoord());
    return Point<decltype(coord), DataType>(p.getData(), coord);
}
} // namespace conversion

template <typename C, typename DataType>
class Point : public CoordinateSystemStorage<TfTrait<C>::IsStatic, CoordinateTag::Other, C> {
public:
    // Specify Eigen Alignment in case DataType is Eigen type.
    // This should be obsolete with c++17.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    explicit Point(DataType data) : data_(std::move(data)) {
        static_assert(TfTrait<C>::IsStatic,
                      "Error in save_coordinate_transform: Point has no coordinate system given, but is dynamic.");
    }

    Point(DataType data, const C& coordinateSystem)
            : CoordinateSystemStorage<TfTrait<C>::IsStatic, CoordinateTag::Other, C>(coordinateSystem),
              data_(std::move(data)) {
    }

    const DataType& getData() const {
        return data_;
    }

    /**
     * @brief getCoord, Getter for coordinate system.
     * Cannot be const C& since coordinate system is a temporary object if coordinate system is static.
     */
    constexpr C getCoord() const {
        return CoordinateSystemStorage<TfTrait<C>::IsStatic, CoordinateTag::Other, C>::get();
    }


private:
    DataType data_; ///> Data of point.
};

/**
 * @brief Typedef for a set of points, operator* of transforms is specialized so that format is compatible to
 * Eigen::Isometry. Size can also be Eigen::Dynamic.
 */
template <typename C, typename ScalarType, int Size>
using Points = Point<C, Eigen::Matrix<ScalarType, 3, Size>>;

/**
 * @brief The TransformException struct
 * Excpetion thrown if coordinate systems do not match.
 */
struct TransformException : public std::exception {
    const char* what() const noexcept override {
        return "Error in coordinate_save: Coordinate systems do not match!";
    }
};


/**
 * @brief The Transform class
 * Use this class to assign coordinate systems to trasnsformations.
 * Depending on the type of coordinate system, the correct concatenation of transforms is checked at compile- or
 * runtime.
 * The type of the transform data (DataType) is flexible, however it must implemented operator* since it is used for
 * concatenation.
 * For usage examples look at unittests.
 */
template <typename ToCoord, typename FromCoord, typename DataType>
class Transform : public CoordinateSystemStorage<TfTrait<ToCoord>::IsStatic, CoordinateTag::To, ToCoord>,
                  CoordinateSystemStorage<TfTrait<FromCoord>::IsStatic, CoordinateTag::From, FromCoord> {
public:
    // Specify Eigen Alignment in case DataType is Eigen type.
    // This should be obsolete with c++17.
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
public:
    /**
     * @brief Transform
     * @param data transform data.
     */
    explicit Transform(const DataType& data) : data_(data) {
        static_assert(TfTrait<FromCoord>::IsStatic && TfTrait<ToCoord>::IsStatic,
                      "Error in save_coordinate_transform: Transform has no coordinate system given, but has dynamic "
                      "coordinate systems.");
    }

    /**
     * @brief Transform
     * @param data transform data.
     * @param toCoord coordinate system to transform to.
     * @param fromCoord coordinate system to trasnform from.
     */
    Transform(const DataType& data, const ToCoord& toCoord, const FromCoord& fromCoord)
            : CoordinateSystemStorage<TfTrait<ToCoord>::IsStatic, CoordinateTag::To, ToCoord>(toCoord),
              CoordinateSystemStorage<TfTrait<FromCoord>::IsStatic, CoordinateTag::From, FromCoord>(fromCoord),
              data_(data) {
    }

    /**
     * @brief operator* Concatenate two transforms.
     * @param otherT Transform to concatenate.
     * @return Concatenated transform.
     */
    template <typename OtherTo, typename OtherFrom>
    Transform<ToCoord, OtherFrom, DataType> operator*(const Transform<OtherTo, OtherFrom, DataType>& otherT) const {
        // Test at compile time of coordinate systems are correct.
        // For dynamic systems of the same type this is always true.
        if constexpr (TfTrait<FromCoord>::IsStatic && TfTrait<OtherTo>::IsStatic) {
            static_assert(TfTrait<FromCoord>::isSame(otherT.getToCoord()));
        } else {
            // Check at runtime if coordinate systems are correct.
            if (!TfTrait<FromCoord>::isSame(getFromCoord(), otherT.getToCoord())) {
                throw TransformException();
            }
        }
        return Transform<ToCoord, OtherFrom, DataType>(
            this->getData() * otherT.getData(), getToCoord(), otherT.getFromCoord());
    }

    /**
     * @brief operator*= Concatenation for transforms with same dynamic type or static convertible type.
     * @param rhs Transform to concatenate.
     * @return Concatenated transform.
     */
    template <typename OtherTo, typename OtherFrom>
    Transform<ToCoord, FromCoord, DataType>& operator*=(const Transform<OtherTo, OtherFrom, DataType>& rhs) {
        // This only works for dynamic types -> convert.
        auto rhsDyn = conversion::toDynamic(rhs);

        assert(TfTrait<ToCoord>::isSame(getToCoord(), rhsDyn.getFromCoord()));
        static_assert(!TfTrait<ToCoord>::IsStatic, "In Transform: Concatenation only valid if lhs is dynamic.");
        static_assert(!TfTrait<FromCoord>::IsStatic, "In Transform: Concatenation only valid if lhs is dynamic.");
        static_assert(
            std::is_same<decltype(rhsDyn.getFromCoord()), FromCoord>::value,
            "In Transform: Concatenation only if FromCoord of rhs has same or castable type to FromCoord of lhs.");

        *this = rhsDyn * (*this);
        return *this;
    }

    /**
     * @brief operator* Transform Point to ToCoord.
     * @param point Point to transform. PointDataType must be compatible with DataType.
     * @return Transformed point.
     */
    template <typename OtherFrom, typename PointDataType>
    Point<ToCoord, PointDataType> operator*(const Point<OtherFrom, PointDataType>& point) const {
        // Test at compile time of coordinate systems are correct.
        // For dynamic systems of the same type this is always true.
        if constexpr (TfTrait<FromCoord>::IsStatic && TfTrait<OtherFrom>::IsStatic) {
            static_assert(TfTrait<FromCoord>::isSame(point.getCoord()));
        } else {
            // Test dynamic system.
            if (!TfTrait<FromCoord>::isSame(getFromCoord(), point.getCoord())) {
                throw TransformException();
            }
        }
        return Point<ToCoord, PointDataType>(this->getData() * point.getData(), getToCoord());
    }

    /**
     * @brief operator* Overload for points to bring them in correct format for multipliying.
     * @param points Set of points to transform. Internally represented as Eigen::Matrix. Size can be Eigen::Dynamic.
     * @return Transformed points.
     */
    template <typename OtherFrom, typename ScalarType, int Size>
    Points<ToCoord, ScalarType, Size> operator*(const Points<OtherFrom, ScalarType, Size>& points) const {
        // Test at compile time of coordinate systems are correct.
        // For dynamic systems of the same type this is always true.
        if constexpr (TfTrait<FromCoord>::IsStatic && TfTrait<OtherFrom>::IsStatic) {
            static_assert(TfTrait<FromCoord>::isSame(points.getCoord()));
        } else {
            // Test dynamic system.
            if (!TfTrait<FromCoord>::isSame(getFromCoord(), points.getCoord())) {
                throw TransformException();
            }
        }
        return Points<ToCoord, ScalarType, Size>(this->getData() * points.getData().colwise().homogeneous(),
                                                 getToCoord());
    }

    /**
     * @brief getData Getter for data.
     * @return Data without coordinate systems.
     */
    const DataType& getData() const {
        return data_;
    }

    /**
     * @brief getToCoord Getter for coordinate system.
     * Cannot be const C& since coordinate system is a temporary object if coordinate system is static.
     * @return coordinate system
     */
    constexpr ToCoord getToCoord() const {
        return CoordinateSystemStorage<TfTrait<ToCoord>::IsStatic, CoordinateTag::To, ToCoord>::get();
    }

    /**
     * @brief getFromCoord Getter for coordinate system.
     * Cannot be const C& since coordinate system is a temporary object if coordinate system is static.
     * @return coordinate system
     */
    constexpr FromCoord getFromCoord() const {
        return CoordinateSystemStorage<TfTrait<FromCoord>::IsStatic, CoordinateTag::From, FromCoord>::get();
    }

private:
    DataType data_; ///> Data of transform.
};
} // namespace ct
