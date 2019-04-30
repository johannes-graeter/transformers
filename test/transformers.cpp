// google test docs
// wiki page: https://code.google.com/p/googletest/w/list
// primer: https://code.google.com/p/googletest/wiki/V1_7_Primer
// FAQ: https://code.google.com/p/googletest/wiki/FAQ
// advanced guide: https://code.google.com/p/googletest/wiki/V1_7_AdvancedGuide
// samples: https://code.google.com/p/googletest/wiki/V1_7_Samples
//
// List of some basic tests fuctions:
// Fatal assertion                      Nonfatal assertion                   Verifies / Description
//-------------------------------------------------------------------------------------------------------------------------------------------------------
// ASSERT_EQ(expected, actual);         EXPECT_EQ(expected, actual);         expected == actual
// ASSERT_NE(val1, val2);               EXPECT_NE(val1, val2);               val1 != val2
// ASSERT_LT(val1, val2);               EXPECT_LT(val1, val2);               val1 < val2
// ASSERT_LE(val1, val2);               EXPECT_LE(val1, val2);               val1 <= val2
// ASSERT_GT(val1, val2);               EXPECT_GT(val1, val2);               val1 > val2
// ASSERT_GE(val1, val2);               EXPECT_GE(val1, val2);               val1 >= val2
//
// ASSERT_FLOAT_EQ(expected, actual);   EXPECT_FLOAT_EQ(expected, actual);   the two float values are almost equal (4
// ULPs) ASSERT_DOUBLE_EQ(expected, actual);  EXPECT_DOUBLE_EQ(expected, actual);  the two double values are almost
// equal
// (4 ULPs) ASSERT_NEAR(val1, val2, abs_error);  EXPECT_NEAR(val1, val2, abs_error);  the difference between val1 and
// val2 doesn't exceed the given absolute error
//
// Note: more information about ULPs can be found here:
// http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
//
// Example of two unit test:
// TEST(Math, Add) {
//    ASSERT_EQ(10, 5+ 5);
//}
//
// TEST(Math, Float) {
//	  ASSERT_FLOAT_EQ((10.0f + 2.0f) * 3.0f, 10.0f * 3.0f + 2.0f * 3.0f)
//}
//=======================================================================================================================================================
#include <save_coordinate_transform.hpp>
#include <topology_graph.hpp>
#include "gtest/gtest.h"

// A google test function (uncomment the next function, add code and
// change the names TestGroupName and TestName)
TEST(TfTraits, all) {
    DEFINE_STATIC_COORDINATE_SYSTEM(C0);
    Eigen::Vector3d p0 = Eigen::Vector3d(1., 0., 0.);
    ct::Point<C0, decltype(p0)> b0(p0);
    auto hash = TfTrait<decltype(b0.getCoord())>::toDynamic(b0.getCoord());

    ct::Point<decltype(hash), decltype(p0)> b0Dyn = ct::conversion::toDynamic(b0);

    ASSERT_EQ(hash, typeid(C0).hash_code());
    ASSERT_EQ(hash, b0Dyn.getCoord());
}

TEST(SaveCoordinateTransformations, concatenation) {

    ///@todo Split tests.

    using Ttd = Eigen::Isometry3d;
    using Ptd = Eigen::Vector3d;
    Ttd a21 = Ttd::Identity();
    a21.translate(Eigen::Vector3d(1., 1.4, 3.));
    a21.rotate(Eigen::AngleAxisd(0.5, Eigen::Vector3d(0., 0., 1.)).toRotationMatrix());
    Ttd a10 = Ttd::Identity();
    a10.translate(Eigen::Vector3d(1., 0., 0.));
    Ttd a30 = a21 * a10;

    Ptd p0 = Eigen::Vector3d(1., 0., 0.);

    Ttd res2;
    res2 = a21 * a10;
    // Static coordinate system system.
    {
        // This should be identical to macro.
        struct C0 : public StaticTf {};

        DEFINE_STATIC_COORDINATE_SYSTEM(C1);
        DEFINE_STATIC_COORDINATE_SYSTEM(C2);

        ASSERT_TRUE(TfTrait<C0>::IsStatic);
        ASSERT_TRUE(TfTrait<C1>::IsStatic);

        ct::Transform<C2, C1, Ttd> trans21(a21);
        ct::Transform<C1, C0, Ttd> trans10(a10);

        ct::Point<C0, Ptd> b0(p0);

        ct::Transform<C2, C0, Ttd> res = trans21 * trans10;
        ct::Point<C1, Ptd> b1 = trans10 * b0;
        auto b0Dyn = ct::conversion::toDynamic(b0);
        ct::Point<C1, Ptd> b2 = trans10 * b0Dyn;
        auto trans10Dyn = ct::conversion::toDynamic(trans10);
        auto b3 = trans10Dyn * b0Dyn;
        auto b4 = trans10Dyn * b0;

        std::cout << res.getData().matrix() << std::endl;

        ASSERT_TRUE(res2.isApprox(res.getData()));
        ASSERT_NEAR(((a10 * p0).eval() - b1.getData()).norm(), 0., 1e-15);
        ASSERT_NEAR((b3.getData() - b2.getData()).norm(), 0., 1e-15);
        ASSERT_NEAR((b4.getData() - b2.getData()).norm(), 0., 1e-15);
    }
    // Dynamic coordinate system.
    {
        ASSERT_FALSE(TfTrait<std::string>::IsStatic);
        ASSERT_FALSE(TfTrait<int>::IsStatic);

        ct::Transform<int, std::string, Ttd> trans21(a21, 1, "c1");
        ct::Transform<std::string, int, Ttd> trans10(a10, "c1", 0);

        ASSERT_EQ(trans21.getToCoord(), 1);
        ASSERT_EQ(trans21.getFromCoord(), "c1");

        ct::Point<int, Ptd> b0(p0, 0);

        ct::Transform<int, int, Ttd> res = trans21 * trans10;
        ct::Point<std::string, Ptd> b1 = trans10 * b0;


        std::cout << res.getData().matrix() << std::endl;

        ASSERT_TRUE(res2.isApprox(res.getData()));
        ASSERT_NEAR(((a10 * p0).eval() - b1.getData()).norm(), 0., 1e-15);
    }

    // Dynamic and static coordinate systems.
    {

        DEFINE_STATIC_COORDINATE_SYSTEM(C0);

        ct::Transform<int, std::string, Ttd> trans21(a21, 1, "c1");
        ct::Transform<std::string, C0, Ttd> trans10(a10, "c1", C0{});

        ct::Point<C0, Ptd> b0(p0);

        ct::Transform<int, C0, Ttd> res = trans21 * trans10;
        ct::Point<std::string, Ptd> b1 = trans10 * b0;

        std::cout << res.getData().matrix() << std::endl;

        ASSERT_TRUE(res2.isApprox(res.getData()));
        ASSERT_NEAR(((a10 * p0).eval() - b1.getData()).norm(), 0., 1e-15);
    }

    {
        DEFINE_STATIC_COORDINATE_SYSTEM(C0);
        DEFINE_STATIC_COORDINATE_SYSTEM(C1);

        ct::Transform<C1, C0, Ttd> trans21(a21);
        ct::Transform<C0, std::string, Ttd> trans10(a10, C0{}, "base");

        ct::Point<std::string, Ptd> b0(p0, "base");

        ct::Transform<C1, std::string, Ttd> res = trans21 * trans10;
        ct::Point<C0, Ptd> b1 = trans10 * b0;

        std::cout << res.getData().matrix() << std::endl;

        ASSERT_TRUE(res2.isApprox(res.getData()));
        ASSERT_NEAR(((a10 * p0).eval() - b1.getData()).norm(), 0., 1e-15);
    }
    {
        DEFINE_STATIC_COORDINATE_SYSTEM(C0);
        DEFINE_STATIC_COORDINATE_SYSTEM(C1);

        ct::Transform<C1, std::string, Ttd> trans21(a21, C1{}, "c2");
        ct::Transform<std::string, int, Ttd> trans10(a10, "c1", 0);

        bool didThrow = false;
        try {
            ct::Transform<C1, int, Ttd> res = trans21 * trans10;
            std::cout << res.getData().matrix() << std::endl;
        } catch (const ct::TransformException& e) {
            didThrow = true;
            std::cout << "Exception thrown as wanted! message: " << e.what() << std::endl;
        }
        ASSERT_TRUE(didThrow);
    }

    // Test concatenation.
    {
        DEFINE_STATIC_COORDINATE_SYSTEM(C0);
        DEFINE_STATIC_COORDINATE_SYSTEM(C1);
        DEFINE_STATIC_COORDINATE_SYSTEM(C2);
        DEFINE_STATIC_COORDINATE_SYSTEM(C3);

        ct::Transform<C1, C0, Ttd> trans10(a10);
        ct::Transform<C2, C1, Ttd> trans21(a21);
        ct::Transform<C3, C2, Ttd> trans30(a30);

        auto transAccDyn = ct::conversion::toDynamic(trans10);
        transAccDyn *= trans21;
        transAccDyn *= trans30;

        ASSERT_TRUE(TfTrait<C3>::isSame(transAccDyn.getToCoord()));
        ASSERT_TRUE(TfTrait<C0>::isSame(transAccDyn.getFromCoord()));
        ASSERT_TRUE(transAccDyn.getData().isApprox(a30 * a21 * a10));
    }
}

TEST(SaveCoordinateTransformations, PointSets) {

    DEFINE_STATIC_COORDINATE_SYSTEM(C0);
    DEFINE_STATIC_COORDINATE_SYSTEM(C1);

    ///@todo Split tests.
    using Ttd = Eigen::Isometry3d;
    using Ptd = Eigen::Vector3d;

    Ttd a10 = Ttd::Identity();
    a10.translate(Eigen::Vector3d(1., 1.4, 3.));
    a10.rotate(Eigen::AngleAxisd(0.5, Eigen::Vector3d(0., 0., 1.)).toRotationMatrix());

    Ptd p0 = Eigen::Vector3d(1., 0., 0.);

    // Dynamic size.
    {
        ct::Transform<C1, C0, Ttd> t(a10);
        using Pointcloud = Eigen::Matrix3Xd;
        Pointcloud pcl(3, 4);
        pcl << p0, (p0 + Eigen::Vector3d(0.2, 1, 5.)), (p0 + Eigen::Vector3d(0.89912, 1.24, 2.)),
            (p0 + Eigen::Vector3d(-124.1, -2.41, 0.23));

        ct::Points<C0, double, Eigen::Dynamic> pts(pcl);
        ct::Points<C1, double, Eigen::Dynamic> pts2 = t * pts;

        Pointcloud transformed = a10 * pcl.colwise().homogeneous();
        ASSERT_TRUE(transformed.isApprox(pts2.getData()));
    }

    // Static size.
    {
        ct::Transform<std::string, int, Ttd> t(a10, "end", 0);
        using Pointcloud = Eigen::Matrix<double, 3, 4>;
        Pointcloud pcl;
        pcl << p0, (p0 + Eigen::Vector3d(0.2, 1, 5.)), (p0 + Eigen::Vector3d(0.89912, 1.24, 2.)),
            (p0 + Eigen::Vector3d(-124.1, -2.41, 0.23));

        ct::Points<int, double, 4> pts(pcl, 0);
        ct::Points<std::string, double, 4> pts2 = t * pts;

        ASSERT_EQ(pts2.getCoord(), "end");

        Pointcloud transformed = a10 * pcl.colwise().homogeneous();
        ASSERT_TRUE(transformed.isApprox(pts2.getData()));
    }
}

TEST(TopologyGraph, get) {

    ct::TopologyGraph<int, std::string> graph;

    //    DEFINE_STATIC_COORDINATE_SYSTEM(C0);
    //    DEFINE_STATIC_COORDINATE_SYSTEM(C1);
    //    DEFINE_STATIC_COORDINATE_SYSTEM(C2);
    //    DEFINE_STATIC_COORDINATE_SYSTEM(C3);

    graph.addEdge(0, 1);
    graph.addEdge(1, 7);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(1, 5);
    graph.addEdge(5, 6);
    ct::TopologyPath<int, std::string> top0 = graph.get(3, 0, 123345);
    top0.print();

    graph.addEdge("s0", "s1");
    graph.addEdge("s0", "s7");
    graph.addEdge("s1", "s2");
    graph.addEdge("s2", "s3");
    graph.addEdge("s3", "s4");
    graph.addEdge("s1", "s5");
    graph.addEdge("s5", "s6");
    graph.print();

    ct::TopologyPath<int, std::string> top2 = graph.get(0, 3, 123345);
    {
        ct::TopologyPath<int, std::string> top1 = graph.get("s3", "s0", 123345);
        top1.print();
        top2.print();
        ASSERT_EQ(top0.getData().size(), 4);
    }
    {
        graph.addEdge(0, "s0");
        ct::TopologyPath<int, std::string> top3 = graph.get("s3", 0, 123345);
        auto topConcat = top3 * top2;
        topConcat.print();
        ASSERT_EQ(topConcat.getData().size(), 8);
    }
    ct::TopologyPath<int, std::string> top4 = graph.get(6, 0, 123345);
    ct::TopologyPath<int, std::string> top5 = graph.get(0, 7, 123345);
    auto topConcat2 = top4 * top5;
    std::cout << "concat normal:" << std::endl;
    top4.print();
    top5.print();
    topConcat2.print();
    ASSERT_EQ(topConcat2.getData().size(), 4);
    {
        ct::TopologyPath<int, std::string> top4 = graph.get(6, 0, 123345);
        ct::TopologyPath<int, std::string> top5 = graph.get(0, 7, 223345);
        auto topConcat2 = top4 * top5;
        std::cout << "concat time diff:" << std::endl;
        top4.print();
        top5.print();
        topConcat2.print();
        ASSERT_EQ(topConcat2.getData().size(), 7);
    }

    {
        bool wasThrown = false;
        ct::TopologyPath<int, std::string> top4 = graph.get(6, 0, 123345);
        ct::TopologyPath<int, std::string> top5 = graph.get(7, 0, 123345);
        try {
            auto unused = top4 * top5;
        } catch (const ct::NotConnectedException& e) {
            std::cout << "Exception thrown as expected! message: " << e.what() << std::endl;
            wasThrown = true;
        }
        ASSERT_TRUE(wasThrown);

        auto topConcat3 = top4 * top5.inverse();
        top4.print();
        top5.print();
        topConcat3.print();
        ASSERT_EQ(topConcat2.getData().size(), 4);
        ASSERT_TRUE(topConcat2 == topConcat3);
    }
}

// TEST(Transformer, basic){

//    DEFINE_STATIC_COORDINATE_SYSTEM(C0);
//    DEFINE_STATIC_COORDINATE_SYSTEM(C2);
//    DEFINE_STATIC_COORDINATE_SYSTEM(C3);

//    using TT=Eigen::Isometry3d;

//    TT trans0, trans1, trans2;
//    trans0.translate(Eigen::Vector3d(1.,0.,1.));
//    trans0.rotate(Eigen::AngleAxisd(0.3, Eigen::Vector3d(1.,0.,1.)));
//    trans1=trans0;
//    trans1.translate(Eigen::Vector3d(0.5,0.1,1.));
//    trans1.rotate(Eigen::AngleAxisd(0.1, Eigen::Vector3d(1.,1.,1.)));
//    trans2=trans1*trans0;

//    ct::Transform<int,C0, TT> t0(trans0,0,C0{});
//    ct::Transform<C2,int, TT> t1(trans1,0);
//    ct::Transform<C3,C2, TT> t2(trans2);

//    ct::Transformer t;
//    t.add(trans0);
//    t.add(trans2);
//    t.add(trans1);

//    // Test data consistency.
//    ct::Transform<int, C0, TT> c10 = t.get(0, C0);
//    ASSERT_TRUE(c10.getData().isApprox(trans0);
//    ASSERT_TRUE((t.get(C2,C0).isApprox(trans0*trans1);

//    // Test relative transform consistency static.
//    ct::Transform<C2,C0, TT> c20 = t.get(C2, C0);  // From-to
//    ct::Transform<C0,C2, TT> c02 = t.get(C0, C2);  // to-from
//    ASSERT_TRUE(c20.getData().isApprox(c02.getData().inverse()));

//    // Test relative transform consistency dynamic.
//    ct::Transform<C2, int, TT> c21 = t.get(C2, 0);  // Dynamic-static
//    ct::Transform<int, C2, TT> c12 = t.get(0, C2);  // static-dynamic
//    ASSERT_TRUE(c21.getData().isApprox(c12.getData().inverse()));
//}
