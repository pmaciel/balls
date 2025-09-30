add_library(balls
    BallTree.cc
    BallTree.h
    Metric.cc
    Metric.h
    Vec3.cc
    Vec3.h
)

target_include_directories(balls PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
