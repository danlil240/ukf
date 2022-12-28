cmake_minimum_required(VERSION 3.24)
project(los_rates_estimator)


find_package(Boost COMPONENTS thread chrono REQUIRED)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_srvs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(ament_index_cpp REQUIRED)
find_package(ros_common REQUIRED)
find_package(ros_common_interfaces REQUIRED)
find_package(gimbal_controller_interfaces REQUIRED)
find_package(eyeit_interfaces REQUIRED)
find_package(target_estimator_interfaces REQUIRED)
find_package(fcu_driver_interfaces REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(diagnostic_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp)


find_package(PythonLibs REQUIRED)


find_package(ament_cmake_gtest REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
include(${ros_common_INCLUDE_DIRS}/cmake/GetVersionInformation.cmake)
get_git_version(${CMAKE_SOURCE_DIR})
list(INSERT CMAKE_MODULE_PATH 0 "${ros_common_INCLUDE_DIRS}/cmake")
find_package(GeographicLib REQUIRED)


# ---------------------------------------------------------------------

# Specify additional locations of header files
include_directories(
        ${ros_common_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
)

add_executable(los_rates_estimator main.cpp UKF.cpp model.cpp test_ukf.cpp)
add_executable(ukf UKF.cpp model.cpp)

target_link_libraries(ukf
        ${PYTHON_LIBRARIES}
        )


add_executable(test_ukf test_ukf.cpp model.cpp  UKF.cpp)

target_link_libraries(test_ukf
        ${PYTHON_LIBRARIES}
        yaml-cpp
        )

target_include_directories(test_ukf PRIVATE ${PYTHON_INCLUDE_DIRS})
