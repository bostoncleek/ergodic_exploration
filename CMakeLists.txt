cmake_minimum_required(VERSION 3.9)
project(ergodic_exploration)

# add compile options
add_compile_options(-Wall -Wextra)

# Compile as C++20
set(CMAKE_CXX_STANDARD 20)
# don't use gnu extensions
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
	geometry_msgs
	nav_msgs
  roscpp
	rospy
)

# find_package(Boost REQUIRED COMPONENTS system)
find_package(Armadillo REQUIRED)


## Generate messages in the 'msg' folder
# add_message_files(
#   FILES
#   Message1.msg
#   Message2.msg
# )

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
# generate_messages(
#   DEPENDENCIES
#   std_msgs  # Or other packages containing msgs
# )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS
	roscpp
	rospy
 # DEPENDS
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
include
  ${catkin_INCLUDE_DIRS}
	${ARMADILLO_INCLUDE_DIRS}
)

# Declare a C++ library
add_library(${PROJECT_NAME}
	include/${PROJECT_NAME}/cart.hpp
	include/${PROJECT_NAME}/integrator.hpp
	include/${PROJECT_NAME}/omni.hpp
	src/${PROJECT_NAME}/collision.cpp
	src/${PROJECT_NAME}/grid.cpp
	src/${PROJECT_NAME}/numerics.cpp
	# src/${PROJECT_NAME}/ergodic_exploration.cpp
)

## Add cmake target dependencies of the library
# add_dependencies(${PROJECT_NAME}
# 	${${PROJECT_NAME}_EXPORTED_TARGETS}
# 	${catkin_EXPORTED_TARGETS}
# )


# Declare a C++ executable
add_executable(exploration src/exploration_node.cpp)

## Rename C++ executable without prefix
## The above recommended prefix causes long target names, the following renames the
## target back to the shorter version for ease of user use
## e.g. "rosrun someones_pkg node" instead of "rosrun someones_pkg someones_pkg_node"
# set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME node PREFIX "")

## Add cmake target dependencies of the executable
# add_dependencies(exploration ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

# Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
	${ARMADILLO_LIBRARIES}
)

target_link_libraries(exploration
  ${catkin_LIBRARIES}
	${PROJECT_NAME}
	${ARMADILLO_LIBRARIES}
)

#############
## Install ##
#############

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
catkin_install_python(PROGRAMS
  scripts/grid_saver
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

# Mark executables for installation
install(TARGETS exploration
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

# Mark libraries for installation
install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)

# Mark cpp header files for installation
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

## Mark other files for installation (e.g. launch and bag files, etc.)
install(DIRECTORY config/ DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
install(DIRECTORY launch/ DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_ergodic_exploration.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

if(CATKIN_ENABLE_TESTING)
    catkin_add_gtest(${PROJECT_NAME}_test
			test/test_cart.cpp
			test/test_grid.cpp
			test/test_integrator.cpp
			test/test_omni.cpp
			)
    target_link_libraries(${PROJECT_NAME}_test
			${catkin_Libraries}
			${ARMADILLO_LIBRARIES}
			${PROJECT_NAME}
 			gtest_main
			)
endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
