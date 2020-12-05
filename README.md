# Ergodic Exploration

- [Motivation](#Motivation) </br>
- [Getting Started](#Getting-Started) </br>
  - [Dependencies](#Dependencies) </br>
  - [Workspace](#Workspace) </br>
  - [Launch](#Launch) </br>
- [Exploration](#Exploration) </br>
  - [Ergodic Control](#Ergodic-Control) </br>
  - [Exploration Stack](#Exploration-Stack) </br>
  - [Collision Detection](#Collision-Detection) </br>
- [Results](#Results) </br>
- [Nodes](#Nodes) </br>
  - [Subscribed Topics](#Subscribed-Topics) </br>
  - [Published Topics](#Published-Topics) </br>
  - [Parameters](#Parameters) </br>
  - [Required tf Transforms](#Required-tf-Transforms) </br>
- [Issues](#Credits) </br>
- [Credits](#Credits) </br>

# Motivation
Currently there are no exploration packages for ROS Noetic. The ROS packages that currently exist perform only frontier exploration. The goal of this project is to provide a robot agnostic information theoretic exploration strategy for known and unknown environments.

This package requires a user specified target distribution representing the expected
information gain. The target distribution is represented as a Gaussian or multiple Gaussians.

As an alternative mutual information can be used as the target distribution. Checkout the `feature/sportdeath_mi` branch. That branch is experimental because the mutual information implementation has not been verified yet.

# Getting Started
## Dependencies
This package was built and tested using [Armadillo](http://arma.sourceforge.net/) 10.1.1. You will need to install Armadillo.

## Workspace
Create a catkin workspace or clone in an existing workspace

```
mkdir -p ~/exploration_ws/src
cd ~/exploration_ws/src
git clone https://github.com/bostoncleek/ergodic_exploration.git
```
```
cd ~/exploration_ws
catkin init
catkin build -DCMAKE_BUID_TYPE=Release
source devel/setup.bash
```

## Launch
Example exploration configurations for an [omni directional](https://github.com/bostoncleek/ergodic_exploration/blob/master/config/explore_omni.yaml) and a [differential
drive](https://github.com/bostoncleek/ergodic_exploration/blob/master/config/explore_cart.yaml) robot are provided.

Spawn either a omni directional or a differential drive robot in Gazebo and start your
favorite SLAM package. Use the [example launch](https://github.com/bostoncleek/ergodic_exploration/blob/master/launch/exploration.launch) file to start the exploration node. Be sure to set the `holonomic` arg to false if using a differential drive robot.

```
roslaunch ergodic_exploration exploration.launch
```

# Exploration
## Ergodic Control

## Exploration Stack

## Collision Detection


# Results


# Nodes
The nodes provided are `explore_cart` and `explore_omni`. The only difference between them
is the motion model. Both use an occupancy grid for collision detection and provide a body twist as the control output.

## Subscribed Topics
- map ([nav_msgs/OccupancyGrid](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/OccupancyGrid.html)): required for collision detection
- odom ([nav_msgs/Odometry](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/Odometry.html)): required for control feedback

## Published Topics
- cmd_vel ([geometry_msgs/Twist](http://docs.ros.org/en/noetic/api/geometry_msgs/html/msg/Twist.html)): body twist control update
- trajectory ([nav_msgs/Path](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/Path.html)): ergodic controller optimized trajectory
- dwa_trajectory ([nav_msgs/Path](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/Path.html)): dynamic window trajectory
- target ([visualization_msgs/MarkerArray](http://docs.ros.org/en/noetic/api/visualization_msgs/html/msg/MarkerArray.html)): target distribution

## Parameters

The following parameters are set to 0 for the cart because of the holonomic constraint:
max_vel_y, min_vel_y, acc_lim_y, and vy_samples. You only need to proved the control weights for the x-velocity component and the rotational velocity.

### General
- map_frame_id (string, default: "map"): map frame id
- base_frame_id (string, default: "base_link"): base link frame id
- frequency (double, default: 10.0): control loop frequency (Hz)
- val_dt (double, default: 0.1): control validation time step for collision detection (s)
- val_horizon (double, default: 0.5): control validation horizon for collision detection (s)
- max_vel_x (double, default: 1.0): max x velcocity (m/s)
- max_vel_y (double, default: 1.0):  max y velcocity (m/s)
- max_rot_vel (double, default: 1.0): max ratotional velcocity (m/s)
- min_vel_x (double, default: -1.0): min x velcocity (m/s)
- min_vel_y (double, default: -1.0): min y velcocity (m/s)
- min_rot_vel (double, default: -1.0): min ratotional velcocity (m/s)
- acc_lim_x (double, default: 1.0): x acceleration limit (m/s^2)
- acc_lim_y (double, default: 1.0): y acceleration limit (m/s^2)
- acc_lim_th (double, default: 1.0): rotational acceleration limit (m/s^2)

### Collision Parameters
- boundary_radius (double, default: 0.7): bounding radius around robot (m)
- search_radius (double, default: 1.0): max search radius for collision detection (m)
- obstacle_threshold (double, default: 0.2): obstacles within radius from boundary are cosidered collisions (m)
- occupied_threshold (double, default: 0.8): occupancy grid cell probability to be considered an obstacle [0 1]

### Ergodic Control Parameters
- ec_dt (double, default: 0.1): time step used in ergodic control integration (s)
- ec_horizon (double, default: 2.0): ergodic control horizon (s)
- target_resolution (double, default: 0.1): target grid resolution (m)
- expl_weight (double, default: 1.0): ergodic exploration weight influences how ergodic the exploration is
- num_basis (unsigned int, default: 10): number of basis functions used in composes the trajectory and spatial fourier coefficients
- buffer_size (unsigned int, default: 1e6): total number of past states stored in memory
- batch_size (unsigned int, default: 100): number of past states randomly sampled in each ergodic control loop and is used to compose the trajectory fourier coefficients
- control_weights (std::vector<double>, default: [1, 1, 1]): weights on twist [vx vy w] used by the ergodic controller

### Dynamic Window Parameters
- dwa_dt (double, default: 0.1): time step used in dynamic window integration (s)
- dwa_horizon (double, default: 1.0): dynamic window control horizon (s)
- acc_dt (double, default: 0.2): time step the acceleration limits are applied (s)
- vx_samples (unsigned int, default: 3): number of x velcocity samples
- vy_samples (unsigned int, default: 8): number of y velcocity samples
- vth_samples (unsigned int, default: 5): number of rotational velcocity samples
- means (double array[], default: none): target x and y means (m)
- sigmas (double array[], default: none): target x and y standard deviations (m)


## Required tf Transforms
- map -> base_link : usually provided in combination by odometry and SLAM systems

# Issues


# Credits
