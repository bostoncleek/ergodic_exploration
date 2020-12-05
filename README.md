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
- [Future Improvements](#Future-Improvements) </br>
- [Citing Ergodic Exploration](#Citing-Ergodic-Exploration) </br>
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
Ergodicity is defined as the fraction of time spent sampling an area should be equal to a metric quantifying the density information in that area. The ergodic metric is the difference between the probability density func­tions representing the spatial distribution and the statistical
representation of the time-averaged trajectory [1]. The objective function includes the ergodic metric and the control cost.

The ergodic controller performs receding horizon trajectory optimization in real time. The real time performance is achieved by integrating the [co-state](https://en.wikipedia.org/wiki/Hamiltonian_(control_theory)#:~:text=%2C%20referred%20to%20as%20costate%20variables,maximize%20the%20Hamiltonian%2C%20for%20all) backwards in time [2].

## Exploration Stack
The ergodic controller cannot guarantee a collision free trajectory. The dynamic window approach is used as the local planner. The next twist from ergodic controller is propagated forward in for a short time to detect collisions. If there is a collision the dynamic window approach is used. The dynamic window approach forward propagates a constant twist for a short time and disregards it if there is a collisions.

In the case where the ergodic control results in a collision the dynamic window approach uses the optimized trajectory as a reference. The dynamic window approach finds a twist that will produce a trajectory most similar to the reference that is collision free. Similarity between the optimized trajectory and the dynamic window approach trajectory is defined as the distance between the robots x and y position and the absolute difference in heading at each time step. The twist produced by the dynamic window approach is followed for a number of steps to ensure the robot is clear from any obstacles.

It is possible that the twist from the dynamic window approach can cause a collision in a dynamic environment. In this case the dynamic window approach will search the velocity space for the most similar collision free twist to the previous twist given by the dynamic window approach. Similarity is defined here as the inner product of the twist. If the dynamic window approach fails in this case a flag is set for the ergodic controller to replan. The ergodic controller is sensitive to the pose of the robot and is capable of optimizing a different trajectory. This strategy prevents the robot from getting stuck.

## Collision Detection
To detect obstacle cells in the occupancy grid [Bresenham's circle](http://members.chello.at/~easyfilter/bresenham.html) algorithm is used. The robot's base is modeled at a circle. If there are obstacles within a threshold of the robot's bounding radius the robot is considered to be in a collision state.


# Results


# Nodes
The nodes provided are `explore_cart` and `explore_omni`. The only difference between them
is the motion model. Both use an occupancy grid for collision detection and provide a body twist as the control output.

## Subscribed Topics
- map ([nav_msgs/OccupancyGrid](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/OccupancyGrid.html)): required for collision detection
- odom ([nav_msgs/Odometry](http://docs.ros.org/en/noetic/api/nav_msgs/html/msg/Odometry.html)): required for control feedback in the dynamic window approach

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

# Future Improvements
- Recovery behaviors for when the dynamic window approach fails or for when the robot crashes.
- Model robot's base as a polygon for collision detection.
- Use a distance field or something else other than Bresenham's circle algorithm for detecting obstacle cells. Bresenham's algorithm is fast but can alias. It is also redundant to preform Bresenham's when searching the velocity space in the dynamic window approach. The twists are propagated forward for a short time therefore the redundancy stems from ray tracing circles multiple times within the same vicinity.
- Implement the [fast shannon mutual information algorithm](https://arxiv.org/pdf/1905.02238.pdf). Provide a node capable of using this with any motion model.

# Citing Ergodic Exploration
## TODO: update version to 1.0.0 and add release tag

```
@software{ergodicexploration2020github,
  author = {Boston Cleek},
  title = {ErgodicExploration: Robot agnostic information theoretic exploration},
  url = {https://github.com/bostoncleek/ergodic_exploration},
  version = {0.0.1},
  year = {2020},
}
```

# Credits
[1] L. M. Miller and T. D. Murphey, “Trajectory optimization for continuous
ergodic exploration,” in American Control Conference, 2013, pp. 4196–4201. [PDF](https://murpheylab.github.io/pdfs/2013ACCMiMu.pdf)

[2] Decentralized ergodic control: Distribution-driven sensing and exploration for multi-agent systems
I. Abraham and T. D. Murphey
IEEE Robotics and Automation Letters, vol. 3, no. 4, pp. 2987–2994, 2018. [PDF](https://murpheylab.github.io/pdfs/2018RALAbMu.pdf), [Video](https://murpheylab.github.io/videos/2018RALAbMu.mp4)
