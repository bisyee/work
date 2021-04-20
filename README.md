# work

Introduction
We use an Extended Kalman Filter (EKF) that fuses data from
an inertial measurement unit (IMU) along with wheel encoder data from the robot.
The work is replacing the
GPS measurement with wheel encoding. We use accelerometer and gyroscope
data from the IMU. We approximate the  model as a differential drive
robot.

Notation
In the algorithm, we estimate the position, velocity, orientation, gyroscope bias,
and accelerometer bias. We denote our state with x:
x = [p, v, b, g, a ]
(1)
In the above, p is position, v is velocity, b is the quaternion describing inertial to
body frame transformation, a is accelerometer bias, and  g is gyroscope bias.

Extended Kalman Filter

We use the extended Kalman filter as follows:

Gyroscope and accelerometer are giving information for the Vehicle Kinematics. From the Vehicle Kinematics we predict the x with Kalman Filter and correct it. To check if the prediction model we use the pose and attitude measurement which are measured by the wheel encoders and accelerometer

