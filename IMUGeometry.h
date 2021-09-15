#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Geometry>
#include<vector>
#include<iostream>

struct MotionData{
    double timestamp;
    Eigen::Matrix3d Rwb;
    Eigen::Vector3d twb;

    Eigen::Vector3d imu_velocity;

    Eigen::Vector3d imu_acc;
    Eigen::Vector3d imu_gyro;

    Eigen::Vector3d imu_acc_bias;
    Eigen::Vector3d imu_gyro_bias;
};

class IMUTest{
    private:
        int imu_frequency = 200;
        int cam_frequency = 30;
        double imu_timestep = 1./imu_frequency;
        double cam_timestep = 1./cam_frequency;

        double t_start = 0.;
        double t_end = 20;

        double gyro_bias_sigma = 1.0e-5;
        double gyro_noise_sigma = 0.015;

        double acc_bias_sigma = 1.0e-4;
        double acc_noise_sigma = 0.019;

        Eigen::Matrix3d euler2Rotation(const Eigen::Vector3d& eulerAngles);
        Eigen::Matrix3d eulerRates2bodyRates(const Eigen::Vector3d& eulerAngles);

    public:
        IMUTest();
        Eigen::Vector3d gyro_bias_;
        Eigen::Vector3d acc_bias_;

        Eigen::Vector3d init_velocity_;
        Eigen::Vector3d init_Rwb_;
        Eigen::Vector3d init_twb_;

        MotionData generateIMUData(double t);
        void addIMUnoise(MotionData& motion_data);
        void testIMUInfer();
        void testIMUPreIntegration();
};
