#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatXX;

struct Camera{
    Eigen::Matrix3d Rwc;
    Eigen::Quaterniond qwc;
    Eigen::Vector3d twc;
    Camera(Eigen::Matrix3d R, Eigen::Vector3d t):Rwc(R),qwc(R),twc(t) {};

    std::vector<Eigen::Vector3d> observed_points;
    Eigen::Matrix3d _K;
};


class test{
    private:
        std::vector<Camera> _cameras;
        std::vector<Eigen::Vector3d> _world_points;


     public:
        //test();
        void generateVisualData(int camNums, int worldPointsNum);
        void generateIMUData();

    public:
        void testEightPointEpipolar();
        void testTriangulation();
        void testPnP();

};
