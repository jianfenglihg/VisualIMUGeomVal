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


class VisualTest{
    private:
        std::vector<Camera> _cameras;
        std::vector<Eigen::Vector3d> _world_points;
        std::vector<Eigen::Vector3d> Triangulation(const Eigen::Matrix3d &R_wc, const Eigen::Vector3d &t_wc, const std::vector<Eigen::Vector3d> &left_pts, const std::vector<Eigen::Vector3d> &right_pts);
        double count(const std::vector<Eigen::Vector3d> &points, const Eigen::Matrix3d &R_wc, const Eigen::Vector3d &t_wc);

     public:
        void generateVisualData(int camNums, int worldPointsNum);
        void testEightPointEpipolar();
        void testTriangulation();
        void testPnPDLT();
        void testPnPBA();
};
