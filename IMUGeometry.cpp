#include"IMUGeometry.h"
#include<random>
#include<fstream>

IMUTest::IMUTest(){
    gyro_bias_ = Eigen::Vector3d::Zero();
    acc_bias_  = Eigen::Vector3d::Zero();
}

Eigen::Matrix3d IMUTest::euler2Rotation(const Eigen::Vector3d& eulerAngles){
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);
    double yaw = eulerAngles(2);
    Eigen::Matrix3d R_x;
    Eigen::Matrix3d R_y;
    Eigen::Matrix3d R_z;
    Eigen::Matrix3d R;
    R_x << 1, 0, 0,
           0, cos(roll), -sin(roll),
           0, sin(roll), cos(roll);
    R_y << cos(pitch), 0, sin(pitch),
           0, 1, 0,
           -sin(pitch), 0, cos(pitch);
    R_z << cos(yaw), -sin(yaw), 0,
           sin(yaw), cos(yaw), 0,
           0, 0, 1;
    R = R_z * R_y * R_x;
    return R;
}

Eigen::Matrix3d IMUTest::eulerRates2bodyRates(const Eigen::Vector3d& eulerAngles){
    double roll = eulerAngles(0);
    double pitch = eulerAngles(1);
    Eigen::Matrix3d R;
    R << 1, 0, -sin(pitch),
         0, cos(roll), sin(roll)*cos(pitch),
         0, -sin(roll), cos(roll)*cos(pitch);
    return R;
}

MotionData IMUTest::generateIMUData(double t){
    MotionData data;
    // param
    float ellipse_x = 15;
    float ellipse_y = 20;
    float z = 1;           // z轴做sin运动
    float K1 = 10;          // z轴的正弦频率是x，y的k1倍
    float K = M_PI/ 10;    // 20 * K = 2pi 　　由于我们采取的是时间是20s, 系数K控制yaw正好旋转一圈，运动一周

    // translation
    // twb:  body frame in world frame
    Eigen::Vector3d position( ellipse_x * cos( K * t) + 5, ellipse_y * sin( K * t) + 5,  z * sin( K1 * K * t ) + 5);
    Eigen::Vector3d dp(- K * ellipse_x * sin(K*t),  K * ellipse_y * cos(K*t), z*K1*K * cos(K1 * K * t));              // position导数　in world frame
    double K2 = K*K;
    Eigen::Vector3d ddp( -K2 * ellipse_x * cos(K*t),  -K2 * ellipse_y * sin(K*t), -z*K1*K1*K2 * sin(K1 * K * t));     // position二阶导数

    // Rotation
    double k_roll = 0.1;
    double k_pitch = 0.2;
    Eigen::Vector3d eulerAngles(k_roll * cos(t) , k_pitch * sin(t) , K*t );   // roll ~ [-0.2, 0.2], pitch ~ [-0.3, 0.3], yaw ~ [0,2pi]
    Eigen::Vector3d eulerAnglesRates(-k_roll * sin(t) , k_pitch * cos(t) , K);      // euler angles 的导数

//    Eigen::Vector3d eulerAngles(0.0,0.0, K*t );   // roll ~ 0, pitch ~ 0, yaw ~ [0,2pi]
//    Eigen::Vector3d eulerAnglesRates(0.,0. , K);      // euler angles 的导数

    Eigen::Matrix3d Rwb = euler2Rotation(eulerAngles);         // body frame to world frame
    Eigen::Vector3d imu_gyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates;   //  euler rates trans to body gyro

    Eigen::Vector3d gn (0,0,-9.81);                                   //  gravity in navigation frame(ENU)   ENU (0,0,-9.81)  NED(0,0,9,81)
    Eigen::Vector3d imu_acc = Rwb.transpose() * ( ddp -  gn );  //  Rbw * Rwn * gn = gs

    data.imu_gyro = imu_gyro;
    data.imu_acc = imu_acc;
    data.Rwb = Rwb;
    data.twb = position;
    data.imu_velocity = dp;
    data.timestamp = t;
    return data;
}

void IMUTest::addIMUnoise(MotionData& data)
{
    std::random_device rd;
    std::default_random_engine generator_(rd());
    std::normal_distribution<double> noise(0.0, 1.0);

    Eigen::Vector3d noise_gyro(noise(generator_),noise(generator_),noise(generator_));
    Eigen::Matrix3d gyro_sqrt_cov = gyro_noise_sigma * Eigen::Matrix3d::Identity();
    data.imu_gyro = data.imu_gyro + gyro_sqrt_cov * noise_gyro / sqrt( imu_timestep ) + gyro_bias_;

    Eigen::Vector3d noise_acc(noise(generator_),noise(generator_),noise(generator_));
    Eigen::Matrix3d acc_sqrt_cov = acc_noise_sigma * Eigen::Matrix3d::Identity();
    data.imu_acc = data.imu_acc + acc_sqrt_cov * noise_acc / sqrt( imu_timestep ) + acc_bias_;

    // gyro_bias update
    Eigen::Vector3d noise_gyro_bias(noise(generator_),noise(generator_),noise(generator_));
    gyro_bias_ += gyro_bias_sigma * sqrt(imu_timestep ) * noise_gyro_bias;
    data.imu_gyro_bias = gyro_bias_;

    // acc_bias update
    Eigen::Vector3d noise_acc_bias(noise(generator_),noise(generator_),noise(generator_));
    acc_bias_ += acc_bias_sigma * sqrt(imu_timestep ) * noise_acc_bias;
    data.imu_acc_bias = acc_bias_;

}

void IMUTest::testIMUInfer(){
    std::cout<<"=============================================================================="<<std::endl;
    std::cout<<"============================== IMU Infer start  =============================="<<std::endl;
    std::cout<<"=============================================================================="<<std::endl;

    std::ofstream filestream_infer("./imu_infer_res.xyz");
    std::ofstream filestream_gt("./imu_gt_res.xyz");
    MotionData init_data = generateIMUData(t_start);
    Eigen::Vector3d P = init_data.twb;
    Eigen::Vector3d V = init_data.imu_velocity;
    Eigen::Quaterniond Q(init_data.Rwb);
    Eigen::Vector3d g(0,0,-9.8);
    filestream_gt << P.x() << " " << P.y() << " " << P.z() << std::endl;
    filestream_infer << P.x() << " " << P.y() << " " << P.z() << std::endl;
    for(double t=t_start+imu_timestep; t<t_end;){
        MotionData prev_data = generateIMUData(t-imu_timestep);
        MotionData curr_data = generateIMUData(t);
        filestream_gt << curr_data.twb.x() << " " << curr_data.twb.y() << " " << curr_data.twb.z() << std::endl;
        addIMUnoise(curr_data);
        addIMUnoise(prev_data);
        Eigen::Vector3d omiga = (prev_data.imu_gyro + curr_data.imu_gyro)/2;
        Eigen::Quaterniond Q_prev = Q;
        Eigen::Quaterniond dq;
        dq.w() = 1;
        dq.x() = 0.5*omiga.x()*imu_timestep;
        dq.y() = 0.5*omiga.y()*imu_timestep;
        dq.z() = 0.5*omiga.z()*imu_timestep;
        Q = Q * dq;
        /* Eigen::Vector3d acc = 0.5*(Q_prev*(prev_data.imu_acc+prev_data.imu_acc_bias)+Q*(curr_data.imu_acc+curr_data.imu_acc_bias))+g; */
        Eigen::Vector3d acc = 0.5*(Q_prev*prev_data.imu_acc+Q*curr_data.imu_acc)+g;
        V = V + acc*imu_timestep;
        P = P + V*imu_timestep + 0.5*acc*imu_timestep*imu_timestep;
        filestream_infer << P.x() << " " << P.y() << " " << P.z() << std::endl;
        t += imu_timestep;
    }
    filestream_gt.close();
    filestream_infer.close();

    std::cout<<"the result saved in imu_infer_res and imu_gt_res!"<<std::endl;
    std::cout<<"=============================================================================="<<std::endl;
    std::cout<<"============================== IMU Infer end ================================="<<std::endl;
    std::cout<<"=============================================================================="<<std::endl;

}

void IMUTest::testIMUPreIntegration(){
    std::cout<<"=============================================================================="<<std::endl;
    std::cout<<"========================= IMU PreIntergration start =========================="<<std::endl;
    std::cout<<"=============================================================================="<<std::endl;



    std::cout<<"=============================================================================="<<std::endl;
    std::cout<<"========================= IMU PreIntergration end ============================"<<std::endl;
    std::cout<<"=============================================================================="<<std::endl;
}
