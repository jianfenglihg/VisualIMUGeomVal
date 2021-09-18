#include "VisualGeometry.h"
#include "IMUGeometry.h"

int main(){
    VisualTest vision_obj;
    vision_obj.generateVisualData(8,100);
    vision_obj.testTriangulation();
    vision_obj.testPnPDLT();
    vision_obj.testEightPointEpipolar();

    IMUTest imu_obj;
    imu_obj.testIMUInfer();
    return 0;
}
