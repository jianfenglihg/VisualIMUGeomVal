#include "VisualGeometry.h"


int main(){
    VisualTest vision_obj;
    vision_obj.generateVisualData(8,100);
    vision_obj.testTriangulation();
    vision_obj.testPnP();
    vision_obj.testEightPointEpipolar();
    return 0;
}
