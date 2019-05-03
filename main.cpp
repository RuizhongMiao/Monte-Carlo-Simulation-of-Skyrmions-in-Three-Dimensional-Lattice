#include "lattice.h"

using namespace std;

int main(int argc, const char * argv[]) {
    srand((unsigned)time(NULL));
    
    int index = 0;
    char fn[20] = "";
    
    lattice spin(30, 30, 30);
    spin.setSpin('r');
    spin.setTemperature(0.1);
    spin.setMagneticFieldUniform(0.4);
    spin.setMagneticFieldWindow(index, 10, 10, 10, 0.1);
    
    for (int i = 0; i < spin.MCS() * 10000; i++) {
        spin.oneFlip();
    }
    
    //sprintf(fn, "%d.txt",index);
    //spin.writeSpin(fn);
    spin.writeSpin("aaa.txt");
    
    /*for (index = 1; index < 20; index++) {
        spin.setMagneticFieldWindow(index, 10, 10, 10, 4);
        spin.setMagneticFieldWindow(index, 10, 10, 10, 0.1);
        
        for (int i = 0; i < spin.MCS() * 10000; i++) {
            spin.oneFlip();
        }
        
        sprintf(fn, "%d.txt",index);
        spin.writeSpin(fn);
    }*/
    
    return 0;
}
