#include "tgaimage.h"




const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);


int main(int argc, char** argv) {


    TGAImage image(100, 100, TGAImage::RGB);//设定渲染窗口

    image.set(52, 41, red); //设置单点像素


    image.flip_vertically(); 
    image.write_tga_file("output.tga");
    return 0;
}