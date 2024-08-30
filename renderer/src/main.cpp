#include <vector>
#include <cstdlib>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Model* model = NULL;
float* shadowbuffer = NULL;

const int width = 800;
const int height = 800;

Vec3f       eye(1.2, -.8, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

TGAImage total(1024, 1024, TGAImage::GRAYSCALE);
TGAImage  occl(1024, 1024, TGAImage::GRAYSCALE);


struct ZShader : public IShader {
    mat<4, 3, float> varying_tri; //三角形顶点坐标经过一些变换之后再转换为齐次坐标之后填充的矩阵

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert)); //P*M*齐次坐标顶点 M只是从模型空间转换到摄像机视角
        varying_tri.set_col(nthvert, gl_Vertex);  //varying_tri是这个由这个三角转换之后的三个顶点坐标填充的4*3矩阵
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor& color) {
        color = TGAColor(255, 255, 255) * ((gl_FragCoord.z + 1.f) / 2.f); //计算zuffer 对应纹理color
        return false;
    }
};

struct Shader : public IShader {
    mat<2, 3, float> varying_uv; // 三角形顶点u,v坐标填充的矩阵
    mat<4, 3, float> varying_tri; //三角形顶点坐标经过一些变换之后再转换为齐次坐标之后填充的矩阵

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert)); // 获取uv坐标
        Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert)); //齐次坐标顶点
        varying_tri.set_col(nthvert, gl_Vertex); //三角形顶点
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor& color) {
        Vec2f uv = varying_uv * bar; //插值得到pixel的u、v值

        if (std::abs(shadowbuffer[int(gl_FragCoord.x + gl_FragCoord.y * width)] - gl_FragCoord.z < 1e-2)) {
            occl.set(uv.x * 1024, uv.y * 1024, TGAColor(255)); //判断是否在阴影中
        }
        color = TGAColor(255, 0, 0);
        return false;
    }
};

struct AOShader : public IShader {
    mat<2, 3, float> varying_uv;
    mat<4, 3, float> varying_tri;
    TGAImage aoimage;

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert)); //uv插值矩阵
        Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex); //三角形插值矩阵
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor& color) {
        Vec2f uv = varying_uv * bar; // uv 坐标获取
        int t = aoimage.get(uv.x * 1024, uv.y * 1024)[0]; //获取aoimage
        color = TGAColor(t, t, t);
        return false;
    }
};

Vec3f rand_point_on_unit_sphere() {
    float u = (float)rand() / (float)RAND_MAX;
    float v = (float)rand() / (float)RAND_MAX;
    float theta = 2.f * M_PI * u;
    float phi = acos(2.f * v - 1.f);
    return Vec3f(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
}
 
//计算最大仰角 采用的是HBAO
float max_elevation_angle(float* zbuffer, Vec2f p, Vec2f dir) {
    float maxangle = 0;
    for (float t = 0.; t < 1000.; t += 1.) {
        Vec2f cur = p + dir * t;
        if (cur.x >= width || cur.y >= height || cur.x < 0 || cur.y < 0) return maxangle;

        float distance = (p - cur).norm();
        if (distance < 1.f) continue;
        float elevation = zbuffer[int(cur.x) + int(cur.y) * width] - zbuffer[int(p.x) + int(p.y) * width];
        maxangle = std::max(maxangle, atanf(elevation / distance));
    }
    return maxangle;
}

int main(int argc, char** argv) {
    if (2 > argc) {
        std::cerr << "Usage: " << argv[0] << "obj/model.obj" << std::endl;
        return 1;
    }

    float* zbuffer = new float[width * height];

    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    model = new Model(argv[1]);

    TGAImage frame(width, height, TGAImage::RGB);
    lookat(eye, center, up); //相机
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4); //NDC空间
    projection(-1.f / (eye - center).norm()); //斜投影

    //第一遍pass 获取zbuffer
    ZShader zshader; 
    for (int i = 0; i < model->nfaces(); i++) {
        for (int j = 0; j < 3; j++) {
            zshader.vertex(i, j);
        }
        triangle(zshader.varying_tri, zshader, frame, zbuffer);
    }
    const int nrenders = 1;

    //屏幕空间AO 
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            if (zbuffer[x + y * width] < -1e5) continue;
            float total = 0;
            for (float a = 0; a < M_PI * 2 - 1e-4; a += M_PI / 4) {
                total += M_PI / 2 - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
            }
            total /= (M_PI / 2) * 8;
            total = pow(total, 100.f);
            frame.set(x, y, TGAColor(total * 255, total * 255, total * 255));
        }
    }

    frame.flip_vertically();
    frame.write_tga_file("framebuffer.tga");
    delete[] zbuffer;
    delete model;
    return 0;
}
