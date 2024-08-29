#include <vector>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model* model = NULL;
float* shadowbuffer = NULL;

const int width = 800;
const int height = 800;

Vec3f light_dir(1, 1, 0);
Vec3f       eye(1, 1, 4);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

struct Shader : public IShader {
    mat<4, 4, float> uniform_M;   //  Projection*ModelView
    mat<4, 4, float> uniform_MIT; // (Projection*ModelView).invert_transpose()
    mat<4, 4, float> uniform_Mshadow; // transform framebuffer screen coordinates to shadowbuffer screen coordinates
    mat<2, 3, float> varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
    mat<3, 3, float> varying_tri; // triangle coordinates before Viewport transform, written by VS, read by FS

    Shader(Matrix M, Matrix MIT, Matrix MS) : uniform_M(M), uniform_MIT(MIT), uniform_Mshadow(MS), varying_uv(), varying_tri() {}

    virtual Vec4f vertex(int iface, int nthvert) {
        //vertex shader 执行了相似的功能
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = Viewport * Projection * ModelView * embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        Vec4f sb_p = uniform_Mshadow * embed<4>(varying_tri * bar); // 得到shadow buffer中对应的坐标
        sb_p = sb_p / sb_p[3];
        int idx = int(sb_p[0]) + int(sb_p[1]) * width; // index in the shadowbuffer array
        float shadow = .3 + .7 * (shadowbuffer[idx] < sb_p[2]); // 避免z-fighting
        Vec2f uv = varying_uv * bar;                 // 插值pixel uv坐标
        Vec3f n = proj<3>(uniform_MIT * embed<4>(model->normal(uv))).normalize(); // 法线
        Vec3f l = proj<3>(uniform_M * embed<4>(light_dir)).normalize(); // 光源
        Vec3f r = (n * (n * l * 2.f) - l).normalize();   // 反射向量
        float spec = pow(std::max(r.z, 0.0f), model->specular(uv)); //高光
        float diff = std::max(0.f, n * l); //漫反射
        TGAColor c = model->diffuse(uv);
        for (int i = 0; i < 3; i++) color[i] = std::min<float>(20 + c[i] * shadow * (1.2 * diff + .6 * spec), 255); 
        return false;
    }
};

struct DepthShader : public IShader {
    mat<3, 3, float> varying_tri;

    DepthShader() : varying_tri() {}

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); //转为齐次坐标

        gl_Vertex = Viewport * Projection * ModelView * gl_Vertex;        //转为屏幕空间坐标  

        varying_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));    

        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor& color) {
        Vec3f p = varying_tri * bar;            //插值得到pixel 坐标
        color = TGAColor(255, 255, 255) * (p.z / depth);
        return false; 
    }
};

int main(int argc, char** argv) {
    if (2 > argc) {
        std::cerr << "Usage: " << argv[0] << "obj/model.obj" << std::endl;
        return 1;
    }

    float* zbuffer = new float[width * height];
    shadowbuffer = new float[width * height];  //shadow map

    for (int i = width * height; --i; ) {
        zbuffer[i] = shadowbuffer[i] = -std::numeric_limits<float>::max();
    }

    model = new Model(argv[1]);
    light_dir.normalize();

    { // 光源方向渲染
        TGAImage depth(width, height, TGAImage::RGB);
        lookat(light_dir, center, up); //视口矩阵
        viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4); //观察矩阵
        projection(0); //透镜

        DepthShader depthshader;
        Vec4f screen_coords[3];
        for (int i = 0; i < model->nfaces(); i++) {
            for (int j = 0; j < 3; j++) {
                screen_coords[j] = depthshader.vertex(i, j);
            }
            triangle(screen_coords, depthshader, depth, shadowbuffer);
        }
        depth.flip_vertically(); 
        depth.write_tga_file("depth.tga");
    }

    Matrix M = Viewport * Projection * ModelView;

    { // 相机方向渲染
        TGAImage frame(width, height, TGAImage::RGB);
        lookat(eye, center, up); 
        viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
        projection(-1.f / (eye - center).norm());

        //Shader(Matrix M, Matrix MIT, Matrix MS) MIT是转移模型空间法线的矩阵，MS是将相机空间的顶点转移到灯源空间
        Shader shader(ModelView, (Projection * ModelView).invert_transpose(), M * (Viewport * Projection * ModelView).invert()); 
        Vec4f screen_coords[3];
        for (int i = 0; i < model->nfaces(); i++) {
            for (int j = 0; j < 3; j++) {
                screen_coords[j] = shader.vertex(i, j);
            }
            triangle(screen_coords, shader, frame, zbuffer);
        }
        frame.flip_vertically(); 
        frame.write_tga_file("framebuffer.tga");
    }

    delete model;
    delete[] zbuffer;
    delete[] shadowbuffer;
    return 0;
}
