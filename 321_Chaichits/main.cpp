#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include <omp.h>
#include <cstring>
#include <unordered_map>
#include "geom.h"
#include "Bitmap.h"

typedef vec<3, float> Vec3f;
typedef vec<3, int> Vec3i;
typedef vec<4, float> Vec4f;

int scene;

struct Material {
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
    bool is_ch;
    Material(const float r, const Vec4f &a, const Vec3f &color, const float spec, bool ch=false) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec), is_ch(ch) {}
    Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L * dir;
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrt(radius * radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

const int DEPTH = 4; // depth of recursion for tree
class Model {
    std::vector<Vec3f> verts;
    std::vector<Vec3i> faces;
    std::vector<std::pair<Vec3f, Vec3f>> paral; // vector of aabb boxes
    std::vector<std::vector<int>> faces_in_bbox; // vector of triangles in each aabb box

public:
	Vec3f n1; // min of aabb box of whole model
    Vec3f n2; // max of aabb box of whole model

    // Loading model from .obj file
    Model(const char *filename) : verts(), faces(), paral(20000000), faces_in_bbox(20000000), n1(Vec3f(0, 0, 0)), n2(Vec3f(0, 0, 0)) {
        std::ifstream in;
        in.open (filename, std::ifstream::in);
        if (in.fail()) {
            std::cerr << "Failed to open " << filename << std::endl;
            return;
        }
        std::string line;
        while (!in.eof()) {
            std::getline(in, line);
            std::istringstream iss(line.c_str());
            char trash;
            if (!line.compare(0, 2, "v ")) {
                iss >> trash;
                Vec3f v;
                for (int i = 0; i < 3; i++) {
                    iss >> v[i];
                    v[i] /= 15; // scaling model
                }
                // making right position
            	v[0] += 7;
           		v[1] += 1.;
            	v[2] -= 13.;
                verts.push_back(v);
            } else if (!line.compare(0, 2, "f ")) {
                Vec3i f;
                int idx, cnt=0;
                iss >> trash;
                while (iss >> idx) {
                    idx--;
                    f[cnt++] = idx;
                }
                if (3 ==  cnt) faces.push_back(f);
            }
        }
        get_bbox(n1, n2); 
        paral[1] = std::make_pair(n1, n2);
        get_bbox_map(n1, n2); // making tree of model
    }

    int nverts() const { // number of verticies
        return (int)verts.size(); 
    }

    int nfaces() const { // number of faces
        return (int)faces.size();
    }

    // fast minimal storage ray triangle intersection
    bool ray_triangle_intersect(const int &fi, const Vec3f &orig, const Vec3f &dir, float &tnear, Vec3f &norm) const {
        Vec3f edge1 = point(vert(fi,1)) - point(vert(fi,0));
        Vec3f edge2 = point(vert(fi,2)) - point(vert(fi,0));
        Vec3f pvec = cross(dir, edge2);
        float det = edge1 * pvec;
        norm = cross(edge1, edge2).normalize();
        if (norm * dir > 0.) { // if normal if inverted inverts it back
        	norm = -norm;
        }
        if (fabs(det) < 1e-5) return false;

        Vec3f tvec = orig - point(vert(fi,0));
        float u = tvec * pvec / det;
        if (u < 0.0 || u > 1.) return false;

        Vec3f qvec = cross(tvec, edge1);
        float v = dir * qvec / det;
        if (v < 0.0 || u + v > 1.) return false;

        tnear = edge2 * qvec * (1. / det);
        return tnear > 1e-5;
    }

    const Vec3f &point(int i) const { // coordinates of the vertex i
        assert(i >= 0 && i < nverts());
        return verts[i];
    }

    Vec3f &point(int i) { // coordinates of the vertex i
        assert(i >= 0 && i < nverts());
        return verts[i];
    }

    int vert(int fi, int li) const { // index of the vertex for the triangle fi and local index li
        assert(fi >= 0 && fi < nfaces() && li >= 0 && li < 3);
        return faces[fi][li];
    }

    void get_bbox(Vec3f &min, Vec3f &max) { // bounding box for all the vertices
        min = max = verts[0];
        for (int i = 1; i < nverts(); i++) {
            for (int j = 0; j < 3; j++) {
                min[j] = std::min(min[j], verts[i][j]);
                max[j] = std::max(max[j], verts[i][j]);
            }
        }
    }

    bool is_inside(const Vec3f &min, const Vec3f &max, const int ind) const { // checks if face i is inside the aabb box
    	Vec3f min_tr = point(vert(ind, 0)), max_tr = point(vert(ind, 0));
    	for (int i = 1; i < 3; i++) {
    		Vec3f v = point(vert(ind, i));
    		for (int j = 0; j < 3; j++) { // making bounding box for the face
    			min_tr[j] = std::min(min_tr[j], v[j]);
    			max_tr[j] = std::max(max_tr[j], v[j]);
    		}
    	}
    	for (int i = 0; i < 3; ++i) { // checking if boxes intersect
			if (min[i] > max_tr[i] || max[i] < min_tr[i]) {
				return false;
			}
    	}
		return true;
    }

    void get_bbox_map(const Vec3f &min, const Vec3f &max, int depth=DEPTH, int ind=1) { // recursively build tree 
    	if (depth == 0) {
   			for (int j = 0; j < nfaces(); j++) {
   				if (is_inside(min, max, j)) {
   					faces_in_bbox[ind].push_back(j);
   				}
   			}
		} else if (depth > 0) {
    		Vec3f mid = (min + max) * 0.5;
	    	Vec3f delta = mid - min;
	    	for (size_t i = 0; i < 8; i++) {
	    		Vec3f add((float)((i >> 2) & 1) * delta[0], (float)((i >> 1) & 1) * delta[1], (float)(i & 1) * delta[2]);
	    		assert(ind * 8 + i < paral.size());
	   			paral[ind * 8 + i] = std::make_pair(min + add, mid + add);
	   			get_bbox_map(min + add, mid + add, depth - 1, ind * 8 + i);
	    	}
		}
    }


    // using built tree find the nearest face for each box
    void find_intersection(const Vec3f &orig, const Vec3f &dir, std::vector<int> &s,int depth=DEPTH, int ind=1) const {
    	if (depth > 0) {
	    	for (size_t i = 0; i < 8; i++) {
	    		if (in_bbox(orig, dir, paral[ind * 8 + i].first, paral[ind * 8 + i].second)) {
	    			find_intersection(orig, dir, s, depth - 1, ind * 8 + i);
	    		}
	    	}
    	} else {
    		Vec3f norm_i;
    		float dist_i = 0, dist = std::numeric_limits<float>::max();
    		int q = -1;
    		for (size_t i = 0; i < faces_in_bbox[ind].size(); i++) {
    			if (ray_triangle_intersect(faces_in_bbox[ind][i], orig, dir, dist_i, norm_i) && dist_i < dist) {
    				dist = dist_i;
    				q = faces_in_bbox[ind][i];
    			}
    		}
    		if (q != -1)
    			s.push_back(q);
    	}
    	return;
    }

    // checks whether ray instersects box
    bool in_bbox(const Vec3f &orig, const Vec3f &dir, const Vec3f &min, const Vec3f &max) const {
    	float low = (min[0] - orig[0]) / dir[0];
		float high = (max[0] - orig[0]) / dir[0];
		float tmin = std::min(low, high), tmax = std::max(low, high);
		
		for (int i = 1; i < 3; i++) {
			low = (min[i] - orig[i]) / dir[i];
			high = (max[i] - orig[i]) / dir[i];
			tmin  = std::max(tmin, std::min(low, high));
			tmax = std::min(tmax, std::max(low, high));
		}
		return (tmin <= tmax) && (tmax > 0.f);
	}
};

struct Light {
    Vec3f position;
    float intensity;
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) { // Fresnel's law
    return I - N * 2.f * (I * N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I * N));
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if (cosi < 0) { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
        cosi = -cosi;
        std::swap(etai, etat);
        n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(0, 0, 0) : I * eta + n * (eta * cosi - sqrtf(k));
}

std::vector<Vec3f> t;

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Model> &models, Vec3f &hit, Vec3f &N, Material &material) {
    // model intersection
    float models_dist = std::numeric_limits<float>::max();
    for (size_t j = 0; j < models.size(); j++) {
        if (models[j].in_bbox(orig, dir, models[j].n1, models[j].n2)) {
        	std::vector<int> s; // vector of nearest faces for each box
        	models[j].find_intersection(orig, dir, s);
    		float dist_i = 0;
        	Vec3f norm;
        	for (auto &i : s) {
        		if (models[j].ray_triangle_intersect(i, orig, dir, dist_i, norm) && dist_i < models_dist) {
    				models_dist = dist_i;
    				N = norm;
    				hit = orig + dir * models_dist;
    				material.refractive_index = 1.0;
		            material.diffuse_color = Vec3f(0.5, 0.5, 0.5);
		            material.diffuse_color = material.diffuse_color * .3;
		            material.specular_exponent = 10.;
		            material.albedo = Vec4f(0.9,  0.1, 0.0, 0.0);
    			}
        	}        	
        	
        }
    }

    // spheres insersection
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist && dist_i < models_dist) {
            spheres_dist = dist_i;
            hit = orig + dir * dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
            if (material.is_ch){ // making checkerball
                material.diffuse_color = (int(1.5 * hit[1] + 100) + int(24.0 * atan((hit[0] -  spheres[i].center[0]) / (hit[2] - spheres[i].center[2])) / M_PI + 100)) & 1 ? Vec3f(1, 1, 1) : Vec3f(0, 0, 0);
            }
        } 
    }


    // checkerboard y = -4
    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir[1]) > 1e-3 && scene != 3)  {
        float d = -(orig[1] + 4)/dir[1]; // distance to the board
        Vec3f pt = orig + dir * d; // hit point
        if (d > 0 && d < spheres_dist && d < models_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0, 1, 0);
            material = Material();
            material.diffuse_color = (int(.5 * hit[0] + 1000) + int(.5 * hit[2])) & 1 ? Vec3f(1, 1, 1) : Vec3f(0, 0, 0);
            material.diffuse_color = material.diffuse_color *.3;
            material.specular_exponent = 50.;
            material.albedo = Vec4f(0.6,  0.3, 0.1, 0.0);
        }
    }

    // board x = -15
    float board_dist = std::numeric_limits<float>::max();
    if (fabs(dir[0]) > 1e-3 && scene != 3)  {
        float d = -(orig[0] + 15) / dir[0];
        Vec3f pt = orig + dir * d;
        if (d > 0 && pt[1] > -4. && pt[1] < 14. && pt[2] < -10. && pt[2] > -42. && d < checkerboard_dist && d < spheres_dist && d < models_dist) {
            board_dist = d;
            hit = pt;
            N = Vec3f(1, 0, 0);
            material = Material();
            material.diffuse_color = t[(1080 - (int)((pt[1] + 4.) * 60.)) * 1920 + (int)((pt[2] + 42.) * 60.)];
            material.specular_exponent = 50.;
            material.albedo = Vec4f(0.6,  0.3, 0.1, 0.0);
        }
    }
    return std::min(spheres_dist, std::min(checkerboard_dist, std::min(models_dist, board_dist))) < 1000; // find the nearest object
}

std::vector<Vec3f> env_room; // 3d envmap vector
int e_w, e_h;

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const std::vector<Model> &models, const Vec3f &back, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth > 3 || !scene_intersect(orig, dir, spheres, models, point, N, material)) { // background
        if (scene == 3) { // 3d envmap
			Sphere env(Vec3f(0, 0, 0), 1000, Material());
			float dist = 0;
			env.ray_intersect(orig, dir, dist);
			Vec3f hit = orig + dir * dist;
			size_t a = (atan2(hit[2], hit[0]) / (2 * M_PI) + .5) * e_w;
			size_t b = acos(hit[1] / 1000) / M_PI * e_h;
			if (a + b * e_w < env_room.size())
				return env_room[a + b * e_w];
			else
				return Vec3f(0.9, 0.0, 0.8);
		} else { // 2d envmap
        	return back;
        }
    }

    Vec3f reflect_dir = reflect(dir, N).normalize(); // direction of reflected ray
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize(); // direction of refracted ray
    Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, models, back, depth + 1); // find color in reflection point
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, models, back, depth + 1); // find color in refraction point
    
    float diffuse_light_intensity = 0, specular_light_intensity = 0;

    for (size_t i = 0; i < lights.size(); i++) { // making light and shadow places
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, models, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;
        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N) * (1. / (lights[i].position - point).norm()) * (1. / (lights[i].position - point).norm());
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.specular_exponent) * lights[i].intensity * (1. / (lights[i].position - point).norm()) * (1. / (lights[i].position - point).norm());
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] + reflect_color * material.albedo[2] + refract_color * material.albedo[3];
}


void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights, const std::vector<Model> &models, std::string name, int threads) {
    const int width = 1920;
    const int height = 1080;
    const float fov = M_PI / 3.;
    
    // 2d envmap
    std::ifstream ifs;
    ifs.open("../env.ppm");
    int w, h;
    std::string tmp;
    
    ifs >> tmp >> w >> h >> tmp;
    std::vector<Vec3f> inp(w*h);

    // loading .ppm file
    unsigned char c;
    ifs >> std::noskipws >> c; 
    for (int i = 0; i < h * w; ++i) {
        std::vector<float> vc(3);
        for (size_t j = 0; j < 3; j++) {
            ifs >> std::noskipws >> c;
            vc[j] = 1.f / 255 * (float)c;
        }
        inp[i] = Vec3f(vc[0], vc[1], vc[2]);
    }
            
    std::vector<Vec3f> framebuffer(width * height);

    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (size_t j = 0; j < height; j++) {
        for (size_t i = 0; i < width; i++) {
            float x = (2 * (i + 0.5) / (float)width  - 1) * tan(fov / 2.) * width / (float)height;
            float y = -1 * (2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            int ind = 2 * j <= height ? j : height - j - 1;
            framebuffer[(width - i - 1) + j * width] = cast_ray(Vec3f(0, 0, 7), dir, spheres, lights, models, inp[i + ind * width]);
        }
    }

    // writing result
    std::vector<uint32_t> image(width * height);
    for(size_t i = 0; i < width * height; i++) {
    	image[i] = 0;
    	for (int j = 2; j >= 0; j--) {
            image[i] |= (uint32_t)(255 * std::max(0.f, std::min(1.f, framebuffer[width * height - i - 1][j])));
            image[i] <<= 8;
        }
        image[i] >>= 8;
    }
    
    SaveBMP(name.c_str(), image.data(), width, height);
}


Vec3f RGB(int r, int g, int b) {
    return Vec3f((float)r / 255., (float)g / 255., (float)b / 255.);
}

int main(int argc, char** argv) {
	std::string name = "zout.bmp";
	int threads = 1;
	std::unordered_map<std::string, std::string> cmdLineParams;

	for (int i = 0; i < argc; i++) {
		std::string key(argv[i]);

		if (key.size() > 0 && key[0] == '-') {
			if (i != argc - 1) {// not last argument
				cmdLineParams[key] = argv[i + 1];
				i++;
			}
			else
				cmdLineParams[key] = "";
		}
	}

	if(cmdLineParams.find("-out") != cmdLineParams.end())
		name = cmdLineParams["-out"];

	if(cmdLineParams.find("-scene") != cmdLineParams.end())
   		scene = atoi(cmdLineParams["-scene"].c_str());

   	if(cmdLineParams.find("-threads") != cmdLineParams.end())
   		threads = atoi(cmdLineParams["-threads"].c_str());

	Material           red(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB(102,  20,  20), 50);
    Material     aluminium(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), RGB( 72,  75,  87), 10.);
    Material        mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), RGB(255, 255, 255), 1425.);
    Material         glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), RGB(153, 179, 204), 125.);
    Material purple_baloon(1.0, Vec4f(0.2,  0.5, 0.1, 0.8), RGB(130,  20, 150), 125.);
    Material       checker(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB(102,  20,  20), 50, true);
    Material         green(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB( 20, 102,  20), 50);
    Material          blue(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB( 20,  20, 102), 50);
    Material        orange(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB(150, 102,  20), 50);
    Material        purple(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), RGB(102,  20, 102), 50);

    std::vector<Sphere> spheres;
    std::vector<Model> models;
    std::vector<Light> lights;
	


	switch (scene) {
		case 1: {
			spheres.push_back(Sphere(Vec3f(  0,    0, -16),   2, red));
			spheres.push_back(Sphere(Vec3f( -1, -1.5, -12),   2, glass));
			spheres.push_back(Sphere(Vec3f( 17, -1.5, -18),   3, aluminium));
			spheres.push_back(Sphere(Vec3f(  7,   -1, -10),   3, checker));
			spheres.push_back(Sphere(Vec3f(-10,   -4, -10),   4, mirror));
			spheres.push_back(Sphere(Vec3f(-10,    5, -22),   5, purple_baloon));

			spheres.push_back(Sphere(Vec3f( -3, -3.5,  -7), 0.5, aluminium));
		    spheres.push_back(Sphere(Vec3f(  5, -3.5,  -3), 0.5, purple));
		    spheres.push_back(Sphere(Vec3f(  2, -3.5, -10), 0.5, glass));
		    spheres.push_back(Sphere(Vec3f(  0, -3.5,  -9), 0.5, green));
		    spheres.push_back(Sphere(Vec3f( -1, -3.5,  -3), 0.5, blue));

		    spheres.push_back(Sphere(Vec3f( -5, -3.5,  -4), 0.5, mirror));
		    spheres.push_back(Sphere(Vec3f(  3, -3.5,  -7), 0.5, red));
		    spheres.push_back(Sphere(Vec3f( -3, -3.5,  -2), 0.5, glass));
		    spheres.push_back(Sphere(Vec3f(  3, -3.5,  -2), 0.5, green));
		    spheres.push_back(Sphere(Vec3f(  2, -3.5,  -5), 0.5, orange));
		    spheres.push_back(Sphere(Vec3f( -8, -3.5,  -4), 0.5, aluminium));
			
		    lights.push_back(Light(Vec3f(-20, 20,  20), 5000));
    		lights.push_back(Light(Vec3f( 30, 50, -10), 7000));
    		lights.push_back(Light(Vec3f( 30, 20,  30), 6000));
			break;
		}
		case 2: {
		    spheres.push_back(Sphere(Vec3f(-10, -1, -10), 3, glass));
		    spheres.push_back(Sphere(Vec3f(10, -1, -5), 2.5, checker));
		    spheres.push_back(Sphere(Vec3f(1.0, -1.5, -10), 2, mirror));
		    spheres.push_back(Sphere(Vec3f(-10, 5, -22), 5, purple_baloon));

			models.push_back((Model("../pf.obj"))); // 2938 faces

			lights.push_back(Light(Vec3f(-20, 20,  20), 5000));
		    lights.push_back(Light(Vec3f( 30, 50, -10), 7000));
		    lights.push_back(Light(Vec3f( 30, 20,  30), 6000));
			break;
		}
		case 3: {
			// 3d envmap
			std::ifstream ifs;
		    ifs.open("../envroom.ppm");
		    std::string tmp;

		    ifs >> tmp >> e_w >> e_h >> tmp;
		    env_room.resize(e_w * e_h);

		    unsigned char c;
		    ifs >> std::noskipws >> c;
		    for (int i = 0; i < e_h * e_w; ++i) {
		        std::vector<float> vc(3);
		        for (size_t j = 0; j < 3; j++) {
		            ifs >> std::noskipws >> c;
		            vc[j] = 1.f / 255 * (float)c;
		        }
		        env_room[i] = Vec3f(vc[0], vc[1], vc[2]);
		    }

		    spheres.push_back(Sphere(Vec3f(  0, 1, -15), 5, mirror));
		    spheres.push_back(Sphere(Vec3f(-12, 1, -15), 5, glass));
		    spheres.push_back(Sphere(Vec3f( 12, 1, -15), 5, blue));

		    lights.push_back(Light(Vec3f( 80, 60, 30), 10000));
		    lights.push_back(Light(Vec3f(-10, 60, 30), 10000));
		    lights.push_back(Light(Vec3f(-90, 60, 30), 10000));
			break;
		}
		default:
			return 0;
	}
    
    // picture for the second checkerboard
    std::ifstream ifs;
    ifs.open("../planet.ppm");
    int w, h;
    std::string tmp;

    ifs >> tmp >> w >> h >> tmp;
    t.resize(w * h);

    unsigned char c;
    ifs >> std::noskipws >> c;
    for (int i = 0; i < h * w; ++i) {
        std::vector<float> vc(3);
        for (size_t j = 0; j < 3; j++) {
            ifs >> std::noskipws >> c;
            vc[j] = 1.f / 255 * (float)c;
        }
        t[i] = Vec3f(vc[0], vc[1], vc[2]);
    }

    render(spheres, lights, models, name, threads);
    return 0;
}
