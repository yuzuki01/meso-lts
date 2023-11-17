#ifndef HEAD_SAMPLE
#define HEAD_SAMPLE

#ifndef HEAD_CORE
#include "core.h"
#endif

#ifndef HEAD_MESH
#include "mesh.h"
#endif

#ifndef HEAD_SOLVER
#include "solver.h"
#endif


class Demo;


class Container {
public:
    std::vector<Demo> CELLS{};
    Container() = default;
};

class Demo {
public:
    int key{};
    std::vector<int> near_key{};
    std::vector<Demo *> near_demo{};
    Demo() = default;
    explicit Demo(int _key) : key(_key) {};
    void set_near(int _key, Container &_container) {
        near_key.push_back(_key);
        near_demo.push_back(&(_container.CELLS[_key]));
    }
};

int sample_vector_test();

#endif //HEAD_SAMPLE
