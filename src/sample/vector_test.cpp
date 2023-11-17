#include "sample.h"

int sample_vector_test() {
    Container c;
    for (int i = 0; i < 10; i++) {
        c.CELLS.emplace_back(i);
    }
    for (int i = 0; i < c.CELLS.size(); i++) {
        int li = i - 1, ri = i + 1;
        auto &now = c.CELLS[i];
        if (li >= 0) now.set_near(li, c);
        if (ri < c.CELLS.size()) now.set_near(ri, c);
    }

    for (auto & now : c.CELLS) {
        now.near_demo.clear();
    }
    for (auto & now : c.CELLS) {
        std::cout << now.key << std::endl;
    }
    return 0;
}
