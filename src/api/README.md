# meso/api

meso 接口

---

## 功能列表

 - [handle_help](#handle_help)
 - [handle_solver](#handle_solver)
 - [handle_parse_mesh](#handle_parse_mesh)

---

### handle_help

输出帮助信息到控制台。

---

### handle_solver

求解器的调用遵循统一格式：

```c++
int handle_solver(ConfigReader &config, ArgParser &parser)  {
    config.info();  // 配置文件信息
    SOLVER solver(config, parser);  // 实例化求解器
    solver.init();  // 求解器初始化
    if (!solver.continue_to_run) {
        warn_println("Solver init failed.");
        return -1;
    }
    solver.info();  // 求解器信息
    solver.do_save();   // 保存，第 0 步流场文件
    while (solver.continue_to_run) {
        solver.do_step();   // 推进 dt
    }
    solver.do_save();   // 保存最终结果
    return 0;
}
```

---

### handle_parse_mesh

读取网格文件，并输出 tecplot 格式的包含网格内单元体积信息的文件，用于检测网格构建的正确性。
