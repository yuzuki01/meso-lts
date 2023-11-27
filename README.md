# Meso

DUGKS 求解器

重构前<a href="https://github.com/yuzuki01/meso">代码</a>已经不再维护，仅供参考

## 重大更新(Nov 14, 2023)

重构了 core 模块和 mesh 模块。

 - core 中内置了新模块
   
    - 启动参数解释器 ArgParser 对象
      
    - 日志输出 Logger 对象 (替换了原来的 pprint)
     
 - mesh 中重构了网格对象
   
    - 固定网格 ```MESH::ListMesh ``` 和可变网格 ```MESH::MapMesh```
    
    - 适用于高阶格式的次相邻网格 ```std::vector<key_type> second_near_cell_key```
    
 - 重构后求解器仅有用于验证的 dugks@incompressible

## bug

### Nov 14,2023 可变网格的输出仍有几何关系上的错误

## 未来计划

 - Lagrange 粒子追踪器 (for DPM)

 - 可变网格的网格单元增减接口 (for VoF)
