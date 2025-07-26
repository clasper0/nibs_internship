"""
1. 读取 cif 文件
2. 将 cif 文件转换为 RDKit 分子对象
3. 使用 RDKit 对分子进行分析，匹配子结构
4. 将提取到的匹配分子对进行标准化处理（参考 csd_workflow.py）
5. 将标准化后的对象保存为 sdf 文件

csd_workflow.py
1. 读取预定义 --frag1 & --frag2，并定义各自的smart, smiles, heavy_atom及其数目
2. 构建模板片段对
3. 搜索 CSD 数据库，提取匹配的分子片段对
4. 将匹配到的 hits 转换为 mol2 文件
5. 矫正氢原子
6. 对齐
7. 去除重复分子
8. 保存对齐后/未能矫正的分子为 sdf 文件
"""