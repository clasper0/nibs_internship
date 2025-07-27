# coding=utf-8
from pymatgen.io.cif import CifParser   # 更完整地读取 cif 文件信息（包括键长等）
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core.periodic_table import Element    # 获取元素半径

parser = CifParser('/home/clasper/2507/cif/1/00/40/1004002.cif')
structure = parser.parse_structures(primitive=True)[0]
cif_data = parser.as_dict()['1004002']

# 获取成键原子信息
label1 = cif_data['_geom_bond_atom_site_label_1']
label2 = cif_data['_geom_bond_atom_site_label_2']
labels = list(zip(label1, label2))
# 获取所有原子的相邻成键原子（但不包括 H ）
bond_dict = {}
for label1, label2 in labels:
    if label1 not in bond_dict:
        bond_dict[label1] = []
    bond_dict[label1].append(label2)
    if label2 not in bond_dict:
        bond_dict[label2] = []
    bond_dict[label2].append(label1)
# print(bond_dict)

# 建立 {label: index} 字典
label_to_index = {site.label: i for i, site in enumerate(structure)}
print(label_to_index)   # Warning: 一个label可能有多个index，这是因为等效的 index 会被分配相同的 label

# analyzer = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(structure)
vnn = VoronoiNN()
for i, site in enumerate(structure):
    site_dict = site.as_dict()
    if site_dict['species'][0]['occu'] != 1:
        continue
    if not site.is_ordered:
        print(site) # 出现异常未排除 site，打印出来
    score = 0
    patt = 'CCCCCC'
    # 通过对 score 和 patt 的双重约束保证子结构的匹配
    # dict['properties'] = []
    # print(dict)
    route = set()
    site_element = site_dict['species'][0]['element']
    site_label = site_dict['label']
    # if site_element == 'C':    # 待改进，从 patt 中获取路线信息
    #     if site_dict['label'] not in route:  # 如果不走回头路，则加到 route 中，可继续前进
    #         route.add(site_label)
    #         neighbors = [neighbor for neighbor in vnn.get_nn_info(structure, i) if neighbor['image'] == (0, 0, 0)]   # 取原胞内的邻居，不一定成键
    #         # print(dict)
    #         bond_neighbors = []
    #         for neighbor in neighbors:
    #             # print(neighbor['site'], site.distance(structure[neighbor['site_index']]))   # 打印邻居和距离
    #             neighbor_dict = structure[neighbor['site_index']].as_dict()
    #             neighbor_element = neighbor_dict['species'][0]['element']
    #             neighbor_label = neighbor_dict['label']
    #             print("Before: ", site_label, neighbor_label, site.distance(structure[neighbor['site_index']]))
    #             if (site_label, neighbor_label) in labels or (neighbor_label, site_label) in labels:
    #                 bond_neighbors.append(neighbor_dict)
    #                 # print(site_label, neighbor_label, site.distance(structure[neighbor['site_index']]))
    #                 print("After: ", site_label, neighbor_label, site.distance(structure[neighbor['site_index']]), site_element, neighbor_element)
    #         print(len(bond_neighbors))
    #         break
    # raw_neighbors = vnn.get_nn_info(structure, i)
    # neighbors = [neighbor for neighbor in raw_neighbors if neighbor['image'] == (0, 0, 0)]  # 只保留原胞内的邻居
    # for neighbor in neighbors:
    #     print(site, neighbor)
    # break