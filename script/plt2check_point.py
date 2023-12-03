import re

fp = "../build/result/cylinder/step_33000.dat"
D = 2

with open(fp, 'r') as file:
    file_content = file.readlines()

# 使用正则表达式匹配VARIABLES后面的内容
variables_match = re.search(r'VARIABLES\s*=\s*("[\w_-]*",?)*', file_content[0])
if variables_match:
    variables_content = variables_match.group()
    print(variables_content)
    # 提取所有变量名
    variables_list = re.findall(r'"(.*?)"', variables_content)
    print('VARIABLES变量名称:', variables_list)

# 使用正则表达式匹配结点个数N
count_match = re.search(r'ZONE\s*N=(\d+),\s*E=(\d+)', file_content[1])
if count_match:
    node_count = int(count_match.group(1))
    element_count = int(count_match.group(2))
    print('结点个数N:', node_count, '单元个数E:', element_count)


count = 0
data = {}
for i in range(len(file_content) - 2):
    line = file_content[i + 2]
    if line == "## cell value\n":
        count += 1
        data[variables_list[count + D - 1]] = list()
        continue
    elif line == "# geom\n":
        break
    if count == 0:
        continue
    # 匹配所有数字的正则表达式
    number_pattern = re.compile(r'([-+]?\d*\.?\d+[[eE]?[-+]?\d+]?)')

    # 查找所有匹配的数字
    matches = number_pattern.findall(line)
    # 打印匹配结果
    # print("匹配到的数字：", matches)
    data[variables_list[count + D - 1]] += [ float(i) for i in matches]

print(data.keys())
for key, value in data.items():
    print(key, len(value))

with open("result.check_point", 'w') as fp2:
    fp2.write("step= 0\n\ncell:\n")
    for i in range(element_count):
        fp2.write(f"{i} {data['density'][i]} {data['temperature'][i]} {data['velocity-x'][i]} {data['velocity-y'][i]} 0.0 {data['heat_flux-x'][i]} {data['heat_flux-y'][i]} 0.0\n")
    fp2.write(":end")
