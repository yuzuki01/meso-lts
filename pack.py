import os
import time
import zipfile
from colorama import init, Fore, Back, Style

skip_dir = {".idea": 999, ".git": 999, "cmake-build-minsizerel": 999, "cmake-build-debug": 999, "__pycache__": 999,
            "build": 0, "pack": 0}


def get_file_tree(path, depth=0):
    result = []
    file = [_x for _x in os.listdir(path) if os.path.isfile("%s/%s" % (path, _x))]
    child_dir = [_x for _x in os.listdir(path) if os.path.isdir("%s/%s" % (path, _x))]
    for _x in file:
        result.append("%s/%s" % (path, _x))
    for _x in child_dir:
        if _x in skip_dir.keys():
            if depth > skip_dir[_x]:
                result += get_file_tree("%s/%s" % (path, _x), depth + 1)
        else:
            result += get_file_tree("%s/%s" % (path, _x), depth + 1)
    return result


if __name__ == '__main__':
    init(autoreset=True)
    root = "."
    time_struct = time.localtime()
    pack_time = "{year:d}{mon:02d}{mday:02d}".format(
        year=time_struct.tm_year, mon=time_struct.tm_mon, mday=time_struct.tm_mday
    )
    zip_name = "meso-mesh-{pack_time}.zip".format(pack_time=pack_time)
    # pack
    if not os.path.isdir("./pack"):
        os.mkdir("./pack")

    with zipfile.ZipFile("./pack/%s" % zip_name, "w") as zf:
        for fp in get_file_tree(root):
            print(Fore.YELLOW + "Write file: %s" % fp + Style.RESET_ALL)
            zf.write(fp)
    print(Fore.GREEN + "ZipFile: %s/pack/%s" % (root, zip_name) + Style.RESET_ALL)
