import os
import shutil


target_name = ("cmake-build-minsizerel", "cmake-build-debug")


def find_cmake(path):
    result = []
    dir_list = [it for it in os.listdir(path) if os.path.isdir("%s/%s" % (path, it))]
    for it in dir_list:
        target = "%s/%s" % (path, it)
        if it in target_name:
            result.append(target)
        else:
            result += find_cmake(target)
    return result


if __name__ == '__main__':
    # remove
    for i in find_cmake("."):
        shutil.rmtree(i)
