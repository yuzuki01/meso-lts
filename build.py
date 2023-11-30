import os
import shutil
from colorama import init, Fore, Back, Style


class CMakeProject:

    output_path = "./build"

    def __init__(self, root, target):
        self.target_name = target
        self.source_path = root
        self.build_path = root + "/cmake-build-minsizerel"
        if os.name == "posix":
            self.generator = "\"CodeBlocks - Unix Makefiles\""
        elif os.name == "nt":
            self.generator = "\"CodeBlocks - MinGW Makefiles\""
        else:
            self.generator = "NULL"
            raise OSError("Unsupported OS: os.name = %s" % os.name)

    def _make(self):
        # if os.path.isdir(self.build_path):
        #     shutil.rmtree(self.build_path)
        os.system(
            "cmake -DCMAKE_BUILD_TYPE=MinSizeRel -G {generator:s} -S {source_path:s} -B {build_path:s}".format(
                source_path=self.source_path, build_path=self.build_path, generator=self.generator
            )
        )

    def _build(self):
        os.system(
            "cmake --build {build_path:s} --target {target_name:s} -- -j 12".format(
                build_path=self.build_path, target_name="all"
            )
        )
        if os.name == "nt":
            for file in [it for it in os.listdir(self.output_path) if it.endswith(".dll.a")]:
                path = "{output_path:s}/{file:s}".format(output_path=self.output_path, file=file)
                print("Delete - %s" % path)
                os.remove(path)

    def _clean(self):
        os.system(
            "cmake --build {build_path:s} --target {target_name:s} -- -j 12".format(
                build_path=self.build_path, target_name="clean"
            )
        )

    def build(self):
        print("[Python] Build cmake project: %s" % self.target_name)
        self._make()
        self._clean()
        self._build()


if __name__ == '__main__':
    Core = CMakeProject('./src/core', 'core')
    Mesh = CMakeProject('./src/mesh', 'mesh')
    Solver = CMakeProject('./src/solver', 'solver')
    API = CMakeProject('./src/api', 'api')
    # 可执行文件
    Meso = CMakeProject('.', 'meso')

    Core.build()
    Mesh.build()
    Solver.build()
    API.build()
    Meso.build()
