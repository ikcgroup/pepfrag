from distutils.core import setup, Extension
import os
import sys

PACKAGE_DIR = "pepfrag"

extra_compiler_args = []
extra_link_args = []
if sys.platform != "win32":
    extra_compiler_args.append("-std=c++14")
    extra_link_args.append("-lstdc++")


cpepfrag = Extension(
    "cpepfrag",
    sources=[
        os.path.join(PACKAGE_DIR, "iongenerator.cpp"),
        os.path.join(PACKAGE_DIR, "converters.cpp"),
        os.path.join(PACKAGE_DIR, "mass.cpp"),
        os.path.join(PACKAGE_DIR, "cpepfrag.cpp"),
    ],
    language="c++11",
    extra_compile_args=extra_compiler_args,
    extra_link_args=extra_link_args
)

setup(
    name="pepfrag",
    version="0.2",
    packages=[
        "pepfrag",
    ],
    ext_modules=[
        cpepfrag,
    ]
)
