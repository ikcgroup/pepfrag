from distutils.core import setup, Extension
import os
import sys

PACKAGE_DIR = "pepfrag"

extra_compiler_args = []
if sys.platform == "darwin":
    extra_compiler_args.append("-std=c++14")


cpepfrag = Extension(
    "cpepfrag",
    sources=[
        os.path.join(PACKAGE_DIR, "iongenerator.cpp"),
        os.path.join(PACKAGE_DIR, "converters.cpp"),
        os.path.join(PACKAGE_DIR, "cpepfrag.cpp"),
    ],
    language="c++11",
    extra_compile_args=extra_compiler_args
)

setup(
    name="pepfrag",
    version="0.1",
    packages=[
        "pepfrag",
    ],
    ext_modules=[
        cpepfrag,
    ]
)
