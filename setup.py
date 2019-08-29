from distutils.core import setup, Extension
import os

PACKAGE_DIR = "pepfrag"

cpepfrag = Extension(
    "cpepfrag",
    sources=[
        os.path.join(PACKAGE_DIR, "iongenerator.cpp"),
        os.path.join(PACKAGE_DIR, "converters.cpp"),
        os.path.join(PACKAGE_DIR, "cpepfrag.cpp"),
    ]
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