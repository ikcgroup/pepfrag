import os
import sys

from setuptools import Extension, setup

PACKAGE_DIR = "pepfrag"

extra_compiler_args = []
extra_link_args = []
if sys.platform != "win32":
    extra_compiler_args.append("-std=c++14")
    extra_link_args.append("-lstdc++")
if sys.platform == "win32":
    extra_compiler_args.append("/utf-8")


# Extract README.md
readme = os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md")
with open(readme, encoding='utf-8') as fh:
    long_description = fh.read()


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
    version="0.4.1",
    packages=[
        "pepfrag",
    ],
    license="MIT",
    description="A library for peptide fragment ion generation",
    author="Daniel Spencer",
    author_email="danielspencer305@hotmail.co.uk",
    url="https://github.com/ikcgroup/pepfrag",
    keywords=[
        "PEPTIDE",
        "MASS SPECTROMETRY",
        "PROTEOMICS"
    ],
    ext_modules=[
        cpepfrag,
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8"
    ],
    long_description_content_type='text/markdown',
    long_description=long_description,
)
