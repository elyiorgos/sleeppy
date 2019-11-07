import setuptools


def readme():
    with open("README.md") as f:
        return f.read()


def requirements():
    with open("requirements.txt", "r") as fh:
        return [x for x in fh.read().split("\n") if x]


def get_version():
    with open("sleeppy/version.py") as f:
        return f.readlines()[-1].split()[-1].strip("\"'")


setuptools.setup(
    name="sleeppy",
    version=get_version(),
    author="Yiorgos Christakis",
    author_email="elyiorgos@gmail.com",
    description="Python package for sleep analysis of raw accelerometer data from GeneActiv wrist watches",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="http://github.com/elyiorgos/sleeppy",
    packages=setuptools.find_packages(),
    keywords=[
        "sleep",
        "wake",
        "sensor",
        "digital",
        "wearable",
        "python",
        "rest",
        "period",
        "christakis",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
    ],
    license="MIT",
    zip_safe=False,
    install_requires=requirements(),
    include_package_data=True,
)
