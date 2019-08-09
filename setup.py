import setuptools


def readme():
    with open('README.md') as f:
        return f.read()


def requirements():
    with open('requirements.txt', "r") as fh:
        return [x for x in fh.read().split('\n') if x]


setuptools.setup(name='sleeppy',
                 version='0.13',
                 author="Yiorgos Christakis",
                 author_email="elyiorgos@gmail.com",
                 description='Python package for sleep analysis of raw accelerometer data from GeneActiv wrist watches',
                 long_description=readme(),
                 long_description_content_type="text/markdown",
                 url='http://github.com/elyiorgos/sleeppy',
                 packages=setuptools.find_packages(),
                 keywords=['sleep', 'wake', 'sensor', 'digital', 'wearable', 'python', 'rest', 'period', 'christakis'],
                 classifiers=["Programming Language :: Python :: 2.7",
                              "License :: OSI Approved :: MIT License"],
                 license='MIT',
                 zip_safe=False,
                 install_requires=requirements(),
                 include_package_data=True)
