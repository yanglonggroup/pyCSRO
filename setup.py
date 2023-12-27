import setuptools

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pycsro", # Replace with your own name
    version="0.1.0",
    description="A python package for CSRO paramater calculation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yanglonggroup/pyCSRO",
    # packages=setuptools.find_packages(),
    packages=['pycsro'],
    package_dir={"": "src"},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    entry_points={'console_scripts': [
        'pycsro = pycsro.main:run_pycsro_pmsro',],
    },
    data_files = [("", ["LICENSE.txt"])],
    python_requires='>=3.9',
    zip_safe=False,
)
