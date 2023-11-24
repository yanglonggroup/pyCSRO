import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="pkgname", # Replace with your own name
    version="0.0.1",
    author="example author",
    author_email="author@example.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    package_dir={"pkgname": "pkgname"},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    entry_points={'console_scripts': [
            'pkgname = pkgname.main:main',],
        },
    data_files = [("", ["LICENSE.txt"])],
    python_requires='>=3.7',
    zip_safe=False,
)
