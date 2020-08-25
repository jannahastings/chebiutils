import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chebiutils",
    version="0.2",
    author="Janna Hastings",
    author_email="janna.hastings@gmail.com",
    description="A set of functions for caching, searching and accessing the ChEBI database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jannahastings/chebiutils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

