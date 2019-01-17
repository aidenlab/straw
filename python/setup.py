import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hic-straw",
    version="0.0.2",
    author="Neva C. Durand",
    description="Extract data quickly from Juicebox hic files via straw",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aidenlab/straw",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
