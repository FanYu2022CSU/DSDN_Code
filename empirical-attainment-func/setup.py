import setuptools


requirements = []
with open("requirements.txt", "r") as f:
    for line in f:
        requirements.append(line.strip())


setuptools.setup(
    packages=setuptools.find_packages(),
    python_requires=">=3.8",
    platforms=["Linux"],
    install_requires=requirements,
    include_package_data=True,
)
