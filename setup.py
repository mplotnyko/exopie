from setuptools import find_packages, setup

setup(
    name="exopie",
    version="1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"exopie": ["Data/*"]},
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "corner",
        "numpy",
        "pandas",
        "scipy"
    ],
)
