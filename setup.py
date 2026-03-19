from setuptools import find_packages, setup


setup(
    name="brain-magnet",
    version="0.1.2",
    description="BRAIN-MAGNET utilities",
    package_dir={"": "src"},
    packages=find_packages("src"),
    install_requires=[
        "numpy>=1.20,<3",
        "torch>=2.0",
        "pyfaidx>=0.7.0",
    ],
    entry_points={
        "console_scripts": [
            "brain-magnet=brain_magnet.main_cli:main",
        ]
    },
    python_requires=">=3.8",
)

