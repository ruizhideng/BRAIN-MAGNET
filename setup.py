from setuptools import find_packages, setup


setup(
    name="brain-magnet",
    version="0.1.2",
    description="BRAIN-MAGNET utilities",
    package_dir={"": "src"},
    packages=find_packages("src"),
    # Keep the base install light:
    # - `prepare_data` can run without heavy deps when users pass `--sequences-fasta`.
    # - training requires `torch` + `numpy`.
    install_requires=[],
    extras_require={
        # Training (CNN) needs torch + numpy.
        "train": ["numpy>=1.20,<3", "torch>=2.0"],
        # Genome-based extraction needs pyfaidx (optional).
        "prepare-genome": ["pyfaidx>=0.9.0.3,<0.10"],
    },
    entry_points={
        "console_scripts": [
            "brain-magnet=brain_magnet.main_cli:main",
        ]
    },
    python_requires=">=3.8",
)

