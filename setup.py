from setuptools import setup, find_packages

setup(
    name="nulloretriever",
    version="0.1.0",
    description="Base package for the needed functions of the nullo-retriever pipeline.",
    url="https://github.com/Masthetheus/nulloretriever",
    author="Matheus Pedron Cassol",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Researchers",
        "Topic :: Bioinformatics :: Nullomer :: Pipelines",
        "License :: GNU License",
    ],
    keywords="bioinformatics, pipeline, nullomer",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "bitarray",
        "numpy",
        "pandas",
    ]
)