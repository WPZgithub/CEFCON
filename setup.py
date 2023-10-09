from setuptools import setup, find_packages


def readme():
    with open('README.md', 'r', encoding='utf-8') as f:
        content = f.read()
    return content


setup(
    name='cefcon',
    version='0.2.0',
    description='Deciphering cell fate control from single-cell RNA-seq data',
    long_description=readme(),
    long_description_content_type='text/markdown',
    author='Peizhuo Wang',
    author_email='wangpeizhuo_37@163.com',
    url='https://github.com/WPZgithub/CEFCON',
    packages=find_packages('.'),
    entry_points={
        "console_scripts": ['cefcon = cefcon.CEFCON:main']
    },
    python_requires="==3.10",  # >=3.10 is recommended
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn',
        'tqdm',
        'torch>=1.13.0,<2.0',
        'torch-geometric>=2.1.0',
        'scanpy>=1.9.0',
        'networkx>=3.0',
        'cvxpy>=1.2.0',
        'gurobi>=10.0.0',
        'pyscenic>=0.12.0',
        'matplotlib',
        'matplotlib-venn',
        'seaborn',
        'biomart',
        'r',
        'rpy2',
    ],
    # package_dir={'cefcon':'cefcon'},
    license='MIT',
)
