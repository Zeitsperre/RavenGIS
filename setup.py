#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', ]

setup_requirements = ['pytest-runner', 'wheel', ]

test_requirements = ['pytest>=3', ]

docs_requirements = [dependency for dependency in open("requirements_docs.txt").readlines()]

dev_requirements = [dependency for dependency in open("requirements_dev.txt").readlines()]

setup(
    author="Trevor James Smith",
    author_email='smith.trevorj@ouranos.ca',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Python geotools supporting spatial operations in the RavenWPS hydrologic modelling service.",
    entry_points={
        'console_scripts': [
            'ravengis=ravengis.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords='ravengis',
    name='ravengis',
    packages=find_packages(include=['ravengis', 'ravengis.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    extras_require={
        "docs": docs_requirements,
        "dev": dev_requirements,
    },
    url='https://github.com/Zeitsperre/ravengis',
    version='0.1.0',
    zip_safe=False,
)
