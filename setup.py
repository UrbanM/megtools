from setuptools import setup


setup(
    name='megtools',
    version='0.21',
    description='A python package for working with MEG data i. e. source localization, preprocessing, vector calculations etc.',
    url='https://github.com/UrbanM/megtools',
    author='Urban Marhl',
    author_email='urban.marhl@imfm.si',
    license='GNU GPL v3',
    packages=['megtools'],
    install_requires=[
        'numpy>=1.17.4',
	'scipy>=1.3.2'
    ],
    test_suite='tests',
    classifiers=[
        'Development Status :: Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
)
