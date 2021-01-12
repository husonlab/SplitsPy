# setup.py
"""Setup splitspy.py

LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""

import setuptools

setuptools.setup(name='SplitsPy',
                 version='0.0.9',
                 packages=setuptools.find_packages(),
                 url='http://github.com/husonlab/SplitsPhy',
                 license='GPL',
                 author='Daniel H. Huson',
                 author_email='daniel.huson@uni-tuebingen.de',
                 description='Phylogenetic outlines',
                 keywords='phylogenetics, network',
                 entry_points={'console_scripts': ['outline=splitspy.outline:main', ], },
                 python_requires='>=3.5, <4',
                 install_requires=['numpy', 'pillow'])
