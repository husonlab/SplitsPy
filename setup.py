# setup.py
"""Setup splitspy.py

LICENSE: This is open-source software released under the terms of the
GPL (http://www.gnu.org/licenses/gpl.html).
"""


from setuptools import setup

setup(
    name='SplitsPy',
    version='0.1',
    packages=['splitspy', 'splitspy.nnet', 'splitspy.graph', 'splitspy.splits', 'splitspy.outline'],
    url='http://github.com/husonlab/SplitsPhy',
    license='GPL',
    author='Daniel H. Huson',
    author_email='daniel.huson@uni-tuebingen.de',
    description='Phylogenetic outlines'
)
